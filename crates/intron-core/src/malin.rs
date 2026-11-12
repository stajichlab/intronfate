//! Reading and writing MALIN intron tables.
//!
//! MALIN (Csűrös 2008) takes a table with one row per taxon and one character per intron *site*:
//!
//! ```text
//! #ORGANISM Hsap  Homo sapiens
//! Hsap    003000030...            0 no intron · ? not observable · 3 / 1 / 2 intron of phase 0 / 1 / 2
//! #MAP    site    kog     alipos  # one line per site
//! #MAP    0       KOG0004 1427
//! ```
//!
//! `alipos` is a nucleotide position in the protein alignment `kog`: `3k` lies between codons `k` and
//! `k+1`, and `3k+1`, `3k+2` lie inside codon `k+1`. That is exactly this project's
//! `(column, phase) = (alipos / 3, alipos % 3)`, so a MALIN site and a [`Site`] are the same thing and
//! conversion is exact in both directions. The phase digit repeats `alipos % 3` (the MALIN guide calls a
//! codon-boundary intron "phase 3"), and a table whose digit contradicts its position is rejected.
//!
//! Tables without `#MAP` lines (for example MALIN's `rogozin_introns_table`) hold `0/1` only: they carry
//! no phase and no alignment, so they import as one matrix with [`PhaseMode::Ignore`] and the site
//! index as its column. Such a matrix cannot be exported again, because its phase is unknown.

use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::matrix::{Bitset, Cell, Site, SiteMatrix};
use crate::project::PhaseMode;

/// MALIN character for a present intron of the given phase.
fn present_code(phase: u8) -> char {
    match phase {
        0 => '3',
        1 => '1',
        _ => '2',
    }
}

/// Phase encoded by a MALIN presence digit.
fn phase_of_code(c: char) -> Option<u8> {
    match c {
        '3' => Some(0),
        '1' => Some(1),
        '2' => Some(2),
        _ => None,
    }
}

/// Write matrices (all over the same taxa, in the same order) as one MALIN table with `#MAP` lines.
///
/// Errors when the matrices disagree about their taxa, or when one lacks phase information
/// ([`PhaseMode::Ignore`]), since MALIN positions are `3 * column + phase`.
pub fn to_malin(matrices: &[&SiteMatrix]) -> Result<String> {
    let first = matrices
        .first()
        .ok_or_else(|| Error::parse("no matrix to export"))?;
    let mut seen_names = std::collections::HashSet::new();
    for m in matrices {
        if !seen_names.insert(m.alignment.replace(char::is_whitespace, "_")) {
            return Err(Error::parse(format!(
                "two matrices are both named {:?}; a MALIN `#MAP` identifies sites by alignment name, so \
                 re-importing would merge them",
                m.alignment
            )));
        }
        if m.taxa != first.taxa {
            return Err(Error::parse(format!(
                "alignments {:?} and {:?} have different taxa; a MALIN table needs one shared taxon list",
                first.alignment, m.alignment
            )));
        }
        if matches!(m.phase_mode, PhaseMode::Ignore) {
            return Err(Error::parse(format!(
                "alignment {:?} was built with --ignore-phase, so its sites have no phase and cannot be \
                 written as MALIN positions",
                m.alignment
            )));
        }
    }
    if matrices.iter().all(|m| m.sites.is_empty()) {
        return Err(Error::parse(
            "no sites to export: a MALIN table needs at least one intron site",
        ));
    }
    let mut rows: Vec<String> = vec![String::new(); first.taxa.len()];
    let mut map = String::new();
    let mut site_no = 0usize;
    for m in matrices {
        for site in &m.sites {
            for (t, row) in rows.iter_mut().enumerate() {
                row.push(match site.cell(t) {
                    Cell::Present => present_code(site.phase),
                    Cell::Absent => '0',
                    Cell::Missing => '?',
                });
            }
            let alipos = site
                .col
                .checked_mul(3)
                .and_then(|v| v.checked_add(u32::from(site.phase)))
                .ok_or_else(|| {
                    Error::parse(format!(
                        "column {} of {:?} is too large to express as a MALIN position",
                        site.col, m.alignment
                    ))
                })?;
            map.push_str(&format!(
                "#MAP\t{site_no}\t{}\t{alipos}\n",
                m.alignment.replace(char::is_whitespace, "_"),
            ));
            site_no += 1;
        }
    }
    let mut out = String::new();
    for name in &first.taxa {
        out.push_str(&format!("#ORGANISM {name}\t{name}\n"));
    }
    for (name, row) in first.taxa.iter().zip(&rows) {
        out.push_str(&format!("{name}\t{row}\n"));
    }
    out.push_str(
        "#MAP\tsite\tkog\talipos\t# intron site (0-based) -> alignment, nucleotide position\n",
    );
    out.push_str(&map);
    Ok(out)
}

/// Parse a MALIN table into one matrix per alignment (or a single binary matrix without `#MAP`).
///
/// `fallback_name` names the matrix of a table that has no `#MAP` block.
pub fn from_malin(text: &str, fallback_name: &str) -> Result<Vec<SiteMatrix>> {
    let mut taxa: Vec<String> = Vec::new();
    let mut rows: Vec<Vec<char>> = Vec::new();
    let mut map: BTreeMap<usize, (String, u32)> = BTreeMap::new();
    for (ln, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if let Some(rest) = line.strip_prefix("#MAP") {
            let f: Vec<&str> = rest.split_whitespace().collect();
            if let (Some(site), Some(aln), Some(pos)) = (
                f.first().and_then(|s| s.parse::<usize>().ok()),
                f.get(1),
                f.get(2).and_then(|s| s.parse::<u32>().ok()),
            ) {
                map.insert(site, (aln.to_string(), pos));
            }
            continue;
        }
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let (Some(name), Some(chars)) = (parts.next(), parts.next()) else {
            return Err(Error::parse(format!(
                "line {}: expected `<taxon> <characters>`",
                ln + 1
            )));
        };
        taxa.push(name.to_string());
        rows.push(chars.chars().collect());
    }
    if taxa.is_empty() {
        return Err(Error::parse("MALIN table has no taxon rows"));
    }
    if let Some(dup) = taxa
        .iter()
        .enumerate()
        .find(|(i, t)| taxa[..*i].contains(t))
        .map(|(_, t)| t)
    {
        return Err(Error::parse(format!(
            "taxon {dup:?} appears twice in the MALIN table"
        )));
    }
    let n_sites = rows[0].len();
    if rows.iter().any(|r| r.len() != n_sites) {
        return Err(Error::parse("MALIN table rows differ in length"));
    }
    if let Some(bad) = rows.iter().flatten().find(|c| !"0123?".contains(**c)) {
        return Err(Error::parse(format!(
            "unexpected character {bad:?} in MALIN table"
        )));
    }

    let cells_at = |site: usize| -> (Bitset, Bitset) {
        let mut presence = Bitset::new(taxa.len());
        let mut observable = Bitset::new(taxa.len());
        for (t, row) in rows.iter().enumerate() {
            match row[site] {
                '?' => {}
                '0' => {
                    observable.set(t);
                }
                _ => {
                    presence.set(t);
                    observable.set(t);
                }
            }
        }
        (presence, observable)
    };

    if map.is_empty() {
        if rows.iter().flatten().any(|c| *c == '2' || *c == '3') {
            return Err(Error::parse(
                "table contains phase codes 2/3 but no #MAP lines, so intron positions cannot be recovered",
            ));
        }
        let mut sites = Vec::new();
        for s in 0..n_sites {
            let (presence, observable) = cells_at(s);
            if presence.count_ones() > 0 {
                sites.push(Site {
                    col: s as u32,
                    phase: 0,
                    presence,
                    observable,
                });
            }
        }
        return Ok(vec![SiteMatrix::from_parts(
            fallback_name,
            taxa,
            sites,
            PhaseMode::Ignore,
        )]);
    }

    if let Some(missing) = (0..n_sites).find(|s| !map.contains_key(s)) {
        return Err(Error::parse(format!(
            "site {missing} has no #MAP line ({} sites, {} #MAP lines)",
            n_sites,
            map.len()
        )));
    }
    let mut per_aln: BTreeMap<String, Vec<Site>> = BTreeMap::new();
    for s in 0..n_sites {
        let (aln, alipos) = &map[&s];
        let (col, phase) = (alipos / 3, (alipos % 3) as u8);
        if let Some(c) = rows
            .iter()
            .map(|r| r[s])
            .find(|c| phase_of_code(*c).is_some_and(|p| p != phase))
        {
            return Err(Error::parse(format!(
                "site {s}: code {c} means phase {} but alipos {alipos} implies phase {phase}",
                phase_of_code(c).unwrap_or(0)
            )));
        }
        let (presence, observable) = cells_at(s);
        per_aln.entry(aln.clone()).or_default().push(Site {
            col,
            phase,
            presence,
            observable,
        });
    }
    Ok(per_aln
        .into_iter()
        .map(|(aln, mut sites)| {
            sites.sort_by_key(|s| (s.col, s.phase));
            sites.dedup_by(|dup, keep| {
                let same = (dup.col, dup.phase) == (keep.col, keep.phase);
                if same {
                    keep.presence.or_assign(&dup.presence);
                    keep.observable.or_assign(&dup.observable);
                }
                same
            });
            SiteMatrix::from_parts(&aln, taxa.clone(), sites, PhaseMode::Strict)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const TABLE: &str = "#ORGANISM A\tAlpha\n#ORGANISM B\tBeta\nA\t31?\nB\t302\n\
#MAP\tsite\tkog\talipos\t# header\n#MAP\t0\tK1\t3\n#MAP\t1\tK1\t7\n#MAP\t2\tK2\t8\n";

    #[test]
    fn import_maps_alipos_to_column_and_phase() {
        let ms = from_malin(TABLE, "x").unwrap();
        assert_eq!(ms.len(), 2);
        let k1 = &ms[0];
        assert_eq!(k1.alignment, "K1");
        assert_eq!(
            k1.sites
                .iter()
                .map(|s| (s.col, s.phase))
                .collect::<Vec<_>>(),
            vec![(1, 0), (2, 1)]
        );
        assert_eq!(k1.pattern(0), "11");
        assert_eq!(k1.pattern(1), "10");
        let k2 = &ms[1];
        assert_eq!((k2.sites[0].col, k2.sites[0].phase), (2, 2));
        assert_eq!(k2.pattern(0), "?1", "A is `?` at that site");
    }

    #[test]
    fn export_then_import_is_lossless() {
        let ms = from_malin(TABLE, "x").unwrap();
        let refs: Vec<&SiteMatrix> = ms.iter().collect();
        let text = to_malin(&refs).unwrap();
        assert!(text.contains("A\t31?"), "{text}");
        assert!(text.contains("B\t302"), "{text}");
        let back = from_malin(&text, "x").unwrap();
        assert_eq!(back.len(), 2);
        for (a, b) in ms.iter().zip(&back) {
            assert_eq!(a.alignment, b.alignment);
            assert_eq!(a.sites.len(), b.sites.len());
            for i in 0..a.sites.len() {
                assert_eq!(
                    (a.sites[i].col, a.sites[i].phase),
                    (b.sites[i].col, b.sites[i].phase)
                );
                assert_eq!(a.pattern(i), b.pattern(i));
            }
        }
    }

    #[test]
    fn a_phase_code_that_contradicts_alipos_is_rejected() {
        let bad = "A\t1\nB\t0\n#MAP\t0\tK\t3\n";
        let err = from_malin(bad, "x").unwrap_err().to_string();
        assert!(err.contains("implies phase"), "{err}");
    }

    #[test]
    fn binary_table_imports_without_phase_and_cannot_be_exported() {
        let ms = from_malin("A\t0110\nB\t0001\n", "bin").unwrap();
        assert_eq!(ms.len(), 1);
        assert_eq!(ms[0].phase_mode, PhaseMode::Ignore);
        assert_eq!(
            ms[0].sites.iter().map(|s| s.col).collect::<Vec<_>>(),
            vec![1, 2, 3]
        );
        assert!(to_malin(&[&ms[0]]).is_err());
    }

    #[test]
    fn phase_codes_without_a_map_are_refused() {
        assert!(from_malin("A\t3\nB\t0\n", "x").is_err());
    }

    #[test]
    fn malformed_tables_error_instead_of_panicking() {
        assert!(from_malin("", "x").is_err());
        assert!(from_malin("A\t01\nB\t0\n", "x").is_err(), "ragged rows");
        assert!(from_malin("A\t0x\nB\t00\n", "x").is_err(), "bad character");
        assert!(
            from_malin("A\t00\nB\t00\n#MAP\t0\tK\t3\n", "x").is_err(),
            "site 1 unmapped"
        );
    }
}
