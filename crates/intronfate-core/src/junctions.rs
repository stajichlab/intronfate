//! RNA-seq splice-junction support as soft cell evidence (SPEC §5.6b).
//!
//! An annotated intron that reads span with the same splice coordinates is probably real; one that no read spans may be real
//! but unexpressed, or a gene-model error. [`soft_evidence`] turns junction read counts into the `site / taxon / p_present` table
//! that `--ml-soft-evidence` reads: each *present* cell of one alignment gets the probability that its intron is real.
//!
//! The probability of an intron with `n` supporting reads is
//!
//! ```text
//! p(n) = 1 − (1 − p0) · k / (n + k)
//! ```
//!
//! `p0` is the probability for an intron with no reads (the annotation is trusted that far and no further) and `k` is the
//! number of reads that halves the remaining doubt. `p` rises monotonically from `p0` at `n = 0` towards 1. A cell backed by
//! several introns of one taxon (paralogs, isoforms) takes the largest `p`: one supported copy is enough.
//!
//! Only cells the annotation already calls present are touched. A junction with reads at coordinates the annotation does
//! not have would need projecting from the genome into the alignment, which this does not do.

use std::collections::HashMap;
use std::fmt::Write as _;

use crate::error::{Error, Result};

/// Junction table layouts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// STAR `SJ.out.tab`: chromosome, first and last intron base (1-based), strand (0 undefined, 1 `+`, 2 `-`), motif,
    /// annotated, unique reads, multi-mapping reads, overhang.
    Star,
    /// BED6: chromosome, start (0-based), end, name, score = reads, strand.
    Bed,
}

impl std::str::FromStr for Format {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s {
            "star" => Ok(Format::Star),
            "bed" => Ok(Format::Bed),
            other => Err(format!(
                "unknown junction format {other:?}: use star or bed"
            )),
        }
    }
}

/// Reads per junction, keyed by `(seqid, first intron base, last intron base)` (1-based, inclusive) and strand
/// (`+1`, `-1`, `0` unknown).
#[derive(Clone, Debug, Default)]
pub struct JunctionReads {
    reads: HashMap<(String, u64, u64, i8), u64>,
}

impl JunctionReads {
    /// Parse a junction table. `count_multi` adds STAR's multi-mapping reads to the unique ones.
    pub fn parse(text: &str, format: Format, count_multi: bool) -> Result<Self> {
        let mut reads: HashMap<(String, u64, u64, i8), u64> = HashMap::new();
        for (n, line) in text.lines().enumerate() {
            let line = line.trim_end_matches('\r');
            if line.trim().is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let bad = |what: &str| Error::parse(format!("junction table line {}: {what}", n + 1));
            let num = |i: usize, what: &str| -> Result<u64> {
                f.get(i)
                    .and_then(|s| s.trim().parse::<u64>().ok())
                    .ok_or_else(|| bad(what))
            };
            let (seq, start, end, strand, n_reads) = match format {
                Format::Star => {
                    if f.len() < 7 {
                        return Err(bad("expected the 9 columns of STAR SJ.out.tab"));
                    }
                    let strand = match f[3].trim() {
                        "1" => 1,
                        "2" => -1,
                        _ => 0,
                    };
                    let mut r = num(6, "unique read count is not a number")?;
                    if count_multi {
                        r += num(7, "multi-mapping read count is not a number")?;
                    }
                    (
                        f[0],
                        num(1, "start is not a number")?,
                        num(2, "end is not a number")?,
                        strand,
                        r,
                    )
                }
                Format::Bed => {
                    if f.len() < 5 {
                        return Err(bad("expected at least 5 BED columns (score = reads)"));
                    }
                    let strand = match f.get(5).map(|s| s.trim()) {
                        Some("+") => 1,
                        Some("-") => -1,
                        _ => 0,
                    };
                    // BED starts are 0-based; the table here is 1-based inclusive.
                    (
                        f[0],
                        num(1, "start is not a number")? + 1,
                        num(2, "end is not a number")?,
                        strand,
                        num(4, "score (read count) is not a number")?,
                    )
                }
            };
            *reads
                .entry((seq.to_string(), start, end, strand))
                .or_insert(0) += n_reads;
        }
        Ok(JunctionReads { reads })
    }

    pub fn len(&self) -> usize {
        self.reads.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    /// Reads spanning exactly this intron: those on its own strand plus those whose strand is undefined.
    pub fn reads_for(&self, seq: &str, start: u64, end: u64, plus: bool) -> u64 {
        let own = if plus { 1 } else { -1 };
        [own, 0]
            .iter()
            .map(|&s| {
                self.reads
                    .get(&(seq.to_string(), start, end, s))
                    .copied()
                    .unwrap_or(0)
            })
            .sum()
    }
}

/// Probability that an annotated intron is real, given `reads` supporting reads (see the module notes).
pub fn p_present(reads: u64, p0: f64, k: f64) -> f64 {
    1.0 - (1.0 - p0) * k / (reads as f64 + k)
}

/// Settings of the read-count to probability map.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Probability of an intron with no supporting reads.
    pub p0: f64,
    /// Reads that halve the doubt.
    pub k: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options { p0: 0.8, k: 3.0 }
    }
}

/// One soft cell.
#[derive(Clone, Debug, PartialEq)]
pub struct SoftCell {
    pub site: String,
    pub taxon: String,
    pub p_present: f64,
    /// Reads of the best-supported intron behind the cell.
    pub reads: u64,
    /// Introns of this taxon behind the cell.
    pub introns: usize,
}

/// What a run did, for the summary line.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    pub cells: usize,
    pub supported: usize,
    pub unsupported: usize,
    /// Present cells of taxa that have a junction table but no row in `introns.tsv` (not softened).
    pub no_intron_record: usize,
    /// Present cells of taxa without a junction table (left hard).
    pub taxa_without_junctions: usize,
}

/// One intron record behind a cell: `(seqid, first base, last base, on + strand)`.
type IntronRecord<'a> = (&'a str, u64, u64, bool);

struct Table<'a> {
    rows: Vec<Vec<&'a str>>,
    head: Vec<&'a str>,
}

impl<'a> Table<'a> {
    fn parse(text: &'a str, what: &str) -> Result<Self> {
        let mut lines = text
            .lines()
            .map(|l| l.trim_end_matches('\r'))
            .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
        let head: Vec<&str> = lines
            .next()
            .ok_or_else(|| Error::parse(format!("{what} is empty")))?
            .split('\t')
            .collect();
        Ok(Table {
            rows: lines.map(|l| l.split('\t').collect()).collect(),
            head,
        })
    }

    fn col(&self, name: &str, what: &str) -> Result<usize> {
        self.head
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| Error::parse(format!("{what}: no `{name}` column")))
    }
}

/// Soft cells for one alignment.
///
/// `introns` is the run's `introns.tsv` (genomic coordinates), `sites` the alignment's `sites.tsv` (which cell each intron fills),
/// and `junctions` maps a taxon name to its junction reads. Cells of taxa not in `junctions` are left out, and stay hard.
pub fn soft_evidence(
    introns: &str,
    sites: &str,
    junctions: &[(String, JunctionReads)],
    opts: &Options,
) -> Result<(Vec<SoftCell>, Summary)> {
    if !(0.0..=1.0).contains(&opts.p0) || opts.p0.is_nan() {
        return Err(Error::parse("p0 must be in [0, 1]"));
    }
    if opts.k <= 0.0 || opts.k.is_nan() {
        return Err(Error::parse("k must be positive"));
    }
    let it = Table::parse(introns, "introns table")?;
    let (ip, isq, istr, ist, ien, ico) = (
        it.col("protein_id", "introns table")?,
        it.col("seqid", "introns table")?,
        it.col("strand", "introns table")?,
        it.col("genomic_start", "introns table")?,
        it.col("genomic_end", "introns table")?,
        it.col("cds_offset", "introns table")?,
    );
    // (protein, cds_offset) -> every intron record behind it (one per isoform that carries the intron)
    let mut by_site: HashMap<(&str, &str), Vec<IntronRecord>> = HashMap::new();
    for r in &it.rows {
        let get = |i: usize| {
            r.get(i)
                .copied()
                .ok_or_else(|| Error::parse("introns table: a row is short"))
        };
        let num = |i: usize| -> Result<u64> {
            get(i)?
                .parse()
                .map_err(|_| Error::parse("introns table: a coordinate is not a number"))
        };
        by_site.entry((get(ip)?, get(ico)?)).or_default().push((
            get(isq)?,
            num(ist)?,
            num(ien)?,
            get(istr)? != "-",
        ));
    }

    let st = Table::parse(sites, "sites table")?;
    let (sa, sx, sp, sc, sl) = (
        st.col("alignment", "sites table")?,
        st.col("taxon", "sites table")?,
        st.col("protein_id", "sites table")?,
        st.col("cds_offset", "sites table")?,
        st.col("site_label", "sites table")?,
    );
    let alignments: std::collections::BTreeSet<&str> =
        st.rows.iter().filter_map(|r| r.get(sa).copied()).collect();
    if alignments.len() > 1 {
        return Err(Error::parse(format!(
            "the sites table covers {} alignments ({}); site labels are per alignment, so give one alignment's sites.tsv",
            alignments.len(),
            alignments.iter().take(3).copied().collect::<Vec<_>>().join(", ")
        )));
    }
    let by_taxon: HashMap<&str, &JunctionReads> =
        junctions.iter().map(|(t, j)| (t.as_str(), j)).collect();

    let mut summary = Summary::default();
    let mut cells: Vec<SoftCell> = Vec::new();
    let mut at: HashMap<(String, String), usize> = HashMap::new();
    for r in &st.rows {
        let get = |i: usize| {
            r.get(i)
                .copied()
                .ok_or_else(|| Error::parse("sites table: a row is short"))
        };
        let taxon = get(sx)?;
        let Some(jr) = by_taxon.get(taxon) else {
            summary.taxa_without_junctions += 1;
            continue;
        };
        let Some(recs) = by_site.get(&(get(sp)?, get(sc)?)) else {
            summary.no_intron_record += 1;
            continue;
        };
        let best = recs
            .iter()
            .map(|&(seq, s, e, plus)| jr.reads_for(seq, s, e, plus))
            .max()
            .unwrap_or(0);
        let p = p_present(best, opts.p0, opts.k);
        let key = (get(sl)?.to_string(), taxon.to_string());
        match at.get(&key) {
            Some(&i) => {
                let c = &mut cells[i];
                c.introns += recs.len();
                if best > c.reads {
                    c.reads = best;
                    c.p_present = p;
                }
            }
            None => {
                at.insert(key.clone(), cells.len());
                cells.push(SoftCell {
                    site: key.0,
                    taxon: key.1,
                    p_present: p,
                    reads: best,
                    introns: recs.len(),
                });
            }
        }
    }
    summary.cells = cells.len();
    summary.supported = cells.iter().filter(|c| c.reads > 0).count();
    summary.unsupported = summary.cells - summary.supported;
    Ok((cells, summary))
}

/// The table `--ml-soft-evidence` reads, with the parameters in a `#` line.
pub fn soft_table(cells: &[SoftCell], opts: &Options, source: &str) -> String {
    let mut s = String::new();
    let _ = writeln!(
        s,
        "# intronfate junctions: p_present = 1 - (1 - p0) * k / (reads + k), p0={}, k={}, {source}",
        opts.p0, opts.k
    );
    s.push_str("site\ttaxon\tp_present\treads\tintrons\n");
    for c in cells {
        let _ = writeln!(
            s,
            "{}\t{}\t{:.6}\t{}\t{}",
            c.site, c.taxon, c.p_present, c.reads, c.introns
        );
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    const INTRONS: &str = "protein_id\tgene_id\ttranscript_id\tseqid\tstrand\tintron_rank\tcoding_rank\tgenomic_start\tgenomic_end\tintron_len\tcds_offset\tphase\tcodon\n\
        pA\tg\tt\tchrA\t+\t0\t0\t113\t200\t88\t12\t0\t4\n\
        pA\tg\tt\tchrA\t+\t1\t1\t210\t300\t91\t21\t0\t7\n\
        pB\tg\tt\tchrB\t-\t0\t0\t500\t600\t101\t12\t0\t4\n";
    const SITES: &str = "alignment\ttaxon\tprotein_id\tcolumn\tphase\tcds_offset\tsite_label\n\
        fam\tA\tpA\t4\t0\t12\t4:0\n\
        fam\tA\tpA\t7\t0\t21\t7:0\n\
        fam\tB\tpB\t4\t0\t12\t4:0\n";

    fn star(rows: &[(&str, u64, u64, u8, u64, u64)]) -> JunctionReads {
        let t: String = rows
            .iter()
            .map(|(c, s, e, st, u, m)| format!("{c}\t{s}\t{e}\t{st}\t1\t1\t{u}\t{m}\t30\n"))
            .collect();
        JunctionReads::parse(&t, Format::Star, false).unwrap()
    }

    #[test]
    fn probability_rises_from_p0_towards_one_and_halves_the_doubt_at_k() {
        let (p0, k) = (0.8, 3.0);
        assert!((p_present(0, p0, k) - p0).abs() < 1e-12);
        assert!(
            (p_present(3, p0, k) - 0.9).abs() < 1e-12,
            "doubt 0.2 halves at k reads"
        );
        let mut last = 0.0;
        for n in 0..50 {
            let p = p_present(n, p0, k);
            assert!(p > last - 1e-15 && p <= 1.0);
            last = p;
        }
        assert!(p_present(1_000_000, p0, k) > 0.9999);
        // an unsupported intron is never less probable than p0: reads cannot lower it
        assert!(p_present(1, p0, k) > p_present(0, p0, k));
    }

    #[test]
    fn star_and_bed_tables_parse_with_their_own_coordinates_and_strands() {
        let j = JunctionReads::parse(
            "chrA\t113\t200\t1\t1\t1\t7\t2\t30\nchrA\t210\t300\t0\t0\t0\t4\t0\t20\n",
            Format::Star,
            false,
        )
        .unwrap();
        assert_eq!(j.len(), 2);
        assert_eq!(j.reads_for("chrA", 113, 200, true), 7);
        assert_eq!(
            j.reads_for("chrA", 113, 200, false),
            0,
            "defined strand must agree"
        );
        assert_eq!(
            j.reads_for("chrA", 210, 300, true),
            4,
            "undefined strand counts for either"
        );
        assert_eq!(j.reads_for("chrA", 210, 300, false), 4);
        let multi = JunctionReads::parse("chrA\t113\t200\t1\t1\t1\t7\t2\t30\n", Format::Star, true)
            .unwrap();
        assert_eq!(multi.reads_for("chrA", 113, 200, true), 9, "unique + multi");
        // BED: 0-based start becomes 1-based, so 112..200 is the intron 113..200
        let b = JunctionReads::parse(
            "track name=x\nchrB\t499\t600\tj1\t5\t-\n",
            Format::Bed,
            false,
        )
        .unwrap();
        assert_eq!(b.reads_for("chrB", 500, 600, false), 5);
        assert!(JunctionReads::parse("chrA\t1\n", Format::Star, false).is_err());
        assert!(JunctionReads::parse("chrA\t1\t2\t.\tx\t4\n", Format::Bed, false).is_err());
    }

    #[test]
    fn supported_and_unsupported_introns_get_their_probabilities() {
        let a = star(&[("chrA", 113, 200, 1, 9, 0)]); // first intron has 9 reads, second none
        let (cells, s) =
            soft_evidence(INTRONS, SITES, &[("A".into(), a)], &Options::default()).unwrap();
        assert_eq!(
            cells.len(),
            2,
            "taxon B has no junction table and stays hard"
        );
        let first = &cells[0];
        assert_eq!(
            (first.site.as_str(), first.taxon.as_str(), first.reads),
            ("4:0", "A", 9)
        );
        assert!((first.p_present - p_present(9, 0.8, 3.0)).abs() < 1e-12);
        let second = &cells[1];
        assert_eq!((second.site.as_str(), second.reads), ("7:0", 0));
        assert!((second.p_present - 0.8).abs() < 1e-12);
        assert_eq!(
            s,
            Summary {
                cells: 2,
                supported: 1,
                unsupported: 1,
                no_intron_record: 0,
                taxa_without_junctions: 1
            }
        );
    }

    #[test]
    fn a_wrong_strand_or_shifted_coordinate_is_not_support() {
        // B's intron is on `-`; reads on `+` at the same coordinates do not support it, a one-base shift neither
        let b = star(&[("chrB", 500, 600, 1, 50, 0), ("chrB", 501, 600, 2, 50, 0)]);
        let (cells, _) =
            soft_evidence(INTRONS, SITES, &[("B".into(), b)], &Options::default()).unwrap();
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].reads, 0);
        let b = star(&[("chrB", 500, 600, 2, 6, 0)]);
        let (cells, _) =
            soft_evidence(INTRONS, SITES, &[("B".into(), b)], &Options::default()).unwrap();
        assert_eq!(cells[0].reads, 6);
    }

    #[test]
    fn several_introns_behind_one_cell_take_the_best_supported() {
        let introns = format!("{INTRONS}pA\tg\tt2\tchrA\t+\t0\t0\t113\t199\t87\t12\t0\t4\n");
        let a = star(&[("chrA", 113, 199, 1, 20, 0), ("chrA", 113, 200, 1, 2, 0)]);
        let (cells, _) =
            soft_evidence(&introns, SITES, &[("A".into(), a)], &Options::default()).unwrap();
        let c = cells.iter().find(|c| c.site == "4:0").unwrap();
        assert_eq!((c.reads, c.introns), (20, 2));
    }

    #[test]
    fn malformed_inputs_and_settings_are_refused() {
        let j = star(&[]);
        let two = SITES.replace("fam\tB", "other\tB");
        let e = soft_evidence(
            INTRONS,
            &two,
            &[("A".into(), j.clone())],
            &Options::default(),
        )
        .unwrap_err();
        assert!(e.to_string().contains("alignments"), "{e}");
        assert!(soft_evidence("a\tb\n1\t2\n", SITES, &[], &Options::default()).is_err());
        assert!(soft_evidence(INTRONS, SITES, &[], &Options { p0: 1.5, k: 3.0 }).is_err());
        assert!(soft_evidence(INTRONS, SITES, &[], &Options { p0: 0.8, k: 0.0 }).is_err());
    }

    #[test]
    fn the_table_is_what_the_soft_evidence_reader_expects() {
        let cells = vec![SoftCell {
            site: "4:0".into(),
            taxon: "A".into(),
            p_present: 0.9,
            reads: 3,
            introns: 1,
        }];
        let t = soft_table(&cells, &Options::default(), "x");
        let mut lines = t.lines();
        assert!(lines.next().unwrap().starts_with('#'));
        assert_eq!(
            lines.next().unwrap(),
            "site\ttaxon\tp_present\treads\tintrons"
        );
        assert_eq!(lines.next().unwrap(), "4:0\tA\t0.900000\t3\t1");
    }
}
