//! The `site × taxon` presence matrix and the sharing statistics computed from it.
//!
//! A site is one alignment column plus a codon phase; a taxon carries the site when its own
//! annotation puts an intron at that position. Absence is only meaningful for taxa that were
//! annotated at all, so [`SiteMatrix`] tracks which columns of the matrix are real observations and
//! which are "no data" — conflating the two is the classic way to invent intron losses.
//!
//! Presence sets are bitsets over taxa, so the operations that matter here (site counts, pairwise
//! intersections, pattern strings) are single-pass and cache-friendly.

use std::fmt::Write as _;
use std::path::Path;

use crate::error::{Error, Result};
use crate::model::Phase;
use crate::project::{cell_observable, AlignmentProjection};

/// A fixed-size 0/1 set over taxon indices.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Bitset {
    words: Vec<u64>,
    len: usize,
    ones: u32,
}

impl Bitset {
    pub fn new(len: usize) -> Self {
        Bitset {
            words: vec![0; len.div_ceil(64)],
            len,
            ones: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.ones == 0
    }

    pub fn count_ones(&self) -> u32 {
        self.ones
    }

    /// Sets bit `i`, returning `false` when it was already set.
    pub fn set(&mut self, i: usize) -> bool {
        if i >= self.len {
            debug_assert!(false, "index {i} out of range for bitset of {}", self.len);
            return false;
        }
        let mask = 1u64 << (i % 64);
        let w = &mut self.words[i / 64];
        if *w & mask != 0 {
            return false;
        }
        *w |= mask;
        self.ones += 1;
        true
    }

    pub fn test(&self, i: usize) -> bool {
        i < self.len && self.words[i / 64] & (1u64 << (i % 64)) != 0
    }

    /// Number of taxa carrying both sites: the pairwise intron-position intersection.
    pub fn and_count(&self, other: &Bitset) -> u32 {
        self.words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| (a & b).count_ones())
            .sum()
    }

    /// Set bit positions in ascending order.
    pub fn iter_ones(&self) -> impl Iterator<Item = usize> + '_ {
        let len = self.len;
        self.words.iter().enumerate().flat_map(move |(w, &word)| {
            let base = w * 64;
            (0..64)
                .filter(move |b| word >> b & 1 == 1)
                .map(move |b| base + b)
                .filter(move |&i| i < len)
        })
    }

    /// In-place union with a bitset of the same length.
    pub fn or_assign(&mut self, other: &Bitset) {
        debug_assert_eq!(self.len, other.len);
        for (a, b) in self.words.iter_mut().zip(other.words.iter()) {
            *a |= *b;
        }
        self.ones = self.words.iter().map(|w| w.count_ones()).sum();
    }

    /// Set positions as a vector; used by table writers, not inner loops.
    pub fn positions(&self) -> Vec<usize> {
        self.iter_ones().collect()
    }
}

/// One site: an alignment column with a phase, and the taxa that carry it.
#[derive(Clone, Debug)]
pub struct Site {
    /// 0-based alignment column.
    pub col: u32,
    /// Codon phase used for grouping (always 0 under `PhaseMode::Ignore`).
    pub phase: Phase,
    /// Taxa carrying an intron here. Always a subset of `observable`.
    pub presence: Bitset,
    /// Taxa for which this site is an *observation*: annotated, and aligned around the site. A taxon
    /// outside this set has a missing (`?`) cell, never an absent (`0`) one.
    pub observable: Bitset,
}

/// State of one `(site, taxon)` cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cell {
    /// Observed, intron present.
    Present,
    /// Observed, no intron.
    Absent,
    /// No observation: unannotated taxon, or not aligned around the site.
    Missing,
}

impl Cell {
    /// Matrix character: `1`, `0` or `?`.
    pub fn symbol(self) -> char {
        match self {
            Cell::Present => '1',
            Cell::Absent => '0',
            Cell::Missing => '?',
        }
    }
}

impl Site {
    pub fn n_observable(&self) -> u32 {
        self.observable.count_ones()
    }

    pub fn cell(&self, taxon: usize) -> Cell {
        if self.presence.test(taxon) {
            Cell::Present
        } else if self.observable.test(taxon) {
            Cell::Absent
        } else {
            Cell::Missing
        }
    }

    pub fn n_taxa(&self) -> u32 {
        self.presence.count_ones()
    }

    pub fn taxa(&self) -> Vec<usize> {
        self.presence.positions()
    }

    pub fn has(&self, taxon: usize) -> bool {
        self.presence.test(taxon)
    }

    /// Label used in tables: `<column>:<phase>`, e.g. `128:1`.
    pub fn label(&self) -> String {
        format!("{}:{}", self.col, self.phase)
    }
}

/// Rows = sites, columns = taxa.
#[derive(Clone, Debug)]
pub struct SiteMatrix {
    pub alignment: String,
    /// Alignment length in columns.
    pub alen: usize,
    pub taxa: Vec<String>,
    /// Per taxon: did the annotation contribute introns for this row?
    pub annotated: Vec<bool>,
    pub sites: Vec<Site>,
    /// Site indices per taxon, ascending in `(column, phase)`.
    pub taxon_sites: Vec<Vec<u32>>,
    /// Every "fraction of taxa" statistic divides by this, not by `taxa.len()`.
    pub n_annotated: usize,
    /// Residues each row contributes to the alignment, for reporting.
    pub aligned_residues: Vec<usize>,
    /// Which phase grouping produced this matrix, so a re-read file is self-describing.
    pub phase_mode: crate::project::PhaseMode,
    /// Flank used for per-cell observability ([`crate::project::cell_observable`]).
    pub flank: u32,
    /// Annotated introns that were *not* recorded as present because the taxon was not aligned
    /// around the site. Reported rather than silently promoted to observed: counting presence while
    /// discarding absence in poorly aligned regions would bias towards gains.
    pub masked_present: u64,
}

/// Everything `pairwise` precomputes once per matrix.
struct PairIndex {
    /// Per taxon, ascending site columns split by phase (index 0..3).
    cols: Vec<[Vec<u32>; 3]>,
    /// Per taxon, all ascending site columns (for phase-blind matching).
    all: Vec<Vec<u32>>,
    /// Per taxon, a bitset over sites: observed (`1` or `0`) at that site.
    observed: Vec<Bitset>,
}

/// Greedy one-to-one matching of two ascending column lists within `tolerance` columns, so one site
/// can never support two shared calls.
fn greedy_match(a: &[u32], b: &[u32], tolerance: u32) -> u32 {
    let (mut i, mut j, mut shared) = (0usize, 0usize, 0u32);
    while i < a.len() && j < b.len() {
        if a[i].abs_diff(b[j]) <= tolerance {
            shared += 1;
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    shared
}

impl SiteMatrix {
    /// Assemble a matrix from finished sites (importers, tests). Taxa that are observed at no site are
    /// treated as unannotated, so their cells stay missing rather than becoming absences.
    pub fn from_parts(
        alignment: &str,
        taxa: Vec<String>,
        sites: Vec<Site>,
        phase_mode: crate::project::PhaseMode,
    ) -> SiteMatrix {
        let n_taxa = taxa.len();
        let annotated: Vec<bool> = (0..n_taxa)
            .map(|t| sites.iter().any(|s| s.observable.test(t)))
            .collect();
        let n_annotated = annotated.iter().filter(|&&a| a).count();
        let mut taxon_sites = vec![Vec::new(); n_taxa];
        for (i, site) in sites.iter().enumerate() {
            for t in site.presence.positions() {
                taxon_sites[t].push(i as u32);
            }
        }
        SiteMatrix {
            alignment: alignment.to_string(),
            alen: 0,
            taxa,
            annotated,
            sites,
            taxon_sites,
            n_annotated,
            aligned_residues: Vec::new(),
            phase_mode,
            flank: 0,
            masked_present: 0,
        }
    }

    /// Collapse projected column sites into unique `(column, phase)` rows.
    ///
    /// Requires `ap.sites` sorted by `(col, phase, taxon)`, which
    /// [`project_into_alignment`](crate::project::project_into_alignment) guarantees.
    pub fn build(ap: &AlignmentProjection) -> SiteMatrix {
        let n_taxa = ap.taxa.len();
        let annotated: Vec<bool> = ap.taxa.iter().map(|t| t.has_annotation).collect();
        let n_annotated = annotated.iter().filter(|&&b| b).count();
        let mut sites: Vec<Site> = Vec::new();
        let mut taxon_sites = vec![Vec::new(); n_taxa];
        let mut masked_present = 0u64;

        for cs in &ap.sites {
            let taxon = cs.taxon as usize;
            if taxon >= n_taxa {
                continue;
            }
            // A carrier is judged by the window of its *own* phase (also under --ignore-phase, where the
            // grouping phase is always 0), and a clamped carrier is exempt: its right-hand flank lies
            // past the end of the row by construction.
            if let (Some(res), false) = (&ap.residue_at, cs.clamped) {
                if !cell_observable(res, taxon, cs.col, Some(cs.site.phase), ap.flank) {
                    masked_present += 1;
                    continue;
                }
            }
            if sites
                .last()
                .is_none_or(|s| (s.col, s.phase) != (cs.col, cs.phase))
            {
                sites.push(Site {
                    col: cs.col,
                    phase: cs.phase,
                    presence: Bitset::new(n_taxa),
                    observable: Bitset::new(n_taxa),
                });
            }
            let idx = sites.len() - 1;
            // `set` returns false only if this taxon already carries the site (possible under
            // TailPolicy::Clamp, where two protein-space sites collapse onto one column).
            if sites[idx].presence.set(taxon) {
                taxon_sites[taxon].push(idx as u32);
            }
        }

        let merged = matches!(ap.phase_mode, crate::project::PhaseMode::Ignore);
        for site in &mut sites {
            // Absences are only observations where the taxon covers the window. A merged site
            // (--ignore-phase) pools introns of every phase, so its window is the union of theirs.
            let window_phase = if merged { None } else { Some(site.phase) };
            for (t, &is_annotated) in annotated.iter().enumerate() {
                let seen = is_annotated
                    && ap.residue_at.as_ref().is_none_or(|res| {
                        cell_observable(res, t, site.col, window_phase, ap.flank)
                    });
                if seen {
                    site.observable.set(t);
                }
            }
            // Carriers already passed their own window above (or are clamped), so they are observed.
            let carriers = site.presence.clone();
            site.observable.or_assign(&carriers);
            debug_assert!(
                site.presence
                    .positions()
                    .iter()
                    .all(|&t| site.observable.test(t)),
                "presence must be a subset of observable"
            );
        }

        SiteMatrix {
            alignment: ap.alignment.clone(),
            alen: ap.alen,
            taxa: ap.taxa.iter().map(|t| t.taxon_name.clone()).collect(),
            annotated,
            sites,
            taxon_sites,
            n_annotated,
            aligned_residues: ap.taxa.iter().map(|t| t.residues_aligned).collect(),
            phase_mode: ap.phase_mode,
            flank: ap.flank,
            masked_present,
        }
    }

    /// Read a `matrix.tsv` written by [`crate::tables::matrix_table`], so scoring can run as a
    /// separate step from projection without re-parsing the annotation.
    ///
    /// The leading `#` metadata line carries alignment name, alignment length, phase mode and which
    /// taxa were unannotated; without it those absences would be indistinguishable from losses.
    pub fn read_tsv(path: impl AsRef<Path>) -> Result<SiteMatrix> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)?;
        let mut meta: Vec<(String, String)> = Vec::new();
        let mut header: Option<Vec<String>> = None;
        let mut rows: Vec<Vec<String>> = Vec::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.trim_end();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("# ") {
                meta.extend(rest.split_whitespace().filter_map(|kv| {
                    kv.split_once('=')
                        .map(|(k, v)| (k.to_string(), v.to_string()))
                }));
                continue;
            }
            let cells: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            if header.is_none() {
                header = Some(cells);
            } else {
                rows.push(cells);
            }
            let _ = i;
        }
        let header =
            header.ok_or_else(|| Error::parse(format!("{}: file is empty", path.display())))?;
        let get = |key: &str| meta.iter().find(|(k, _)| k == key).map(|(_, v)| v.clone());
        // Taxon columns follow the `pattern` column. Find it by name so adding a summary column does
        // not silently shift every taxon by one.
        let fixed = header
            .iter()
            .position(|h| h == "pattern")
            .map(|i| i + 1)
            .ok_or_else(|| {
                Error::parse(format!(
                    "{}: header has no `pattern` column",
                    path.display()
                ))
            })?;
        let col_ix = header.iter().position(|h| h == "column").ok_or_else(|| {
            Error::parse(format!("{}: header has no `column` column", path.display()))
        })?;
        let phase_ix = header.iter().position(|h| h == "phase").ok_or_else(|| {
            Error::parse(format!("{}: header has no `phase` column", path.display()))
        })?;
        let taxa: Vec<String> = header[fixed..].to_vec();
        if let Some(dup) = taxa
            .iter()
            .enumerate()
            .find(|(i, t)| taxa[..*i].contains(t))
        {
            return Err(Error::parse(format!(
                "{}: taxon {:?} appears twice in the header",
                path.display(),
                dup.1
            )));
        }
        let n_taxa = taxa.len();
        if n_taxa == 0 {
            return Err(Error::parse(format!(
                "{}: matrix has no taxon columns",
                path.display()
            )));
        }
        let annotated_names: Vec<String> = get("unannotated")
            .filter(|v| v != "-")
            .map(|s| s.split(',').map(crate::util::decode_list_item).collect())
            .unwrap_or_default();
        let annotated: Vec<bool> = taxa
            .iter()
            .map(|t| !annotated_names.iter().any(|n| n == t))
            .collect();
        let n_annotated = annotated.iter().filter(|&&b| b).count();
        let mut parsed: Vec<Site> = Vec::with_capacity(rows.len());
        for (r, row) in rows.iter().enumerate() {
            if row.len() != header.len() {
                return Err(Error::parse(format!(
                    "{}: row {} has {} cells, expected {}",
                    path.display(),
                    r + 2,
                    row.len(),
                    header.len()
                )));
            }
            let col: u32 = row[col_ix].parse().map_err(|_| {
                Error::parse(format!("{}: bad column {:?}", path.display(), row[col_ix]))
            })?;
            let phase: u8 = row[phase_ix]
                .parse()
                .ok()
                .filter(|p| *p <= 2)
                .ok_or_else(|| {
                    Error::parse(format!(
                        "{}: bad phase {:?} (expected 0, 1 or 2)",
                        path.display(),
                        row[phase_ix]
                    ))
                })?;
            let mut presence = Bitset::new(n_taxa);
            let mut observable = Bitset::new(n_taxa);
            for t in 0..n_taxa {
                match row[fixed + t].as_str() {
                    "1" => {
                        if !annotated[t] {
                            return Err(Error::parse(format!(
                                "{}: taxon {:?} is listed as unannotated but carries site {}:{}",
                                path.display(),
                                taxa[t],
                                col,
                                phase
                            )));
                        }
                        presence.set(t);
                        observable.set(t);
                    }
                    // `0` in an unannotated row is the pre-`?` format's way of writing "no data".
                    "0" => {
                        if annotated[t] {
                            observable.set(t);
                        }
                    }
                    "?" => {}
                    other => {
                        return Err(Error::parse(format!(
                            "{}: cell for taxon {:?} is {other:?}, expected 1, 0 or ?",
                            path.display(),
                            taxa[t]
                        )))
                    }
                }
            }
            parsed.push(Site {
                col,
                phase,
                presence,
                observable,
            });
        }
        // Scoring and pairwise matching assume ascending, unique (column, phase); a hand-edited or
        // converted file need not honour that.
        parsed.sort_by_key(|s| (s.col, s.phase));
        let before = parsed.len();
        parsed.dedup_by(|dup, keep| {
            let same = (dup.col, dup.phase) == (keep.col, keep.phase);
            if same {
                keep.presence.or_assign(&dup.presence);
                keep.observable.or_assign(&dup.observable);
            }
            same
        });
        let _ = before;
        let sites = parsed;
        let mut taxon_sites = vec![Vec::new(); n_taxa];
        for (r, site) in sites.iter().enumerate() {
            for t in site.presence.positions() {
                taxon_sites[t].push(r as u32);
            }
        }
        Ok(SiteMatrix {
            alignment: get("alignment").unwrap_or_else(|| path.display().to_string()),
            alen: get("alen").and_then(|v| v.parse().ok()).unwrap_or(0),
            taxa,
            annotated,
            n_annotated,
            sites,
            taxon_sites,
            aligned_residues: get("residues")
                .map(|s| s.split(',').filter_map(|v| v.parse().ok()).collect())
                .unwrap_or_default(),
            phase_mode: match get("phase_mode").as_deref() {
                Some("ignore") => crate::project::PhaseMode::Ignore,
                _ => crate::project::PhaseMode::Strict,
            },
            flank: get("flank").and_then(|v| v.parse().ok()).unwrap_or(0),
            masked_present: get("masked_present")
                .and_then(|v| v.parse().ok())
                .unwrap_or(0),
        })
    }

    pub fn n_sites(&self) -> usize {
        self.sites.len()
    }

    pub fn n_taxa(&self) -> usize {
        self.taxa.len()
    }

    /// Intron positions carried by each taxon.
    pub fn per_taxon_counts(&self) -> Vec<u32> {
        self.taxon_sites.iter().map(|v| v.len() as u32).collect()
    }

    /// `1`/`0` per taxon in matrix order — the quick eyeball column of every site table.
    pub fn pattern(&self, site: usize) -> String {
        let mut s = String::with_capacity(self.taxa.len());
        for t in 0..self.taxa.len() {
            s.push(self.sites[site].cell(t).symbol());
        }
        s
    }

    pub fn taxa_names_at(&self, site: usize) -> Vec<&str> {
        self.sites[site]
            .taxa()
            .into_iter()
            .filter_map(|t| self.taxa.get(t).map(|s| s.as_str()))
            .collect()
    }

    /// Histogram: `spectrum[k]` = number of sites present in exactly `k` taxa.
    ///
    /// This is the shape people actually look at: a pile-up at `k = 1` means most positions are
    /// lineage-specific, a long right tail is the signal of shared ancestral introns.
    pub fn conservation_spectrum(&self) -> Vec<u64> {
        let mut h = vec![0u64; self.taxa.len() + 1];
        for s in &self.sites {
            h[s.n_taxa() as usize] += 1;
        }
        h
    }

    /// Sites carried by every annotated taxon.
    pub fn ubiquitous_sites(&self) -> Vec<usize> {
        (0..self.sites.len())
            .filter(|&i| {
                self.n_annotated > 0 && self.sites[i].n_taxa() as usize == self.n_annotated
            })
            .collect()
    }

    /// Sites unique to one taxon: candidates for lineage-specific gain (or everyone else's loss).
    pub fn singleton_sites(&self) -> Vec<usize> {
        (0..self.sites.len())
            .filter(|&i| self.sites[i].n_taxa() == 1)
            .collect()
    }

    /// Mean intron positions per annotated taxon.
    pub fn mean_sites_per_annotated_taxon(&self) -> f64 {
        if self.n_annotated == 0 {
            return 0.0;
        }
        self.taxon_sites
            .iter()
            .zip(&self.annotated)
            .filter(|(_, a)| **a)
            .map(|(v, _)| v.len())
            .sum::<usize>() as f64
            / self.n_annotated as f64
    }

    /// Pairwise sharing over the taxon upper triangle.
    ///
    /// `tolerance` is in alignment columns: `0` requires the identical column, `1` also accepts a
    /// site one column away, which is how indel-adjacent position drift is normally absorbed.
    pub fn pairwise(&self, tolerance: u32, match_phase: bool) -> Vec<PairStat> {
        use rayon::prelude::*;
        let t = self.n_taxa();
        // Per-taxon column lists and per-taxon observability over sites are built once, so each pair
        // costs two sorted merges and S/64 word ANDs instead of rebuilding lists and testing S bits.
        let idx = self.pair_index();
        (0..t)
            .into_par_iter()
            .flat_map_iter(|a| {
                let idx = &idx;
                ((a + 1)..t).map(move |b| self.pair_from_index(idx, a, b, tolerance, match_phase))
            })
            .collect()
    }

    fn pair_index(&self) -> PairIndex {
        let t = self.n_taxa();
        let mut cols = vec![[Vec::new(), Vec::new(), Vec::new()]; t];
        let mut all = vec![Vec::new(); t];
        for (tx, ids) in self.taxon_sites.iter().enumerate() {
            for &k in ids {
                let s = &self.sites[k as usize];
                cols[tx][usize::from(s.phase.min(2))].push(s.col);
                all[tx].push(s.col);
            }
        }
        let mut observed = vec![Bitset::new(self.sites.len()); t];
        for (i, s) in self.sites.iter().enumerate() {
            for tx in s.observable.iter_ones() {
                observed[tx].set(i);
            }
        }
        PairIndex {
            cols,
            all,
            observed,
        }
    }

    fn pair_from_index(
        &self,
        idx: &PairIndex,
        a: usize,
        b: usize,
        tolerance: u32,
        match_phase: bool,
    ) -> PairStat {
        let shared = if match_phase {
            (0..3)
                .map(|p| greedy_match(&idx.cols[a][p], &idx.cols[b][p], tolerance))
                .sum()
        } else {
            greedy_match(&idx.all[a], &idx.all[b], tolerance)
        };
        self.make_pair(a, b, shared, idx.observed[a].and_count(&idx.observed[b]))
    }

    fn make_pair(&self, a: usize, b: usize, shared: u32, n_both_observable: u32) -> PairStat {
        let (na, nb) = (
            self.taxon_sites[a].len() as u32,
            self.taxon_sites[b].len() as u32,
        );
        PairStat {
            a,
            b,
            taxon_a: self.taxa[a].clone(),
            taxon_b: self.taxa[b].clone(),
            n_a: na,
            n_b: nb,
            shared,
            jaccard: match na + nb {
                0 => 0.0,
                t => shared as f64 / (t - shared) as f64,
            },
            fraction_of_smaller: match na.min(nb) {
                0 => 0.0,
                m => shared as f64 / m as f64,
            },
            both_annotated: self.annotated[a] && self.annotated[b],
            n_both_observable,
        }
    }

    /// Dense `n × n` matrix of shared-position counts (symmetric; diagonal = own count).
    pub fn pairwise_matrix(&self, tolerance: u32, match_phase: bool) -> Vec<u32> {
        self.matrix_from_pairs(&self.pairwise(tolerance, match_phase))
    }

    /// Dense symmetric `shared` matrix built from already-computed pairs, so a caller that also writes
    /// `pairwise.tsv` does not compute every pair twice.
    pub fn matrix_from_pairs(&self, pairs: &[PairStat]) -> Vec<u32> {
        let n = self.n_taxa();
        let mut m = vec![0u32; n * n];
        for p in pairs {
            m[p.a * n + p.b] = p.shared;
            m[p.b * n + p.a] = p.shared;
        }
        for i in 0..n {
            m[i * n + i] = self.taxon_sites[i].len() as u32;
        }
        m
    }

    /// Greedy two-pointer matching over column order, so one site supports at most one shared call.
    pub fn pair_stat(&self, a: usize, b: usize, tolerance: u32, match_phase: bool) -> PairStat {
        let sa = &self.taxon_sites[a];
        let sb = &self.taxon_sites[b];
        let cols = |ids: &[u32], phase: Option<u8>| -> Vec<u32> {
            ids.iter()
                .map(|&k| &self.sites[k as usize])
                .filter(|x| phase.is_none_or(|p| x.phase == p))
                .map(|x| x.col)
                .collect()
        };
        // When phase must agree, match each phase separately: a pooled pass over columns cannot tell
        // "A has 9:0 and 9:1, B only 9:1" apart from a mismatch and would drop the shared 9:1.
        let shared = if match_phase {
            (0u8..3)
                .map(|p| greedy_match(&cols(sa, Some(p)), &cols(sb, Some(p)), tolerance))
                .sum()
        } else {
            greedy_match(&cols(sa, None), &cols(sb, None), tolerance)
        };
        let n_both = self
            .sites
            .iter()
            .filter(|s| s.observable.test(a) && s.observable.test(b))
            .count() as u32;
        self.make_pair(a, b, shared, n_both)
    }

    pub fn summary(&self) -> MatrixSummary {
        MatrixSummary {
            n_sites: self.sites.len() as u64,
            n_taxa: self.taxa.len() as u64,
            n_annotated: self.n_annotated as u64,
            n_introns_projected: self.taxon_sites.iter().map(|v| v.len() as u64).sum(),
            singletons: self.singleton_sites().len() as u64,
            ubiquitous: self.ubiquitous_sites().len() as u64,
            mean_per_annotated_taxon: self.mean_sites_per_annotated_taxon(),
            alen: self.alen,
        }
    }

    /// Compact human-readable one-liner.
    pub fn describe(&self) -> String {
        let s = self.summary();
        format!(
            "{} sites x {} taxa ({} annotated), {} projected introns, mean {:.1}/annotated taxon",
            s.n_sites, s.n_taxa, s.n_annotated, s.n_introns_projected, s.mean_per_annotated_taxon
        )
    }
}

/// Run-level totals for a matrix.
#[derive(Clone, Copy, Debug, Default)]
pub struct MatrixSummary {
    pub n_sites: u64,
    pub n_taxa: u64,
    pub n_annotated: u64,
    pub n_introns_projected: u64,
    pub singletons: u64,
    pub ubiquitous: u64,
    pub mean_per_annotated_taxon: f64,
    pub alen: usize,
}

/// One pairwise comparison of two taxa's intron positions.
#[derive(Clone, Debug)]
pub struct PairStat {
    pub a: usize,
    pub b: usize,
    pub taxon_a: String,
    pub taxon_b: String,
    pub n_a: u32,
    pub n_b: u32,
    pub shared: u32,
    /// `shared / (n_a + n_b - shared)`.
    pub jaccard: f64,
    /// `shared / min(n_a, n_b)` — the usual "fraction of conserved positions" figure.
    pub fraction_of_smaller: f64,
    /// `false` when either taxon has no annotation: the comparison is missing data, not zero.
    pub both_annotated: bool,
    /// Sites where *both* taxa were observed (not `?`). `shared`, `jaccard` and `fraction_of_smaller`
    /// are only comparable between pairs with a similar overlap, so read them together.
    pub n_both_observable: u32,
}

/// Sites as a NEXUS discrete-character block, for external tree programs.
pub fn characters_nexus(m: &SiteMatrix) -> String {
    let mut s = String::new();
    s.push_str("#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=");
    let _ = write!(
        s,
        "{} NCHAR={};\n  FORMAT SYMBOLS=\"0 1\" MISSING=?;\n  MATRIX\n",
        m.taxa.len(),
        m.sites.len()
    );
    for (t, name) in m.taxa.iter().enumerate() {
        let _ = write!(s, "{name}\t");
        for site in &m.sites {
            s.push(site.cell(t).symbol());
        }
        s.push('\n');
    }
    s.push_str("  ;\n;\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::ProteinSite;
    use crate::project::{project_into_alignment, ColumnSite, TaxonProjection};

    /// Build a matrix directly from `(taxon, column, phase)` triples.
    fn matrix(spec: &[(&str, u32, Phase)], taxa: &[&str], annotated: &[bool]) -> SiteMatrix {
        let mut ap = AlignmentProjection {
            alignment: "t".into(),
            alen: 100,
            ..Default::default()
        };
        for (t, id) in taxa.iter().enumerate() {
            ap.taxa.push(TaxonProjection {
                taxon: t,
                row_id: id.to_string(),
                taxon_name: id.to_string(),
                protein_idx: Some(t),
                sites: Vec::new(),
                report: Default::default(),
                residues_aligned: 0,
                has_annotation: annotated[t],
                unannotated_reason: None,
                in_sequence_set: annotated[t],
                zero_introns: false,
            });
        }
        for (id, col, phase) in spec {
            let t = taxa.iter().position(|x| x == id).unwrap();
            ap.sites.push(ColumnSite {
                taxon: t as u32,
                col: *col,
                phase: *phase,
                site: ProteinSite::from_cds_offset(col * 3 + *phase as u32),
                clamped: false,
            });
        }
        ap.sites.sort_unstable_by_key(|c| (c.col, c.phase, c.taxon));
        SiteMatrix::build(&ap)
    }

    #[test]
    fn bitset_basics() {
        let mut b = Bitset::new(130);
        assert!(b.set(0) && b.set(64) && b.set(129));
        assert!(!b.set(64), "already set");
        assert_eq!(b.count_ones(), 3);
        assert!(b.test(64) && !b.test(65));
        assert_eq!(b.iter_ones().collect::<Vec<_>>(), vec![0, 64, 129]);
        let mut c = Bitset::new(130);
        c.set(64);
        c.set(129);
        assert_eq!(b.and_count(&c), 2);
        assert_eq!(Bitset::new(0).iter_ones().count(), 0);
    }

    #[test]
    fn collapses_duplicate_sites_and_keeps_order() {
        let m = matrix(
            &[("A", 4, 0), ("B", 4, 0), ("A", 6, 1)],
            &["A", "B"],
            &[true, true],
        );
        assert_eq!(m.n_sites(), 2);
        assert_eq!(m.sites[0].label(), "4:0");
        assert_eq!(m.sites[0].n_taxa(), 2);
        assert_eq!(m.sites[1].n_taxa(), 1);
        assert_eq!(m.pattern(0), "11");
        assert_eq!(m.pattern(1), "10");
        assert_eq!(m.per_taxon_counts(), vec![2, 1]);
        assert_eq!(m.alen, 100);
    }

    #[test]
    fn same_column_different_phase_is_two_sites() {
        // (9,0) and (9,1) are one codon apart but distinct positions, so two sites — not one.
        let m = matrix(
            &[("A", 9, 0), ("A", 9, 1), ("B", 9, 0)],
            &["A", "B"],
            &[true, true],
        );
        assert_eq!(m.n_sites(), 2);
        assert_eq!(m.taxa_names_at(0), vec!["A", "B"]);
        assert_eq!(m.taxa_names_at(1), vec!["A"]);
        assert_eq!(m.per_taxon_counts(), vec![2, 1]);
        assert_eq!(
            m.pair_stat(0, 1, 0, true).shared,
            1,
            "only the same column *and* phase matches"
        );
        assert_eq!(
            m.pair_stat(0, 1, 0, false).shared,
            1,
            "position-only matching still finds one"
        );
    }

    #[test]
    fn unannotated_taxa_are_not_counted_as_absent_evidence() {
        let m = matrix(&[("A", 4, 0)], &["A", "B", "C"], &[true, false, false]);
        assert_eq!(m.n_annotated, 1);
        assert_eq!(
            m.ubiquitous_sites(),
            vec![0],
            "the only annotated taxon has it, so it is universal"
        );
        assert_eq!(m.summary().n_annotated, 1);
        assert_eq!(m.mean_sites_per_annotated_taxon(), 1.0);
    }

    #[test]
    fn spectrum_and_singletons() {
        let m = matrix(
            &[
                ("A", 1, 0),
                ("B", 1, 0),
                ("C", 1, 0),
                ("A", 2, 0),
                ("B", 2, 0),
                ("A", 3, 1),
            ],
            &["A", "B", "C"],
            &[true, true, true],
        );
        let h = m.conservation_spectrum();
        assert_eq!((h[1], h[2], h[3]), (1, 1, 1), "one site per sharing level");
        assert_eq!(m.singleton_sites().len(), 1);
        assert_eq!(m.ubiquitous_sites(), vec![0]);
        assert_eq!(m.summary().n_sites, 3);
    }

    #[test]
    fn pairwise_exact_matching_and_jaccard() {
        let m = matrix(
            &[
                ("A", 4, 0),
                ("A", 8, 1),
                ("B", 4, 0),
                ("B", 9, 1),
                ("C", 4, 0),
                ("C", 8, 1),
            ],
            &["A", "B", "C"],
            &[true, true, true],
        );
        let ab = m.pair_stat(0, 1, 0, true);
        assert_eq!((ab.n_a, ab.n_b, ab.shared), (2, 2, 1));
        assert!((ab.jaccard - 1.0 / 3.0).abs() < 1e-12, "{ab:?}");
        assert!((ab.fraction_of_smaller - 0.5).abs() < 1e-12);
        let ac = m.pair_stat(0, 2, 0, true);
        assert_eq!(ac.shared, 2);
        assert!(
            (ac.jaccard - 1.0).abs() < 1e-12,
            "identical sets give Jaccard 1"
        );
    }

    #[test]
    fn tolerance_absorbs_position_drift() {
        let m = matrix(
            &[("A", 10, 0), ("B", 11, 0), ("C", 20, 0)],
            &["A", "B", "C"],
            &[true, true, true],
        );
        assert_eq!(
            m.pair_stat(0, 1, 0, true).shared,
            0,
            "one column apart is not the same site"
        );
        assert_eq!(
            m.pair_stat(0, 1, 1, true).shared,
            1,
            "±1 tolerance calls it shared"
        );
        assert_eq!(
            m.pair_stat(0, 2, 5, true).shared,
            0,
            "ten columns apart is beyond any drift"
        );
    }

    #[test]
    fn phase_is_respected_when_required() {
        let m = matrix(&[("A", 5, 0), ("B", 5, 2)], &["A", "B"], &[true, true]);
        assert_eq!(m.pair_stat(0, 1, 0, true).shared, 0);
        assert_eq!(
            m.pair_stat(0, 1, 0, false).shared,
            1,
            "position-only matching ignores phase"
        );
    }

    #[test]
    fn each_site_supports_at_most_one_shared_call() {
        // A's single site must not be counted twice against two nearby B sites.
        let m = matrix(
            &[("A", 10, 0), ("B", 10, 0), ("B", 11, 0)],
            &["A", "B"],
            &[true, true],
        );
        assert_eq!(m.pair_stat(0, 1, 1, true).shared, 1);
        assert_eq!(m.pair_stat(0, 1, 0, true).shared, 1);
    }

    #[test]
    fn dense_matrix_is_symmetric_with_self_counts_on_the_diagonal() {
        let m = matrix(
            &[("A", 1, 0), ("B", 1, 0), ("C", 2, 0)],
            &["A", "B", "C"],
            &[true, true, true],
        );
        let d = m.pairwise_matrix(0, true);
        let at = |i: usize, j: usize| d[i * 3 + j];
        assert_eq!(at(0, 1), 1);
        assert_eq!(at(1, 0), 1, "symmetric");
        assert_eq!(at(0, 0), 1, "diagonal is the taxon's own count");
        assert_eq!(at(2, 2), 1);
        assert_eq!(at(0, 2), 0);
        assert_eq!(m.pairwise(0, true).len(), 3, "upper triangle only");
    }

    #[test]
    fn empty_matrix_has_no_sites() {
        let m = matrix(&[], &["A", "B"], &[true, true]);
        assert_eq!(m.n_sites(), 0);
        assert_eq!(m.summary().mean_per_annotated_taxon, 0.0);
        assert!(m.ubiquitous_sites().is_empty());
    }

    #[test]
    fn nexus_block_is_a_valid_character_matrix() {
        let m = matrix(
            &[("A", 1, 0), ("B", 1, 0), ("B", 2, 1)],
            &["A", "B"],
            &[true, true],
        );
        let s = characters_nexus(&m);
        assert!(s.starts_with("#NEXUS"), "{s}");
        assert!(s.contains("NTAX=2 NCHAR=2"), "{s}");
        assert!(s.contains("A\t10"), "A lacks the second site\n{s}");
        assert!(s.contains("B\t11"), "{s}");
    }

    #[test]
    fn end_to_end_from_annotations_through_matrix() {
        use crate::model::{AnnotatedIntron, Strand, TranscriptAnnotation};
        use crate::msa::Alignment;
        use crate::project::{join_annotations, JoinOptions};
        use crate::seqio::ProteinSet;
        use std::io::Cursor;

        let mut ps = ProteinSet::new();
        ps.read_one(
            Cursor::new(b">A\nMKVLGAWEQP\n>B\nAWEQP\n"),
            std::path::Path::new("m"),
        )
        .unwrap();
        let txs = vec![
            TranscriptAnnotation {
                transcript_id: "A.t".into(),
                gene_id: None,
                protein_id: Some("A".into()),
                seqid: "c".into(),
                strand: Strand::Plus,
                cds_len: 33,
                cds_segments: 2,
                introns: vec![AnnotatedIntron {
                    start: 1,
                    end: 2,
                    rank: 0,
                    cds_offset: Some(18),
                }],
                frame_mismatches: 0,
                source: "t".into(),
            },
            TranscriptAnnotation {
                transcript_id: "B.t".into(),
                gene_id: None,
                protein_id: Some("B".into()),
                seqid: "c".into(),
                strand: Strand::Plus,
                cds_len: 18,
                cds_segments: 2,
                introns: vec![AnnotatedIntron {
                    start: 3,
                    end: 4,
                    rank: 0,
                    cds_offset: Some(3),
                }],
                frame_mismatches: 0,
                source: "t".into(),
            },
        ];
        let proj = join_annotations(&txs, &ps, &JoinOptions::default()).unwrap();
        let a = Alignment::from_fasta_str(">A\nMKVLGAWEQP\n>B\n-----AWEQP\n", "fam").unwrap();
        let ap = project_into_alignment(&a, &proj, &ps, &Default::default());
        let m = SiteMatrix::build(&ap);
        assert_eq!(
            m.n_sites(),
            1,
            "the two species' introns land on the same column"
        );
        assert_eq!(m.pattern(0), "11");
        assert_eq!(m.pair_stat(0, 1, 0, true).shared, 1);
        assert!(
            m.describe().starts_with("1 sites x 2 taxa"),
            "{}",
            m.describe()
        );
    }

    // ---- proposed: matrix.tsv and MALIN round trips on random matrices ----
    #[test]
    fn random_matrices_round_trip_through_matrix_tsv_and_malin() {
        let mut s: u64 = 0xDEAD_BEEF_CAFE_F00D;
        let mut rnd = move |n: u64| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s % n
        };
        let dir = std::env::temp_dir().join(format!("ip_rt_{}", std::process::id()));
        for case in 0..200 {
            let n_taxa = 1 + rnd(7) as usize;
            let taxa: Vec<String> = (0..n_taxa).map(|i| format!("t{i}")).collect();
            let mut used = std::collections::BTreeSet::new();
            let mut sites = Vec::new();
            for _ in 0..rnd(12) {
                let (col, phase) = (rnd(40) as u32, rnd(3) as u8);
                if !used.insert((col, phase)) {
                    continue;
                }
                let mut presence = Bitset::new(n_taxa);
                let mut observable = Bitset::new(n_taxa);
                for t in 0..n_taxa {
                    match rnd(3) {
                        0 => {
                            presence.set(t);
                            observable.set(t);
                        }
                        1 => {
                            observable.set(t);
                        }
                        _ => {}
                    }
                }
                sites.push(Site {
                    col,
                    phase,
                    presence,
                    observable,
                });
            }
            sites.sort_by_key(|s| (s.col, s.phase));
            let m = SiteMatrix::from_parts(
                &format!("aln{case}"),
                taxa,
                sites,
                crate::project::PhaseMode::Strict,
            );
            let path = dir.join(format!("m{case}.tsv"));
            crate::tables::matrix_table(&m).write(&path).unwrap();
            let back = SiteMatrix::read_tsv(&path).unwrap();
            assert_eq!(back.taxa, m.taxa, "case {case}");
            assert_eq!(back.annotated, m.annotated, "case {case}");
            assert_eq!(back.n_sites(), m.n_sites(), "case {case}");
            for i in 0..m.n_sites() {
                assert_eq!(
                    (back.sites[i].col, back.sites[i].phase),
                    (m.sites[i].col, m.sites[i].phase)
                );
                assert_eq!(back.pattern(i), m.pattern(i), "case {case} site {i}");
            }
            assert_eq!(back.taxon_sites, m.taxon_sites, "case {case}");
            // MALIN keeps the same cells; a matrix without sites has nothing to export.
            if m.n_sites() == 0 {
                assert!(crate::malin::to_malin(&[&m]).is_err(), "case {case}");
                continue;
            }
            let text = crate::malin::to_malin(&[&m]).unwrap();
            let mb = crate::malin::from_malin(&text, "x")
                .unwrap_or_else(|e| panic!("case {case}: {e}\n{text}"));
            for i in 0..m.n_sites() {
                assert_eq!(mb[0].pattern(i), m.pattern(i), "case {case} site {i}");
            }
        }
    }

    #[test]
    fn parallel_pairwise_equals_the_per_pair_reference_on_random_matrices() {
        let mut state: u64 = 0xA5A5_1234_9876_5432;
        let mut rnd = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        for case in 0..60 {
            let n_taxa = 2 + rnd(9) as usize;
            let taxa: Vec<String> = (0..n_taxa).map(|i| format!("t{i}")).collect();
            let mut used = std::collections::BTreeSet::new();
            let mut sites = Vec::new();
            for _ in 0..rnd(40) {
                let (col, phase) = (rnd(25) as u32, rnd(3) as u8);
                if !used.insert((col, phase)) {
                    continue;
                }
                let (mut presence, mut observable) = (Bitset::new(n_taxa), Bitset::new(n_taxa));
                for t in 0..n_taxa {
                    match rnd(4) {
                        0 | 1 => {
                            presence.set(t);
                            observable.set(t);
                        }
                        2 => {
                            observable.set(t);
                        }
                        _ => {}
                    }
                }
                sites.push(Site {
                    col,
                    phase,
                    presence,
                    observable,
                });
            }
            sites.sort_by_key(|s| (s.col, s.phase));
            let m = SiteMatrix::from_parts("x", taxa, sites, crate::project::PhaseMode::Strict);
            for (tol, phase) in [(0u32, true), (0, false), (2, true), (1, false)] {
                let fast = m.pairwise(tol, phase);
                let mut slow = Vec::new();
                for a in 0..n_taxa {
                    for b in (a + 1)..n_taxa {
                        slow.push(m.pair_stat(a, b, tol, phase));
                    }
                }
                assert_eq!(fast.len(), slow.len(), "case {case}");
                for (f, g) in fast.iter().zip(&slow) {
                    assert_eq!(
                        (f.a, f.b, f.shared, f.n_a, f.n_b, f.n_both_observable),
                        (g.a, g.b, g.shared, g.n_a, g.n_b, g.n_both_observable),
                        "case {case} tol {tol} phase {phase}"
                    );
                    assert_eq!(f.jaccard.to_bits(), g.jaccard.to_bits());
                    assert_eq!(
                        f.fraction_of_smaller.to_bits(),
                        g.fraction_of_smaller.to_bits()
                    );
                }
            }
        }
    }
}
