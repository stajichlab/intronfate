//! Joining annotation to protein sequences, then projecting intron positions into alignment
//! columns.
//!
//! Two decisions live here, and both are the kind that quietly corrupt results if left implicit:
//!
//! **Which protein does this transcript describe?** Annotation and sequence databases name the same
//! molecule differently. [`join_annotations`] resolves the annotation's protein id against the
//! sequence set (exact, then version-insensitive), reports both sides of every failure, and
//! cross-checks CDS length against protein length — because a join that lands on the wrong paralog
//! produces *plausible-looking* tables, which is worse than no tables at all.
//!
//! **Where does the position land in the alignment?** `residue → column` is a lookup, but a site can
//! anchor past the residues the alignment row actually covers, and several transcripts can describe
//! one protein. Both are explicit policies on [`ProjectOptions`], not silent behaviour.

use rustc_hash::FxHashMap;

use crate::error::Result;
use crate::matrix::Bitset;
use crate::model::{Phase, ProteinSite, Strand, TranscriptAnnotation};
use crate::msa::{Alignment, AlignmentRow};
use crate::report::Diagnostics;
use crate::seqio::ProteinSet;

/// How strictly a site must agree in codon phase with the site it is matched against.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum PhaseMode {
    /// Position **and** phase must agree: two introns one nucleotide apart inside the same codon are
    /// different sites. This is the standard comparative-genomics notion of a shared position.
    #[default]
    Strict,
    /// Only the codon matters; phase is ignored when grouping sites.
    Ignore,
}

/// What to do with an intron anchored to a codon the alignment row does not cover.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum TailPolicy {
    /// Leave the site out and count it. Correct when the alignment really lacks that residue.
    #[default]
    Drop,
    /// Snap to the row's last column. Useful for trim-to-core alignments whose edges were cut, at
    /// the cost of piling unrelated C-terminal sites onto one column.
    Clamp,
}

/// Controls for [`join_annotations`] (annotation ↔ sequence).
#[derive(Clone, Debug)]
pub struct JoinOptions {
    /// Flag proteins whose length disagrees with the CDS-derived expectation by more than this many
    /// residues. Annotations differ on whether the stop codon is part of the CDS, so the default of
    /// 3 absorbs that convention without hiding a genuine mis-join.
    pub length_slack: u64,
    /// Note when transcripts encoding one protein disagree about its intron positions.
    pub check_isoform_consistency: bool,
    /// List protein sequences that no annotation record joined to.
    pub report_unmatched_proteins: bool,
}

impl Default for JoinOptions {
    fn default() -> Self {
        JoinOptions {
            length_slack: 3,
            check_isoform_consistency: true,
            report_unmatched_proteins: true,
        }
    }
}

/// Controls for projection into alignment columns.
#[derive(Clone, Debug)]
pub struct ProjectOptions {
    pub phase_mode: PhaseMode,
    pub tail: TailPolicy,
    /// Skip sites whose stored phase and offset disagree (only reachable from hand-built inputs).
    pub skip_incoherent: bool,
    /// Flanking residues a taxon must have on each side of a site before its absence counts as an
    /// observation (see [`cell_observable`]); otherwise the cell is missing (`?`), not `0`.
    /// `0` means "gap at the anchoring column only". Phase-0 introns always need one residue on each
    /// side, as MALIN requires.
    pub flank: u32,
    /// Alignment row id -> taxon name. Real runs align one protein per species per gene family, so the
    /// row ids differ in every alignment while the tree's tips are species; this maps one to the other.
    /// Rows without an entry keep their own id.
    pub taxon_names: FxHashMap<String, String>,
}

impl Default for ProjectOptions {
    fn default() -> Self {
        ProjectOptions {
            phase_mode: PhaseMode::Strict,
            tail: TailPolicy::Drop,
            skip_incoherent: false,
            flank: 1,
            taxon_names: FxHashMap::default(),
        }
    }
}

/// One intron, tied to the protein it was projected onto.
#[derive(Clone, Debug)]
pub struct SiteRow {
    /// Index into [`Projection::tx_meta`]: the transcript's id, gene, sequence and strand live there once
    /// per transcript, not once per intron. Read them through the `Projection` accessors.
    pub tx: u32,
    /// Index into the [`ProteinSet`]; the protein's id is `Projection::protein_id(&row)`.
    pub protein_idx: usize,
    /// 0-based rank among all introns of the transcript, UTR introns included.
    pub intron_rank: u32,
    /// 0-based rank among the transcript's *coding* introns.
    pub coding_rank: u32,
    pub genomic_start: u64,
    pub genomic_end: u64,
    pub intron_len: u64,
    pub site: ProteinSite,
    pub protein_len: usize,
    /// `false` when the residue the site anchors to is beyond the protein's end.
    pub in_range: bool,
}

/// What a [`SiteRow`] needs to say about its transcript, stored once per transcript.
#[derive(Clone, Debug)]
pub struct TxMeta {
    pub transcript_id: String,
    pub gene_id: Option<String>,
    pub seqid: String,
    pub strand: Strand,
}

/// An intron position that failed to join, with the reason.
#[derive(Clone, Debug)]
pub struct UnmatchedRecord {
    pub transcript_id: String,
    pub gene_id: Option<String>,
    pub attempted_id: String,
    pub reason: &'static str,
    pub source: String,
}

/// The annotation↔sequence view of a gene family: every projectable intron of every protein.
#[derive(Clone, Debug, Default)]
pub struct Projection {
    pub rows: Vec<SiteRow>,
    /// One entry per transcript that produced at least one row; indexed by [`SiteRow::tx`].
    pub tx_meta: Vec<TxMeta>,
    /// Id of every protein that produced a row, indexed by protein index (rows hold only the index).
    protein_names: Vec<Option<Box<str>>>,
    /// Per protein index: where its rows sit in `rows` (`start`, `len`). Valid once the join is finished,
    /// when rows are grouped by protein. Proteins are indexed densely, so this is a plain array rather
    /// than a map with one entry (and one set) per protein.
    row_ranges: Vec<(u32, u32)>,
    /// How many transcripts joined to each protein, indexed by protein index. A protein with a count above
    /// zero is **annotated, even with zero coding introns**: a gene whose introns have all been lost is an
    /// observation of absence, while a gene with no annotation at all is no observation. Collapsing the
    /// two — which an intron-count test cannot do — discards precisely the negative evidence
    /// intron-loss analyses are built on.
    pub n_transcripts_of: Vec<u32>,
    pub unmatched: Vec<UnmatchedRecord>,
    pub proteins_without_annotation: Vec<String>,
    pub diagnostics: Diagnostics,
    pub joined_proteins: usize,
    pub total_introns: u64,
    pub utr_introns: u64,
}

impl Projection {
    /// Transcript metadata of a row.
    pub fn tx(&self, r: &SiteRow) -> &TxMeta {
        &self.tx_meta[r.tx as usize]
    }

    pub fn transcript_id(&self, r: &SiteRow) -> &str {
        &self.tx(r).transcript_id
    }

    pub fn gene_id(&self, r: &SiteRow) -> Option<&str> {
        self.tx(r).gene_id.as_deref()
    }

    pub fn seqid(&self, r: &SiteRow) -> &str {
        &self.tx(r).seqid
    }

    pub fn strand(&self, r: &SiteRow) -> Strand {
        self.tx(r).strand
    }

    /// The id of the protein a row was projected onto.
    pub fn protein_id(&self, r: &SiteRow) -> &str {
        self.protein_names
            .get(r.protein_idx)
            .and_then(|n| n.as_deref())
            .unwrap_or("")
    }

    /// The distinct protein-space sites of a protein, ascending; `None` when it has no coding intron.
    pub fn sites_of(&self, protein_idx: usize) -> Option<Vec<ProteinSite>> {
        let &(start, len) = self.row_ranges.get(protein_idx)?;
        if len == 0 {
            return None;
        }
        let mut v: Vec<ProteinSite> = self.rows[start as usize..(start + len) as usize]
            .iter()
            .map(|r| r.site)
            .collect();
        v.sort_unstable();
        v.dedup();
        Some(v)
    }

    /// `true` when some transcript joined to this protein, whatever its intron count.
    pub fn is_annotated(&self, protein_idx: usize) -> bool {
        self.n_transcripts_of
            .get(protein_idx)
            .is_some_and(|&n| n > 0)
    }

    /// Proteins with at least one coding intron.
    pub fn protein_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.row_ranges
            .iter()
            .enumerate()
            .filter(|(_, r)| r.1 > 0)
            .map(|(i, _)| i)
    }

    /// Every protein with annotation, including those with zero coding introns. Use this — not
    /// [`protein_indices`](Self::protein_indices) — whenever the question is "did we look?".
    pub fn annotated_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.n_transcripts_of
            .iter()
            .enumerate()
            .filter(|(_, &n)| n > 0)
            .map(|(i, _)| i)
    }

    /// Total number of distinct protein-space positions across all joined proteins.
    pub fn n_distinct_sites(&self) -> usize {
        self.protein_indices()
            .map(|i| self.sites_of(i).map_or(0, |s| s.len()))
            .sum()
    }
}

/// Join transcripts to proteins and collect every coding intron in protein space.
pub fn join_annotations(
    transcripts: &[TranscriptAnnotation],
    proteins: &ProteinSet,
    opts: &JoinOptions,
) -> Result<Projection> {
    let mut joiner = Joiner::new(proteins, opts);
    for tx in transcripts {
        joiner.add(tx);
    }
    joiner.finish()
}

/// Incremental form of [`join_annotations`]: transcripts are added one at a time, so a caller that
/// streams annotation (see [`crate::annotation::AnnotationSink`]) never has to hold every transcript of
/// every genome in memory — only what the join keeps.
#[derive(Debug)]
pub struct Joiner<'a> {
    p: Projection,
    proteins: &'a ProteinSet,
    opts: &'a JoinOptions,
}

impl<'a> Joiner<'a> {
    pub fn new(proteins: &'a ProteinSet, opts: &'a JoinOptions) -> Self {
        Joiner {
            p: Projection {
                diagnostics: Diagnostics::new(500),
                protein_names: vec![None; proteins.len()],
                row_ranges: vec![(0, 0); proteins.len()],
                n_transcripts_of: vec![0; proteins.len()],
                ..Default::default()
            },
            proteins,
            opts,
        }
    }

    /// Join one transcript to its protein and record its coding introns.
    pub fn add(&mut self, tx: &TranscriptAnnotation) {
        let p = &mut self.p;
        let proteins = self.proteins;
        let opts = self.opts;
        p.total_introns += tx.introns.len() as u64;
        p.utr_introns += tx.introns.iter().filter(|i| !i.is_coding()).count() as u64;

        let Some(pid) = tx.protein_id.as_deref() else {
            p.diagnostics.warn(
                "no_protein_id",
                &tx.transcript_id,
                "transcript carries no protein identifier",
            );
            p.unmatched.push(UnmatchedRecord {
                transcript_id: tx.transcript_id.clone(),
                gene_id: tx.gene_id.clone(),
                attempted_id: String::new(),
                reason: "no_protein_id",
                source: tx.source.clone(),
            });
            return;
        };
        let idx = match proteins.lookup(pid) {
            Some(i) => i,
            None => {
                p.diagnostics.warn(
                    "protein_not_in_fasta",
                    pid,
                    format!(
                        "transcript {:?} names protein {pid:?}, which is absent from the sequence set",
                        tx.transcript_id
                    ),
                );
                p.unmatched.push(UnmatchedRecord {
                    transcript_id: tx.transcript_id.clone(),
                    gene_id: tx.gene_id.clone(),
                    attempted_id: pid.to_string(),
                    reason: "protein_not_in_fasta",
                    source: tx.source.clone(),
                });
                return;
            }
        };
        let seq = &proteins.seqs[idx];

        if tx.cds_len == 0 {
            p.diagnostics.warn(
                "zero_length_cds",
                pid,
                "joined protein has a transcript with no coding length",
            );
        }

        // Length cross-check: the cheapest signal that a join landed on the wrong molecule.
        if opts.length_slack < u64::MAX && tx.cds_len > 0 {
            let expect_lo = tx.expected_protein_len_min();
            let expect_hi = tx.expected_protein_len_max();
            let got = seq.len as u64;
            if got + opts.length_slack < expect_lo || got > expect_hi + opts.length_slack {
                p.diagnostics.warn(
                    "length_mismatch",
                    pid,
                    format!(
                        "protein is {got} residues but a {} nt CDS implies {expect_lo}..{expect_hi} (CDS of {:?})",
                        tx.cds_len, tx.transcript_id
                    ),
                );
            }
        }

        p.n_transcripts_of[idx] += 1;

        let mut coding_rank = 0u32;
        let mut tx_index: Option<u32> = None;
        for intron in &tx.introns {
            let Some(off) = intron.cds_offset else {
                continue;
            };
            let site = ProteinSite::from_cds_offset(off);
            let in_range = (site.codon_index0() as usize) < seq.len;
            if !in_range {
                p.diagnostics.warn(
                    "site_beyond_protein_end",
                    pid,
                    format!(
                        "intron rank {} at codon {} lies past the {} residues of {pid:?}",
                        intron.rank,
                        site.codon(),
                        seq.len
                    ),
                );
            }
            let tx_i = *tx_index.get_or_insert_with(|| {
                p.tx_meta.push(TxMeta {
                    transcript_id: tx.transcript_id.clone(),
                    gene_id: tx.gene_id.clone(),
                    seqid: tx.seqid.clone(),
                    strand: tx.strand,
                });
                (p.tx_meta.len() - 1) as u32
            });
            p.protein_names[idx].get_or_insert_with(|| seq.id.as_str().into());
            p.rows.push(SiteRow {
                tx: tx_i,
                protein_idx: idx,
                intron_rank: intron.rank,
                coding_rank,
                genomic_start: intron.start,
                genomic_end: intron.end,
                intron_len: intron.length(),
                site,
                protein_len: seq.len,
                in_range,
            });
            coding_rank += 1;
        }
    }

    /// Finish: protein-level checks, isoform consistency, deterministic row order.
    pub fn finish(self) -> Result<Projection> {
        let Joiner {
            mut p,
            proteins,
            opts,
        } = self;
        p.joined_proteins = p.n_transcripts_of.iter().filter(|&&n| n > 0).count();

        if opts.report_unmatched_proteins {
            for (i, seq) in proteins.seqs.iter().enumerate() {
                if p.n_transcripts_of[i] == 0 {
                    p.proteins_without_annotation.push(seq.id.clone());
                }
            }
            if !p.proteins_without_annotation.is_empty() {
                p.diagnostics.warn(
                    "protein_without_annotation",
                    "sequence set",
                    format!(
                        "{} protein(s) had no matching annotation record",
                        p.proteins_without_annotation.len()
                    ),
                );
            }
        }

        // Deterministic row order: by protein id, then transcript id, then coding rank.
        let mut rows = std::mem::take(&mut p.rows);
        {
            let (names, meta) = (&p.protein_names, &p.tx_meta);
            let key = |r: &SiteRow| {
                (
                    names[r.protein_idx].as_deref().unwrap_or(""),
                    meta[r.tx as usize].transcript_id.as_str(),
                    r.coding_rank,
                )
            };
            rows.sort_by(|a, b| key(a).cmp(&key(b)));
        }
        p.rows = rows;

        // Rows of one protein are now contiguous: index them, then check isoforms group by group.
        let mut start = 0usize;
        while start < p.rows.len() {
            let idx = p.rows[start].protein_idx;
            let mut end = start + 1;
            while end < p.rows.len() && p.rows[end].protein_idx == idx {
                end += 1;
            }
            p.row_ranges[idx] = (start as u32, (end - start) as u32);
            start = end;
        }

        if opts.check_isoform_consistency {
            // Only a protein with several transcripts can have conflicting isoforms.
            let mut idxs: Vec<usize> = (0..p.row_ranges.len())
                .filter(|&i| p.n_transcripts_of[i] > 1 && p.row_ranges[i].1 > 0)
                .collect();
            idxs.sort_unstable();
            for idx in idxs {
                let (s0, len) = p.row_ranges[idx];
                let group = &p.rows[s0 as usize..(s0 + len) as usize];
                let mut per_tx: FxHashMap<u32, usize> = FxHashMap::default();
                for r in group {
                    *per_tx.entry(r.tx).or_default() += 1;
                }
                let mut distinct: Vec<ProteinSite> = group.iter().map(|r| r.site).collect();
                distinct.sort_unstable();
                distinct.dedup();
                // Any transcript whose set is smaller than the union differs from another isoform.
                if per_tx.len() > 1 && per_tx.values().any(|&n| n != distinct.len()) {
                    p.diagnostics.info(
                        "isoform_intron_conflict",
                        proteins.id(idx).unwrap_or("?"),
                        format!(
                            "{} transcripts encode this protein with different intron sets ({} distinct positions)",
                            per_tx.len(),
                            distinct.len()
                        ),
                    );
                }
            }
        }
        Ok(p)
    }
}

/// An [`AnnotationSink`](crate::annotation::AnnotationSink) that joins each transcript to its protein the
/// moment it arrives and keeps only what the join keeps. Reading annotation straight into this, rather than
/// collecting a `Vec<TranscriptAnnotation>` first, is what keeps memory flat in the number of genomes: the
/// per-transcript structures of a million-transcript annotation are several times the size of its introns.
///
/// Annotation-level diagnostics (frame mismatches, dropped transcripts, id collisions, ...) are collected
/// in [`diagnostics`](Self::diagnostics) so the caller can report them next to the join's own.
#[derive(Debug)]
pub struct JoiningSink<'a> {
    joiner: Joiner<'a>,
    pub diagnostics: Diagnostics,
    /// Transcripts received (kept ones), whether or not they joined.
    pub n_transcripts: usize,
}

impl<'a> JoiningSink<'a> {
    pub fn new(proteins: &'a ProteinSet, opts: &'a JoinOptions) -> Self {
        JoiningSink {
            joiner: Joiner::new(proteins, opts),
            diagnostics: Diagnostics::new(500),
            n_transcripts: 0,
        }
    }

    /// The projection, the annotation diagnostics, and the number of transcripts seen.
    pub fn finish(self) -> Result<(Projection, Diagnostics, usize)> {
        Ok((self.joiner.finish()?, self.diagnostics, self.n_transcripts))
    }
}

impl crate::annotation::AnnotationSink for JoiningSink<'_> {
    fn push_transcript(&mut self, tx: TranscriptAnnotation) -> Result<()> {
        self.n_transcripts += 1;
        self.joiner.add(&tx);
        Ok(())
    }

    fn diagnostic(
        &mut self,
        level: crate::report::Level,
        kind: &'static str,
        subject: String,
        message: String,
    ) {
        self.diagnostics.push(level, kind, subject, message);
    }
}

/// What one alignment row contributed: its placed sites, and how many could not be placed.
#[derive(Clone, Debug)]
pub struct TaxonProjection {
    /// Index into [`Alignment::rows`].
    pub taxon: usize,
    /// The alignment row's id (the protein id the annotation is joined on).
    pub row_id: String,
    /// Taxon name used in matrices and tables: `row_id` unless `ProjectOptions::taxon_names` renames it.
    pub taxon_name: String,
    /// Index into the [`ProteinSet`], when the row joined to a sequence.
    pub protein_idx: Option<usize>,
    /// Sites in `(column, phase)` order.
    pub sites: Vec<ColumnSite>,
    pub report: ProjectionReport,
    /// Residues this row actually aligns (gaps and Stockholm inserts excluded).
    pub residues_aligned: usize,
    /// `false` when no transcript joined to this row's protein, so its absences carry no
    /// information. Note this is *not* the same as having zero introns: an annotated gene that lost
    /// every intron is annotated, and its zeros are evidence.
    pub has_annotation: bool,
    pub unannotated_reason: Option<&'static str>,
    /// `true` when the row's id names a sequence in the protein set (whether or not it is annotated).
    pub in_sequence_set: bool,
    /// The row is annotated but projected no intron site (the gene has none, or every one was dropped
    /// by the tail policy). Its absences are still observations.
    pub zero_introns: bool,
}

/// Outcome of projecting one protein's sites onto one alignment row.
#[derive(Clone, Copy, Debug, Default)]
pub struct ProjectionReport {
    /// Sites successfully placed on a column of their own.
    pub placed: usize,
    /// Sites whose anchor residue has no column in the row.
    pub tail_sites: usize,
    /// Of those, how many were snapped to the last column by [`TailPolicy::Clamp`].
    pub clamped: usize,
}

impl ProjectionReport {
    pub fn n_sites(&self) -> usize {
        self.placed + self.clamped
    }
}

/// One projected site: an intron position expressed as an alignment column.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ColumnSite {
    /// Taxon index within [`Alignment::rows`].
    pub taxon: u32,
    /// 0-based alignment column.
    pub col: u32,
    /// Phase used for grouping (always 0 under [`PhaseMode::Ignore`]).
    pub phase: Phase,
    /// The original protein-space position, kept for audit tables.
    pub site: ProteinSite,
    /// Placed on the row's last column by [`TailPolicy::Clamp`] although its residue is not there. Such
    /// a site is approximate by construction, so it is exempt from the observability window (whose
    /// right-hand flank would otherwise always fall past the end of the row and mask it).
    pub clamped: bool,
}

/// Every taxon's projected sites, ready to be turned into a matrix.
#[derive(Clone, Debug, Default)]
pub struct AlignmentProjection {
    pub alignment: String,
    /// Alignment length in columns, carried through for table headers.
    pub alen: usize,
    pub taxa: Vec<TaxonProjection>,
    /// All sites, sorted by `(column, phase, taxon)`.
    pub sites: Vec<ColumnSite>,
    pub coverage: Coverage,
    /// Recorded so downstream tables state how sites were grouped.
    pub phase_mode: PhaseMode,
    /// Per alignment column, which taxa have a **residue** there; `None` for hand-built projections
    /// that carry no alignment, which then fall back to "observable = annotated".
    ///
    /// This is what makes missing data per cell rather than per taxon: a taxon gapped around column
    /// *c* cannot show an intron at *c*, so its cell is `?` (not observable) rather than `0`
    /// (observed absent). Treating gapped columns as absence invents losses in exactly the regions
    /// where the alignment is least trustworthy.
    pub residue_at: Option<Vec<Bitset>>,
    /// Flank used when deciding observability ([`ProjectOptions::flank`]).
    pub flank: u32,
}

/// Window of columns a taxon must cover for a site at `col`: `phase` `Some(0)` is between residues,
/// `Some(1|2)` inside the residue at `col`, and `None` a site merged from several phases
/// (`--ignore-phase`), which needs the union of both windows so that no carrier's window is wider
/// than the window used for absences.
pub fn observability_window(col: u32, phase: Option<Phase>, flank: u32) -> (i64, i64) {
    let col = i64::from(col);
    let m = i64::from(flank.max(1));
    let n = i64::from(flank);
    match phase {
        Some(0) => (col - m, col + m - 1),
        Some(_) => (col - n, col + n),
        None => (col - m, col + n),
    }
}

/// Is `taxon`'s cell for a site at `(col, phase)` an observation (so `0` means absent) rather than
/// missing data?
///
/// MALIN's rule: absence is only informative when the taxon is aligned around the site. A phase-0
/// intron lies *between* residues, so it needs residues on both sides (columns `col-n..col` and
/// `col..col+n`, with `n = max(flank, 1)`); because those must be adjacent columns, a boundary
/// sitting beside an indel is missing rather than anchored to a different column per taxon. A phase
/// 1/2 intron lies inside the residue at `col` and needs `col-flank ..= col+flank`.
pub fn cell_observable(
    residue_at: &[Bitset],
    taxon: usize,
    col: u32,
    phase: Option<Phase>,
    flank: u32,
) -> bool {
    let has =
        |c: i64| -> bool { c >= 0 && residue_at.get(c as usize).is_some_and(|b| b.test(taxon)) };
    let (lo, hi) = observability_window(col, phase, flank);
    (lo..=hi).all(has)
}

/// Project every joined protein onto its alignment row.
pub fn project_into_alignment(
    aln: &Alignment,
    proj: &Projection,
    proteins: &ProteinSet,
    opts: &ProjectOptions,
) -> AlignmentProjection {
    let mut out = AlignmentProjection {
        alignment: aln.name.clone(),
        alen: aln.len,
        phase_mode: opts.phase_mode,
        taxa: Vec::with_capacity(aln.n_rows()),
        residue_at: Some(observability(aln)),
        flank: opts.flank,
        ..Default::default()
    };
    for (t, row) in aln.rows.iter().enumerate() {
        let pi = proteins.lookup(&row.id);
        let in_sequence_set = pi.is_some();
        // Annotation and "has introns" are different facts: a gene that lost every intron is
        // annotated, and its absences are evidence.
        let mut has_annotation = pi.is_some_and(|i| proj.is_annotated(i));
        // An alignment row must *be* the protein the annotation was joined to. A trimmed or edited
        // row (trimAl/Gblocks output, another isoform) would shift every projected position with no
        // sign of trouble, so such a row is not an observation at all.
        let mut row_mismatch = false;
        if let (true, Some(i)) = (has_annotation, pi) {
            if !row_matches_protein(row, proteins.residues(i).unwrap_or(&[])) {
                row_mismatch = true;
                has_annotation = false;
                out.coverage.row_mismatch.push(row.id.clone());
            }
        }
        let sites_in = pi.filter(|_| has_annotation).and_then(|i| proj.sites_of(i));
        let (sites, report) = match sites_in {
            Some(sites) => project_row(t, row, &sites, opts),
            None => (Vec::new(), ProjectionReport::default()),
        };
        let reason = if !in_sequence_set {
            Some("protein_not_in_sequence_set")
        } else if row_mismatch {
            Some("row_sequence_mismatch")
        } else if !has_annotation {
            Some("protein_not_annotated")
        } else {
            None
        };
        if has_annotation {
            out.coverage.annotated_taxa.push(t);
        } else {
            out.coverage.unannotated_taxa.push(row.id.clone());
        }
        // Computed before `sites` is moved into the literal below.
        let zero_introns = has_annotation && sites.is_empty();
        out.taxa.push(TaxonProjection {
            taxon: t,
            row_id: row.id.clone(),
            taxon_name: opts
                .taxon_names
                .get(&row.id)
                .cloned()
                .unwrap_or_else(|| row.id.clone()),
            protein_idx: pi,
            sites,
            report,
            residues_aligned: row.n_residues(),
            has_annotation,
            unannotated_reason: reason,
            in_sequence_set,
            zero_introns,
        });
    }
    out.sites = out
        .taxa
        .iter()
        .flat_map(|t| t.sites.iter().copied())
        .collect();
    out.sites
        .sort_unstable_by_key(|c| (c.col, c.phase, c.taxon));

    let present: FxHashMap<&str, ()> = aln.rows.iter().map(|r| (r.id.as_str(), ())).collect();
    // Counted, with a few examples: a run over many families and species would otherwise build one
    // string per annotated protein per alignment (tens of millions) to say something that is true of
    // almost every protein in almost every alignment.
    for pi in proj.annotated_indices() {
        if let Some(id) = proteins.id(pi) {
            if !present.contains_key(id) {
                out.coverage.n_absent_from_alignment += 1;
                let ex = &mut out.coverage.absent_examples;
                if ex.len() < ABSENT_EXAMPLES || ex.last().is_some_and(|l| id < l.as_str()) {
                    let at = ex.partition_point(|e| e.as_str() < id);
                    ex.insert(at, id.to_string());
                    ex.truncate(ABSENT_EXAMPLES);
                }
            }
        }
    }
    out
}

/// Which taxa have a residue in each alignment column.
pub fn observability(aln: &Alignment) -> Vec<Bitset> {
    let mut cols = vec![Bitset::new(aln.n_rows()); aln.len];
    for (t, row) in aln.rows.iter().enumerate() {
        for &c in &row.res2col {
            if let Some(bits) = cols.get_mut(c as usize) {
                bits.set(t);
            }
        }
    }
    cols
}

/// Does the row's residue sequence equal the protein over the window the row covers?
/// Case-insensitive; `X` matches anything, since alignment tools and annotation disagree on it.
fn row_matches_protein(row: &AlignmentRow, protein: &[u8]) -> bool {
    let start = row.residue_offset;
    if start + row.numbered > protein.len() || row.res2col.len() < row.numbered {
        return false;
    }
    (0..row.numbered).all(|q| {
        let a = row.aligned[row.res2col[q] as usize].to_ascii_uppercase();
        let b = protein[start + q].to_ascii_uppercase();
        a == b || a == b'X' || b == b'X'
    })
}

/// Which alignment taxa have annotation, and which annotated proteins lack an alignment row.
#[derive(Clone, Debug, Default)]
pub struct Coverage {
    /// Alignment rows whose residues differ from the protein they were joined to; they are scored as
    /// missing data because projecting through them would misplace every site.
    pub row_mismatch: Vec<String>,
    pub annotated_taxa: Vec<usize>,
    pub unannotated_taxa: Vec<String>,
    /// Annotated proteins with no row in this alignment, and the (lexicographically) first few of them.
    pub n_absent_from_alignment: usize,
    pub absent_examples: Vec<String>,
}

/// How many absent proteins `Coverage::absent_examples` keeps.
pub const ABSENT_EXAMPLES: usize = 5;

impl Coverage {
    pub fn n_annotated(&self) -> usize {
        self.annotated_taxa.len()
    }
}

/// Project one protein's sites onto one alignment row.
///
/// `residue → column` is the only mapping used, so a site lands exactly where its codon does in the
/// alignment — which is what makes two species' positions comparable.
pub fn project_row(
    taxon: usize,
    row: &AlignmentRow,
    sites: &[ProteinSite],
    opts: &ProjectOptions,
) -> (Vec<ColumnSite>, ProjectionReport) {
    let mut out = Vec::with_capacity(sites.len());
    let mut rep = ProjectionReport::default();
    // Clamp target is the last *covered* residue, which is not necessarily the last residue in the
    // row when Stockholm SF/SE excludes part of the source from the alignment.
    let last_col = row.last_column();
    for &site in sites {
        if opts.skip_incoherent && !site.is_coherent() {
            continue;
        }
        let residue = site.codon_index0() as usize;
        let mut was_clamped = false;
        let col = match row.column_of(residue) {
            Some(c) => {
                rep.placed += 1;
                c
            }
            None => {
                rep.tail_sites += 1;
                // Clamp only snaps a genuine C-terminal overhang onto the last covered column. A
                // residue before the Stockholm window, or a row with no residues, is not a tail.
                let past_end = residue >= row.residue_offset + row.numbered;
                match (opts.tail, last_col, past_end) {
                    (TailPolicy::Clamp, Some(c), true) => {
                        rep.clamped += 1;
                        was_clamped = true;
                        c
                    }
                    _ => continue,
                }
            }
        };
        out.push(ColumnSite {
            taxon: taxon as u32,
            col: col as u32,
            phase: match opts.phase_mode {
                PhaseMode::Strict => site.phase,
                PhaseMode::Ignore => 0,
            },
            site,
            clamped: was_clamped,
        });
    }
    (out, rep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::{AnnotationSource, GffSource};
    use crate::model::AnnotatedIntron;
    use std::io::Cursor;

    fn tx(protein: &str, offsets: &[u32]) -> TranscriptAnnotation {
        TranscriptAnnotation {
            transcript_id: format!("{protein}.t1"),
            gene_id: Some(format!("{protein}.g")),
            protein_id: Some(protein.to_string()),
            seqid: "chr1".into(),
            strand: Strand::Plus,
            // long enough to look real without implying a specific protein length
            cds_len: offsets.last().map(|o| *o as u64 + 3).unwrap_or(300),
            cds_segments: offsets.len() as u32 + 1,
            introns: offsets
                .iter()
                .enumerate()
                .map(|(i, &o)| AnnotatedIntron {
                    start: 100 + i as u64 * 50,
                    end: 120 + i as u64 * 50,
                    rank: i as u32,
                    cds_offset: Some(o),
                })
                .collect(),
            frame_mismatches: 0,
            source: "test".into(),
        }
    }

    fn proteins(fasta: &str) -> ProteinSet {
        let mut s = ProteinSet::new();
        s.read_one(Cursor::new(fasta.as_bytes()), std::path::Path::new("m"))
            .unwrap();
        s
    }

    fn aln(text: &str) -> Alignment {
        Alignment::from_fasta_str(text, "fam").unwrap()
    }

    /// Most fixtures here are about projection, not sequence length, so skip that check unless a
    /// test is specifically about it.
    fn no_len_check() -> JoinOptions {
        JoinOptions {
            length_slack: u64::MAX,
            ..Default::default()
        }
    }

    #[test]
    fn joins_and_collects_sites() {
        let ps = proteins(">A\nMKVLGAWE\n>B\nSHORTSEQ\n");
        let p = join_annotations(
            &[tx("A", &[15, 3]), tx("B", &[6]), tx("C", &[3])],
            &ps,
            &no_len_check(),
        )
        .unwrap();
        assert_eq!(p.joined_proteins, 2);
        assert_eq!(p.rows.len(), 3);
        assert_eq!(
            p.sites_of(0)
                .unwrap()
                .iter()
                .map(|s| s.cds_offset)
                .collect::<Vec<_>>(),
            vec![3, 15],
            "sites are stored sorted by offset, not in annotation order"
        );
        assert_eq!(p.unmatched.len(), 1);
        assert_eq!(p.unmatched[0].reason, "protein_not_in_fasta");
        assert_eq!(p.unmatched[0].attempted_id, "C");
        assert_eq!(p.n_distinct_sites(), 3);
    }

    #[test]
    fn flags_out_of_range_sites_and_length_mismatch() {
        let ps = proteins(">A\nMKVLGAWE\n"); // 8 residues
        let p = join_annotations(&[tx("A", &[120])], &ps, &JoinOptions::default()).unwrap();
        assert_eq!(p.diagnostics.count_of("site_beyond_protein_end"), 1);
        assert_eq!(p.diagnostics.count_of("length_mismatch"), 1);
        assert!(!p.rows[0].in_range);
    }

    #[test]
    fn length_slack_tolerates_the_stop_codon_convention() {
        let mut t = tx("A", &[15]);
        t.cds_len = 30; // 10 codons including a stop -> protein of 9 or 10 residues
        let ps = proteins(">A\nMKVLGAWEK\n");
        let p = join_annotations(&[t], &ps, &JoinOptions::default()).unwrap();
        assert_eq!(
            p.diagnostics.count_of("length_mismatch"),
            0,
            "{:?}",
            p.diagnostics.retained()
        );
    }

    #[test]
    fn version_insensitive_join_keeps_the_sequence_set_id() {
        let ps = proteins(">P1.2\nMKVLMKVL\n");
        let p = join_annotations(&[tx("P1", &[3])], &ps, &no_len_check()).unwrap();
        assert_eq!(p.joined_proteins, 1, "annotation carries no version");
        assert_eq!(
            p.protein_id(&p.rows[0]),
            "P1.2",
            "tables use the identifier that was actually matched"
        );
    }

    #[test]
    fn unannotated_and_unaligned_taxa_are_listed() {
        let ps = proteins(">A\nMKVLGAWE\n>B\nMKVLGAWE\n>C\nMKVLGAWE\n");
        let p = join_annotations(&[tx("A", &[3]), tx("C", &[3])], &ps, &no_len_check()).unwrap();
        let a = aln(">A\nMKVL--GAWE\n>B\nMKVL--GAWE\n");
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        assert_eq!(ap.coverage.annotated_taxa, vec![0]);
        assert_eq!(ap.coverage.unannotated_taxa, vec!["B"]);
        assert_eq!(ap.coverage.absent_examples, vec!["C"]);
        assert_eq!(ap.coverage.n_absent_from_alignment, 1);
        assert_eq!(ap.taxa[1].unannotated_reason, Some("protein_not_annotated"));
    }

    #[test]
    fn projection_maps_codon_to_column() {
        let ps = proteins(">A\nMKVLGAWE\n"); // 8 residues
        let p = join_annotations(&[tx("A", &[3, 12, 13])], &ps, &no_len_check()).unwrap();
        // The row gaps the first two residues, so A's residue i sits at column i + 2.
        let a = aln(">A\n--MKVLGAWE\n");
        let sites = project_into_alignment(&a, &p, &ps, &ProjectOptions::default()).sites;
        let cols: Vec<_> = sites.iter().map(|s| (s.col, s.phase)).collect();
        // cds_offset 3 -> codon idx 1 -> col 3; 12 and 13 -> codon idx 4 -> col 6, phases 0 and 1.
        assert_eq!(cols, vec![(3, 0), (6, 0), (6, 1)], "sorted by (col, phase)");
    }

    #[test]
    fn tail_policy_drops_or_clamps_sites_without_a_column() {
        let ps = proteins(">A\nMKVLGAWE\n");
        // cds_offset 21 -> 0-based codon 7, the 8th residue, which the protein has but the row lacks.
        let p = join_annotations(&[tx("A", &[21])], &ps, &no_len_check()).unwrap();
        let a = aln(">A\nMKVLGA--\n"); // only 6 residues aligned
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        assert!(
            ap.sites.is_empty(),
            "past-the-end sites are dropped by default"
        );
        assert_eq!(ap.taxa[0].report.tail_sites, 1);
        assert_eq!(ap.taxa[0].report.clamped, 0);

        let ac = project_into_alignment(
            &a,
            &p,
            &ps,
            &ProjectOptions {
                tail: TailPolicy::Clamp,
                ..Default::default()
            },
        );
        assert_eq!(ac.sites.len(), 1);
        assert_eq!(ac.sites[0].col, 5, "clamped to the row's last column");
        assert_eq!(ac.taxa[0].report.clamped, 1);
    }

    #[test]
    fn phase_mode_groups_by_column_only() {
        let ps = proteins(">A\nMKVLGAWE\n");
        let p = join_annotations(&[tx("A", &[12, 13, 14])], &ps, &no_len_check()).unwrap();
        let a = aln(">A\nMKVLGAWE\n");
        let strict = project_into_alignment(&a, &p, &ps, &ProjectOptions::default()).sites;
        assert_eq!(
            strict.iter().map(|s| s.phase).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        let ignore = project_into_alignment(
            &a,
            &p,
            &ps,
            &ProjectOptions {
                phase_mode: PhaseMode::Ignore,
                ..Default::default()
            },
        )
        .sites;
        assert!(
            ignore.iter().all(|s| s.phase == 0),
            "phases collapse onto one site per column"
        );
        assert_eq!(
            ignore
                .iter()
                .map(|s| s.col)
                .collect::<std::collections::HashSet<_>>()
                .len(),
            1
        );
    }

    #[test]
    fn homologous_sites_share_a_column_despite_renumbering() {
        // B carries a 5-residue deletion upstream, so its intron sits at residue 2 while A's sits at
        // residue 7 — different numbers, same alignment column.
        let ps = proteins(">A\nMKVLGAWEQP\n>B\nAWEQP\n");
        let mut b = tx("B", &[3]);
        b.transcript_id = "B.t1".into();
        let p = join_annotations(&[tx("A", &[18]), b], &ps, &no_len_check()).unwrap();
        let a = aln(">A\nMKVLGAWEQP\n>B\n-----AWEQP\n");
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        let by_taxon: Vec<_> = ap
            .sites
            .iter()
            .map(|s| (a.rows[s.taxon as usize].id.clone(), s.col, s.phase))
            .collect();
        assert_eq!(
            by_taxon,
            vec![("A".to_string(), 6, 0), ("B".to_string(), 6, 0)]
        );
        assert_eq!(
            ap.sites
                .iter()
                .map(|s| s.site.cds_offset)
                .collect::<Vec<_>>(),
            vec![18, 3],
            "the two species' own numbering differs by a whole indel"
        );
    }

    #[test]
    fn multiple_transcripts_union_their_sites() {
        let ps = proteins(">A\nMKVLGAWEQPRRRR\n");
        let mut t2 = tx("A", &[6]);
        t2.transcript_id = "A.t2".into();
        let p = join_annotations(&[tx("A", &[3]), t2], &ps, &no_len_check()).unwrap();
        assert_eq!(
            p.sites_of(0).unwrap().len(),
            2,
            "union of both isoforms' sites"
        );
        assert_eq!(p.n_transcripts_of[0], 2);
        assert!(p.diagnostics.count_of("isoform_intron_conflict") >= 1);
    }

    #[test]
    fn annotated_gene_with_no_introns_is_observed_absent_not_missing() {
        // The distinction that decides whether intron loss is measurable at all: species "lost" has
        // annotation and zero coding introns, so it must contribute absences. Species "missing" is in
        // the alignment but has no annotation record, so it must contribute nothing.
        let ps = proteins(">kept\nMKVLGAWEQP\n>lost\nMKVLGAWEQP\n>missing\nMKVLGAWEQP\n");
        let mut none = tx("lost", &[]);
        none.cds_len = 36;
        let p = join_annotations(&[tx("kept", &[3]), none], &ps, &no_len_check()).unwrap();
        assert!(
            p.is_annotated(ps.lookup("lost").unwrap()),
            "a transcript joined, so we know about it"
        );

        let a = aln(">kept\nMKVLGAWEQP\n>lost\nMKVLGAWEQP\n>missing\nMKVLGAWEQP\n");
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        let by_id = |id: &str| ap.taxa.iter().find(|t| t.row_id == id).unwrap();

        let lost = by_id("lost");
        assert!(
            lost.has_annotation,
            "annotation exists even with zero introns"
        );
        assert!(lost.zero_introns);
        assert_eq!(lost.unannotated_reason, None);
        assert!(lost.sites.is_empty());

        let gone = by_id("missing");
        assert!(!gone.has_annotation);
        assert_eq!(gone.unannotated_reason, Some("protein_not_annotated"));
        assert!(
            gone.in_sequence_set,
            "the sequence is known, only its annotation is absent"
        );

        assert_eq!(
            ap.coverage.annotated_taxa,
            vec![0, 1],
            "kept and lost both count"
        );
        assert_eq!(ap.coverage.unannotated_taxa, vec!["missing"]);
    }

    #[test]
    fn row_absent_from_the_sequence_set_is_flagged_separately() {
        let ps = proteins(">kept\nMKVLGAWEQP\n");
        let p = join_annotations(&[tx("kept", &[3])], &ps, &no_len_check()).unwrap();
        let a = aln(">kept\nMKVLGAWEQP\n>unsampled\nMKVLGAWEQP\n");
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        let gone = ap.taxa.iter().find(|t| t.row_id == "unsampled").unwrap();
        assert!(!gone.in_sequence_set, "no such protein in the FASTA at all");
        assert_eq!(gone.unannotated_reason, Some("protein_not_in_sequence_set"));
    }

    #[test]
    fn annotated_intronless_protein_absent_from_alignment_is_still_reported() {
        // "lost" has annotation and zero introns, and has no row here. It must appear in
        // absent_from_alignment: we have information about it that the alignment fails to use.
        let ps = proteins(">kept\nMKVLGAWEQP\n>lost\nMKVLGAWEQP\n");
        let mut none = tx("lost", &[]);
        none.cds_len = 36;
        let p = join_annotations(&[tx("kept", &[3]), none], &ps, &no_len_check()).unwrap();
        let a = aln(">kept\nMKVLGAWEQP\n");
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        assert_eq!(ap.coverage.absent_examples, vec!["lost"]);
        assert_eq!(ap.coverage.n_absent_from_alignment, 1);
    }

    #[test]
    fn observability_marks_columns_a_taxon_does_not_cover() {
        use crate::matrix::Bitset;
        let a = aln(">A\nMK--VL\n>B\n--MKVL\n");
        let obs = observability(&a);
        assert_eq!(obs.len(), a.len);
        // A holds residues in columns 0,1,4,5; B in 2,3,4,5.
        assert_eq!(obs[0].positions(), vec![0]);
        assert_eq!(obs[2].positions(), vec![1]);
        assert_eq!(obs[4].positions(), vec![0, 1], "both cover column 4");
        let empty = Bitset::new(2);
        assert_eq!(empty.positions(), Vec::<usize>::new());
    }

    #[test]
    fn identical_isoforms_are_not_flagged_as_conflicting() {
        let ps = proteins(">A\nMKVLGAWEQPRRRR\n");
        let mut t2 = tx("A", &[3]);
        t2.transcript_id = "A.t2".into();
        let p = join_annotations(&[tx("A", &[3]), t2], &ps, &no_len_check()).unwrap();
        assert_eq!(p.diagnostics.count_of("isoform_intron_conflict"), 0);
        assert_eq!(p.sites_of(0).unwrap().len(), 1, "the set deduplicates");
    }

    #[test]
    fn real_gff_round_trips_into_sites_and_columns() {
        let dir = std::env::temp_dir().join(format!("ip_join_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gff = dir.join("a.gff3");
        std::fs::write(
            &gff,
            "chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t1;protein_id=P1\n\
             chr1\tt\tCDS\t200\t221\t.\t+\t0\tID=c2;Parent=t1;protein_id=P1\n\
             chr1\tt\tCDS\t300\t320\t.\t+\t1\tID=c3;Parent=t1;protein_id=P1\n",
        )
        .unwrap();
        let set = GffSource::path(gff).load().unwrap();
        // 58 nt CDS -> 18..20 residues; 19 chosen so the 13:1 site (codon index 12) is in range.
        let ps = proteins(">P1\nMKVLGAWETSIGKDELAAQ\n");
        let p = join_annotations(&set.transcripts, &ps, &JoinOptions::default()).unwrap();
        assert_eq!(p.rows.len(), 2);
        assert_eq!((p.rows[0].site.codon(), p.rows[0].site.phase), (5, 0));
        assert_eq!((p.rows[1].site.codon(), p.rows[1].site.phase), (13, 1));
        assert_eq!(p.rows[1].genomic_start, 222);
        assert_eq!(p.rows[1].genomic_end, 299);
        assert_eq!(
            p.diagnostics.count_of("length_mismatch"),
            0,
            "{:?}",
            p.diagnostics.retained()
        );
        assert_eq!(p.diagnostics.count_of("site_beyond_protein_end"), 0);

        let a = aln(">P1\nMKVLGAWETSIGKDELAAQ\n");
        let ap = project_into_alignment(&a, &p, &ps, &ProjectOptions::default());
        assert_eq!(
            ap.sites
                .iter()
                .map(|s| (s.col, s.phase))
                .collect::<Vec<_>>(),
            vec![(5, 0), (12, 1)]
        );
    }

    // ---- proposed: hand-written truth table for cell_observable ----
    #[test]
    fn cell_observable_needs_exactly_the_documented_columns() {
        use crate::matrix::Bitset;
        // One taxon with a residue in every column of 0..12 except `hole`.
        let with_hole = |hole: usize| -> Vec<Bitset> {
            (0..12)
                .map(|c| {
                    let mut b = Bitset::new(1);
                    if c != hole {
                        b.set(0);
                    }
                    b
                })
                .collect()
        };
        let col = 6u32;
        // (phase, flank, columns that must hold a residue) -- written out by hand from the doc comment.
        let cases: &[(u8, u32, &[usize])] = &[
            (0, 0, &[5, 6]), // boundary: residues on both sides, flank 0 still needs 1
            (0, 1, &[5, 6]),
            (0, 2, &[4, 5, 6, 7]),
            (1, 0, &[6]), // inside a residue: only that residue
            (1, 1, &[5, 6, 7]),
            (2, 1, &[5, 6, 7]),
            (2, 2, &[4, 5, 6, 7, 8]),
        ];
        for &(phase, flank, needed) in cases {
            for hole in 0..12 {
                let want = !needed.contains(&hole);
                let got = cell_observable(&with_hole(hole), 0, col, Some(phase), flank);
                assert_eq!(
                    got, want,
                    "phase {phase} flank {flank}: gap at column {hole}"
                );
            }
        }
        // Column 0 boundary site has no left neighbour: never observable.
        assert!(!cell_observable(&with_hole(99), 0, 0, Some(0), 1));
        // A site merged from several phases (--ignore-phase) needs the union: c-max(n,1) ..= c+n.
        for (flank, needed) in [
            (0u32, vec![5usize, 6]),
            (1, vec![5, 6, 7]),
            (2, vec![4, 5, 6, 7, 8]),
        ] {
            for hole in 0..12 {
                let got = cell_observable(&with_hole(hole), 0, col, None, flank);
                assert_eq!(
                    got,
                    !needed.contains(&hole),
                    "merged, flank {flank}: gap at {hole}"
                );
            }
        }
    }
}
