//! Turning annotation sources into per-transcript splicing structure.
//!
//! [`AnnotationSource`] is the seam between "where features come from" and "what we compute".
//! [`GffSource`] implements it today; a Chado/Ensembl/SQLite reader can implement the same trait and
//! every downstream stage — projection, matrices, scoring — works unchanged.
//!
//! The trait is visitor-shaped (`stream_into(&mut self, sink)`) rather than iterator-shaped because
//! both real sources are naturally push-based: GFF finalises transcripts at EOF, and a SQL cursor
//! hands over rows as they stream.
//!
//! ## How a protein-space intron position is derived
//!
//! CDS features of a transcript are ordered along the mRNA (genomic start descending for minus
//! strand), contiguous and overlapping pieces are merged, and each remaining gap is an intron whose
//! `cds_offset` is the number of coding bases preceding it. That count is invariant to UTR length,
//! intron length and strand, and it is checked two ways:
//!
//! 1. against the exon structure, when `exon` features exist — every coding gap must also be an exon
//!    gap, and exon gaps that are not coding gaps are UTR introns (kept, but not projectable);
//! 2. against each CDS feature's `phase`/`frame` attribute, which by spec is the number of bases to
//!    skip from the feature's start to reach the next codon boundary, i.e.
//!    `(3 - coding_bases_before % 3) % 3` (see [`expected_frame`]). Disagreement means the annotation is internally
//!    inconsistent, so the transcript is flagged rather than silently trusted.

use std::io::BufRead;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::error::{Error, Result};
use crate::gff::GffReader;
use crate::model::{AnnotatedIntron, Strand, TranscriptAnnotation};
use crate::report::{Diagnostic, Diagnostics, Level};
use crate::util;

/// Receives transcripts and diagnostics as a source produces them.
pub trait AnnotationSink {
    fn push_transcript(&mut self, tx: TranscriptAnnotation) -> Result<()>;
    fn diagnostic(&mut self, level: Level, kind: &'static str, subject: String, message: String);
    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

/// Anything that can supply transcript-level splicing structure.
pub trait AnnotationSource {
    /// One-line description for logs and report headers.
    fn describe(&self) -> String;
    fn stream_into(&mut self, sink: &mut dyn AnnotationSink) -> Result<()>;

    /// Convenience for tests and small inputs.
    fn load(&mut self) -> Result<AnnotationSet> {
        let mut set = AnnotationSet::default();
        self.stream_into(&mut set)?;
        Ok(set)
    }
}

/// Accumulated transcripts plus the diagnostics raised while reading them.
#[derive(Clone, Debug, Default)]
pub struct AnnotationSet {
    pub transcripts: Vec<TranscriptAnnotation>,
    pub diagnostics: Diagnostics,
}

impl AnnotationSink for AnnotationSet {
    fn push_transcript(&mut self, tx: TranscriptAnnotation) -> Result<()> {
        self.transcripts.push(tx);
        Ok(())
    }

    fn diagnostic(&mut self, level: Level, kind: &'static str, subject: String, message: String) {
        self.diagnostics.push(level, kind, subject, message);
    }
}

/// An already-materialised source — useful for tests and as a cache in front of a slow source.
#[derive(Clone, Debug, Default)]
pub struct MemorySource {
    pub transcripts: Vec<TranscriptAnnotation>,
    pub label: String,
}

impl MemorySource {
    pub fn new(transcripts: Vec<TranscriptAnnotation>) -> Self {
        MemorySource {
            transcripts,
            label: "mem".into(),
        }
    }
}

impl AnnotationSource for MemorySource {
    fn describe(&self) -> String {
        format!(
            "memory:{}",
            if self.label.is_empty() {
                "transcripts"
            } else {
                self.label.as_str()
            }
        )
    }

    fn stream_into(&mut self, sink: &mut dyn AnnotationSink) -> Result<()> {
        for tx in self.transcripts.iter().cloned() {
            sink.push_transcript(tx)?;
        }
        sink.finish()
    }
}

/// Which feature types and attribute names to recognise.
///
/// Every list is matched case-insensitively and tried in order, because annotation pipelines
/// disagree hopelessly on conventions. The defaults cover GFF3-canonical, Ensembl GTF, NCBI RefSeq,
/// MAKER, Augustus/Braker and protein2genome output.
#[derive(Clone, Debug)]
pub struct GffOptions {
    pub cds_types: Vec<String>,
    pub exon_types: Vec<String>,
    pub intron_types: Vec<String>,
    pub transcript_types: Vec<String>,
    pub gene_types: Vec<String>,
    pub protein_types: Vec<String>,
    /// Attribute order used to name the transcript a CDS/exon/intron belongs to.
    pub transcript_attrs: Vec<String>,
    /// Attribute order used to name the protein a transcript encodes.
    pub protein_attrs: Vec<String>,
    /// Attribute order used to name the parent gene. Deliberately excludes `ID`/`Name`, which on a
    /// CDS row name the CDS itself rather than its gene.
    pub gene_attrs: Vec<String>,
    /// Explicit `annotation_id -> protein_id` rewrites, applied after other resolution.
    pub id_map: FxHashMap<String, String>,
    /// Drop trailing numeric versions (`.1`, `-2`) when matching ids.
    pub strip_version: bool,
    /// Remove a leading `rna:` / `transcript:` / … namespace from ids.
    pub strip_id_namespace: bool,
    /// Keep transcripts with no coding structure (their introns can never be projected).
    pub keep_noncoding: bool,
    /// When a transcript has no CDS features, treat its exons as coding, assuming translation starts
    /// at the first exon base. Needed for protein2genome/Exonerate-style GFFs.
    pub derive_cds_from_exons: bool,
    /// Fail on malformed lines instead of counting and skipping them.
    pub strict: bool,
    /// Cap on retained diagnostic rows.
    pub diagnostic_cap: usize,
}

impl Default for GffOptions {
    fn default() -> Self {
        let v = |items: &[&str]| items.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        GffOptions {
            cds_types: v(&["cds"]),
            exon_types: v(&["exon"]),
            intron_types: v(&["intron"]),
            transcript_types: v(&["mrna", "transcript", "pseudogenic_transcript", "cdna"]),
            gene_types: v(&["gene", "pseudogene"]),
            protein_types: v(&["protein", "polypeptide", "translation"]),
            transcript_attrs: v(&["Parent", "transcript_id", "transcript", "rna_id"]),
            protein_attrs: v(&["protein_id", "protein"]),
            gene_attrs: v(&["gene_id", "gene", "locus_tag"]),
            id_map: FxHashMap::default(),
            strip_version: false,
            strip_id_namespace: true,
            keep_noncoding: false,
            derive_cds_from_exons: false,
            strict: false,
            diagnostic_cap: 500,
        }
    }
}

/// Ids that lead a GFF `ID`/`Parent` value as a namespace, e.g. Ensembl's `rna:ENST000001`.
const ID_NAMESPACES: &[&str] = &[
    "rna",
    "transcript",
    "mrna",
    "gene",
    "cds",
    "exon",
    "protein",
    "polypeptide",
    "ncrna",
];

fn normalize_id(id: &str, opts: &GffOptions) -> String {
    let mut id = id.trim();
    if opts.strip_id_namespace {
        if let Some((prefix, rest)) = id.split_once(':') {
            if !rest.is_empty() && ID_NAMESPACES.iter().any(|p| p.eq_ignore_ascii_case(prefix)) {
                id = rest;
            }
        }
    }
    // Some tools append a bracketed list or a parenthetical to Parent values.
    id = id.split(['[', '(']).next().unwrap_or(id).trim();
    if opts.strip_version {
        id = util::strip_version(id);
    }
    id.to_string()
}

/// Record a diagnostic. Counts are never truncated here; the sink's `Diagnostics` keeps full counts and
/// caps only the retained rows, so one pathological input cannot flood the report or hide its own size.
fn note(
    diags: &mut Vec<Diagnostic>,
    level: Level,
    kind: &'static str,
    subject: &str,
    message: impl Into<String>,
) {
    diags.push(Diagnostic {
        level,
        kind,
        subject: subject.to_string(),
        message: message.into(),
    });
}

/// Reads one or more GFF3/GTF files.
#[derive(Clone, Debug)]
pub struct GffSource {
    pub paths: Vec<PathBuf>,
    pub opts: GffOptions,
}

impl GffSource {
    pub fn new(paths: Vec<PathBuf>, opts: GffOptions) -> Self {
        GffSource { paths, opts }
    }

    pub fn path(path: impl Into<PathBuf>) -> Self {
        GffSource {
            paths: vec![path.into()],
            opts: GffOptions::default(),
        }
    }

    /// Parse an `annotation_id,protein_id` mapping file (TSV or CSV, optional `#` header).
    pub fn load_id_map(path: impl AsRef<Path>) -> Result<FxHashMap<String, String>> {
        let path = path.as_ref();
        let mut map = FxHashMap::default();
        let mut reader = util::open_buffered(path, 1 << 16)?;
        let mut raw = String::new();
        let mut i = 0usize;
        while reader.read_line(&mut raw)? > 0 {
            i += 1;
            let line = raw.trim().to_string();
            raw.clear();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let sep = if line.contains('\t') { '\t' } else { ',' };
            let mut parts = line.split(sep).map(str::trim);
            let (Some(from), Some(to)) = (parts.next(), parts.next()) else {
                return Err(Error::Parse {
                    message: format!("expected two columns, got {line:?}"),
                    location: Some(format!("{}:{}", path.display(), i + 1)),
                });
            };
            map.insert(from.to_string(), to.to_string());
        }
        Ok(map)
    }

    fn classify(&self, ty: &str) -> Kind {
        let has = |list: &[String]| list.iter().any(|t| t.eq_ignore_ascii_case(ty));
        if has(&self.opts.cds_types) {
            Kind::Cds
        } else if has(&self.opts.exon_types) {
            Kind::Exon
        } else if has(&self.opts.intron_types) {
            Kind::Intron
        } else if has(&self.opts.transcript_types) {
            Kind::Transcript
        } else if has(&self.opts.protein_types) {
            Kind::Protein
        } else if has(&self.opts.gene_types) {
            Kind::Gene
        } else {
            Kind::Other
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Cds,
    Exon,
    Intron,
    Transcript,
    Gene,
    Protein,
    Other,
}

#[derive(Clone, Copy, Debug)]
struct Seg {
    start: u64,
    end: u64,
    frame: Option<u8>,
}

/// A CDS segment plus the frame it was annotated with, retained after merging so the reading-frame
/// check sees the original value.
#[derive(Clone, Copy, Debug)]
struct AnnSeg {
    seg: Seg,
    orig_frame: Option<u8>,
}

#[derive(Clone, Debug, Default)]
struct TxAcc {
    transcript_id: String,
    gene: Option<String>,
    protein: Option<String>,
    seqid: Option<String>,
    strand: Option<Strand>,
    cds: Vec<Seg>,
    exons: Vec<(u64, u64)>,
    introns: Vec<(u64, u64)>,
    /// `ID` values of this transcript's CDS features, used to resolve `protein` child features.
    cds_ids: Vec<String>,
    n_cds_records: u32,
    frame_mismatches: u32,
    source: String,
}

#[derive(Default, Debug)]
struct Stats {
    lines: usize,
    features: u64,
    records: u64,
    malformed: u64,
    directives: u64,
}

impl AnnotationSource for GffSource {
    fn describe(&self) -> String {
        format!(
            "gff3:[{}]",
            self.paths
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    }

    fn stream_into(&mut self, sink: &mut dyn AnnotationSink) -> Result<()> {
        // Files are independent (one species each): read and finalise them in parallel, so the raw
        // per-transcript features of only one file per worker thread are ever held, and results are merged
        // in file order, which keeps the output deterministic.
        let this = &*self;
        let multi = self.paths.len() > 1;
        // Only a 64-bit hash of each id is kept for the cross-file collision check (an exact set of ten
        // million strings would cost ~0.6 GB); a false positive is a spurious warning, never a lost record.
        let mut seen: FxHashSet<u64> = FxHashSet::default();
        let mut stats = Stats::default();
        let mut kept = 0usize;
        // Files are processed in batches the width of the thread pool and each batch is merged into the sink
        // before the next is read, so a sink that consumes transcripts as they arrive (see
        // `JoiningSink`) never has the transcripts of more than one batch in memory.
        let width = rayon::current_num_threads().max(1);
        for batch in self.paths.chunks(width) {
            let results: Vec<Result<FileResult>> =
                batch.par_iter().map(|path| this.read_file(path)).collect();
            for result in results {
                let fr = result?;
                stats.features += fr.stats.features;
                stats.records += fr.stats.records;
                stats.lines += fr.stats.lines;
                stats.malformed += fr.stats.malformed;
                stats.directives += fr.stats.directives;
                for (level, kind, subject, message) in fr.file_diagnostics {
                    sink.diagnostic(level, kind, subject, message);
                }
                for d in fr.diagnostics {
                    sink.diagnostic(d.level, d.kind, d.subject, d.message);
                }
                for (tx, diags) in fr.built {
                    // Diagnostics are forwarded whether or not the transcript survived.
                    for d in diags {
                        sink.diagnostic(d.level, d.kind, d.subject, d.message);
                    }
                    let Some(tx) = tx else { continue };
                    kept += 1;
                    if multi {
                        use std::hash::{Hash, Hasher};
                        let mut h = rustc_hash::FxHasher::default();
                        tx.transcript_id.hash(&mut h);
                        if !seen.insert(h.finish()) {
                            sink.diagnostic(
                            Level::Warn,
                            "transcript_id_collision",
                            tx.transcript_id.clone(),
                            format!(
                                "id also appears in an earlier file; both records are kept (this one from {:?})",
                                tx.source
                            ),
                        );
                        }
                    }
                    sink.push_transcript(tx)?;
                }
            }
        }
        sink.diagnostic(
            Level::Info,
            "source_summary",
            self.describe(),
            format!(
                "{kept} transcripts retained from {} feature rows over {} lines ({} malformed, {} directives)",
                stats.features, stats.lines, stats.malformed, stats.directives
            ),
        );
        sink.finish()
    }
}

/// What reading one annotation file produced.
struct FileResult {
    /// Notes about the file as a whole (invalid UTF-8, no records).
    file_diagnostics: Vec<(Level, &'static str, String, String)>,
    /// Notes raised while grouping features (missing keys, ...).
    diagnostics: Vec<Diagnostic>,
    built: Vec<(Option<TranscriptAnnotation>, Vec<Diagnostic>)>,
    stats: Stats,
}

impl GffSource {
    /// Read, group and finalise one annotation file, independently of every other file.
    fn read_file(&self, path: &PathBuf) -> Result<FileResult> {
        let mut accs: FxHashMap<String, TxAcc> = FxHashMap::default();
        let mut protein_children: FxHashMap<String, String> = FxHashMap::default();
        let mut diagnostics: Vec<Diagnostic> = Vec::new();
        let mut stats = Stats::default();
        let mut file_diagnostics = Vec::new();

        let label = path.display().to_string();
        let mut reader = GffReader::open(path)?;
        reader.strict = self.opts.strict;
        while let Some(rec) = reader.next_record()? {
            let kind = self.classify(rec.feature_type);
            if kind == Kind::Other || kind == Kind::Gene {
                continue;
            }
            stats.features += 1;

            // Every kept branch below stores copies; nothing borrowed from `rec` survives the
            // iteration, so the next `next_record` call cannot alias it.
            match kind {
                Kind::Protein => {
                    let Some(pid) = rec.attr("ID").or_else(|| rec.attr("protein_id")) else {
                        continue;
                    };
                    let pid = normalize_id(pid, &self.opts);
                    let mut linked = false;
                    for p in rec.parents() {
                        protein_children.insert(normalize_id(p, &self.opts), pid.clone());
                        linked = true;
                    }
                    if !linked {
                        protein_children.insert(pid.clone(), pid);
                    }
                    continue;
                }
                Kind::Transcript => {
                    let Some(id) = rec.attr("ID") else { continue };
                    let key = normalize_id(id, &self.opts);
                    let acc = entry(&mut accs, &label, key);
                    if acc.gene.is_none() {
                        acc.gene = rec.parents().next().map(str::to_string);
                    }
                    if let Some(p) = rec.attr_first_of(&self.opts.protein_attrs) {
                        acc.protein
                            .get_or_insert_with(|| normalize_id(p, &self.opts));
                    }
                    acc.seqid.get_or_insert_with(|| rec.seqid.to_string());
                    if rec.strand.is_known() {
                        acc.strand.get_or_insert(rec.strand);
                    }
                    continue;
                }
                _ => {}
            }

            // Which attribute names the transcript? `Parent` is tried first (GFF3), so NCBI-style files,
            // whose exon rows carry `transcript_id=NM_x` while CDS rows only carry `Parent=rna-NM_x`,
            // key both under the same transcript. A multi-valued `Parent` feeds every parent.
            let found = self
                .opts
                .transcript_attrs
                .iter()
                .find_map(|k| rec.attr(k).map(|v| (k.as_str(), v)));
            let Some((key_attr, raw_key)) = found else {
                note(
                    &mut diagnostics,
                    Level::Warn,
                    "no_transcript_key",
                    &format!("{label}:{}", rec.line_no),
                    format!(
                        "{} has no {} attribute",
                        rec.summary(),
                        self.opts.transcript_attrs.join("|")
                    ),
                );
                continue;
            };
            let keys: Vec<String> = if key_attr.eq_ignore_ascii_case("Parent") {
                rec.parents().map(|p| normalize_id(p, &self.opts)).collect()
            } else {
                vec![normalize_id(raw_key, &self.opts)]
            };
            let strand = rec.strand;
            let (start, end) = (rec.start, rec.end);
            let frame = rec.frame;
            let seqid = rec.seqid.to_string();
            let cds_id = rec.attr("ID").map(|id| normalize_id(id, &self.opts));
            let gene_attr = rec.attr_first_of(&self.opts.gene_attrs).map(str::to_string);
            let protein_attr = rec
                .attr_first_of(&self.opts.protein_attrs)
                .map(|p| normalize_id(p, &self.opts));

            for key in keys {
                let acc = entry(&mut accs, &label, key);
                acc.seqid.get_or_insert_with(|| seqid.clone());
                if strand.is_known() {
                    acc.strand.get_or_insert(strand);
                }
                if acc.gene.is_none() {
                    acc.gene = gene_attr.clone();
                }
                if acc.protein.is_none() {
                    acc.protein = protein_attr.clone();
                }
                match kind {
                    Kind::Cds => {
                        acc.n_cds_records += 1;
                        acc.cds.push(Seg { start, end, frame });
                        if let Some(id) = cds_id.clone() {
                            if acc.cds_ids.len() < 8 {
                                acc.cds_ids.push(id);
                            }
                        }
                    }
                    Kind::Exon => acc.exons.push((start, end)),
                    Kind::Intron => acc.introns.push((start, end)),
                    _ => {}
                }
            }
        }
        stats.records += reader.records;
        stats.lines = reader.line_no();
        stats.malformed += reader.malformed;
        stats.directives += reader.directives;
        if reader.lossy > 0 {
            file_diagnostics.push((
                Level::Warn,
                "invalid_utf8",
                label.clone(),
                format!(
                    "{} line(s) held invalid UTF-8 and were read lossily",
                    reader.lossy
                ),
            ));
        }
        if reader.records == 0 {
            file_diagnostics.push((
                Level::Warn,
                "empty_file",
                label.clone(),
                "no feature lines were parsed".into(),
            ));
        }

        // Deterministic order regardless of hash order, and parallel finalisation because the
        // per-transcript work (sort, merge, gap derivation, cross-checks) dominates here.
        let mut entries: Vec<TxAcc> = accs.into_values().collect();
        entries.sort_unstable_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
        let opts = &self.opts;
        let children = &protein_children;
        let finalised: Vec<(Option<TranscriptAnnotation>, Vec<Diagnostic>, bool)> = entries
            .into_par_iter()
            .map(|acc| finalize_tx(acc, opts, children))
            .collect();
        let n_fallback = finalised.iter().filter(|f| f.0.is_some() && f.2).count();
        if n_fallback > 0 {
            // One record per file: this is the normal case for funannotate, whose GFF has no protein ids.
            file_diagnostics.push((
                Level::Info,
                "protein_id_fallback_transcript",
                label.clone(),
                format!("{n_fallback} transcript(s) carried no protein id and were joined on their transcript id"),
            ));
        }
        let built: Vec<(Option<TranscriptAnnotation>, Vec<Diagnostic>)> =
            finalised.into_iter().map(|(t, d, _)| (t, d)).collect();
        Ok(FileResult {
            file_diagnostics,
            diagnostics,
            built,
            stats,
        })
    }
}

/// Fetch or create a transcript accumulator. Files are read independently, so a key belongs to one file;
/// ids repeated *across* files are reported when the per-file results are merged.
fn entry<'a>(accs: &'a mut FxHashMap<String, TxAcc>, label: &str, key: String) -> &'a mut TxAcc {
    accs.entry(key.clone()).or_insert_with(|| TxAcc {
        transcript_id: key,
        source: label.to_string(),
        ..Default::default()
    })
}

/// Coding structure of one transcript → a [`TranscriptAnnotation`], or `None` when there is nothing
/// projectable and `keep_noncoding` is off.
fn finalize_tx(
    mut acc: TxAcc,
    opts: &GffOptions,
    protein_children: &FxHashMap<String, String>,
) -> (Option<TranscriptAnnotation>, Vec<Diagnostic>, bool) {
    let mut diags: Vec<Diagnostic> = Vec::new();
    // Counted per file by the caller (one record), not noted per transcript: a genome set has ~10⁷.
    let mut used_fallback = false;
    let id = acc.transcript_id.as_str();

    if acc.cds.is_empty() && !acc.exons.is_empty() && opts.derive_cds_from_exons {
        acc.cds = acc
            .exons
            .iter()
            .map(|&(s, e)| Seg {
                start: s,
                end: e,
                frame: None,
            })
            .collect();
        acc.n_cds_records = acc.cds.len() as u32;
        note(
            &mut diags,
            Level::Warn,
            "cds_from_exons",
            id,
            "no CDS features: exons were treated as coding",
        );
    }

    let strand = match acc.strand {
        Some(s) if s.is_known() => s,
        _ => {
            let inferred = infer_strand(&acc);
            note(
                &mut diags,
                Level::Warn,
                "strand_unknown",
                id,
                format!("strand not annotated, inferred {inferred}"),
            );
            inferred
        }
    };

    // --- coding structure -------------------------------------------------------------
    let mut coding_gaps: Vec<(u64, u64, u32)> = Vec::new(); // (lo, hi, cds_offset)
    let mut cds_len = 0u64;
    if !acc.cds.is_empty() {
        let mut segs = order_segments(&acc.cds, strand);
        merge_contiguous(&mut segs, strand, &mut diags, id);
        cds_len = segs.iter().map(|s| s.seg.len()).sum();
        // A 5'-partial model starts mid-codon: its first feature's `phase` says how many leading bases
        // to skip before the first complete codon. Protein space begins there, so every offset is
        // counted from the first full codon, not from the first CDS base.
        let skip = leading_skip(&segs);
        cds_len = cds_len.saturating_sub(skip);
        let mut cum = 0u64;
        for i in 0..segs.len() {
            if let Some(f) = segs[i].orig_frame {
                if f != expected_frame(cum as i64 - skip as i64) {
                    acc.frame_mismatches += 1;
                }
            }
            cum += segs[i].seg.len();
            if i + 1 < segs.len() {
                let (lo, hi) = gap_between(strand, segs[i].seg, segs[i + 1].seg);
                if hi >= lo {
                    match cum.checked_sub(skip) {
                        Some(off) => coding_gaps.push((lo, hi, off as u32)),
                        None => note(
                            &mut diags,
                            Level::Warn,
                            "intron_before_first_codon",
                            id,
                            "intron lies inside the skipped 5' partial codon, not projectable",
                        ),
                    }
                }
            }
        }
    }

    // --- every intron of the transcript, in transcript order --------------------------
    let mut ordered_gaps: Vec<(u64, u64)> = if !acc.exons.is_empty() {
        exon_gaps(&acc.exons, strand, &mut diags, id)
    } else if !acc.introns.is_empty() {
        let mut g = acc.introns.clone();
        g.sort_unstable();
        if strand == Strand::Minus {
            g.reverse();
        }
        g
    } else {
        coding_gaps.iter().map(|&(lo, hi, _)| (lo, hi)).collect()
    };

    // A coding gap that is not an exon gap means CDS and exon features disagree; keep it so the
    // reading frame stays complete, and say so.
    if !coding_gaps.is_empty() {
        let mut missing = 0;
        for &(lo, hi, _) in &coding_gaps {
            if !ordered_gaps.contains(&(lo, hi)) {
                missing += 1;
                ordered_gaps.push((lo, hi));
            }
        }
        if missing > 0 {
            note(
                &mut diags,
                Level::Warn,
                "cds_gap_not_exon_gap",
                id,
                format!(
                    "{missing} coding gap(s) are not exon gaps; CDS and exon features disagree"
                ),
            );
            ordered_gaps.sort_unstable();
            ordered_gaps.dedup();
            if strand == Strand::Minus {
                ordered_gaps.reverse();
            }
        }
    }

    if !acc.introns.is_empty() && !coding_gaps.is_empty() {
        let annotated: FxHashSet<(u64, u64)> = acc.introns.iter().copied().collect();
        let unmatched = coding_gaps
            .iter()
            .filter(|g| !annotated.contains(&(g.0, g.1)))
            .count();
        if unmatched > 0 {
            note(
                &mut diags,
                Level::Warn,
                "intron_feature_mismatch",
                id,
                format!(
                    "{unmatched} coding gap(s) have no matching `intron` feature ({} explicit)",
                    acc.introns.len()
                ),
            );
        }
    }

    if acc.frame_mismatches > 0 {
        note(
            &mut diags,
            Level::Warn,
            "frame_mismatch",
            id,
            format!(
                "{} CDS features whose `phase` disagrees with the derived reading frame",
                acc.frame_mismatches
            ),
        );
    }

    let have_coding = !acc.cds.is_empty();
    if !have_coding && !opts.keep_noncoding {
        note(
            &mut diags,
            Level::Warn,
            "no_cds",
            id,
            "no CDS features: introns cannot be projected, dropped",
        );
        return (None, diags, false);
    }

    let gap_offsets: FxHashMap<(u64, u64), u32> = coding_gaps
        .iter()
        .map(|&(lo, hi, o)| ((lo, hi), o))
        .collect();
    let introns: Vec<AnnotatedIntron> = ordered_gaps
        .iter()
        .enumerate()
        .map(|(rank, &(lo, hi))| AnnotatedIntron {
            start: lo,
            end: hi,
            rank: rank as u32,
            cds_offset: gap_offsets.get(&(lo, hi)).copied(),
        })
        .collect();

    // --- protein identity -------------------------------------------------------------
    let protein_id = acc.protein.clone().or_else(|| {
        acc.cds_ids
            .iter()
            .chain(std::iter::once(&acc.transcript_id))
            .find_map(|k| protein_children.get(k).cloned())
    });
    let protein_id = protein_id.or_else(|| {
        used_fallback = true;
        Some(acc.transcript_id.clone())
    });
    let protein_id = protein_id.map(|p| match opts.id_map.get(&p) {
        Some(mapped) => mapped.clone(),
        None if opts.strip_version => util::strip_version(&p).to_string(),
        None => p,
    });

    (
        Some(TranscriptAnnotation {
            transcript_id: acc.transcript_id,
            gene_id: acc.gene,
            protein_id,
            seqid: acc.seqid.unwrap_or_default(),
            strand,
            cds_len,
            cds_segments: acc.n_cds_records,
            introns,
            frame_mismatches: acc.frame_mismatches,
            source: acc.source,
        }),
        diags,
        used_fallback,
    )
}

impl Seg {
    const fn len(self) -> u64 {
        self.end - self.start + 1
    }
}

fn order_segments(segs: &[Seg], strand: Strand) -> Vec<AnnSeg> {
    let mut v: Vec<AnnSeg> = segs
        .iter()
        .map(|&s| AnnSeg {
            seg: s,
            orig_frame: s.frame,
        })
        .collect();
    v.sort_unstable_by_key(|s| s.seg.start);
    if strand == Strand::Minus {
        v.reverse();
    }
    v
}

/// Merge overlapping and directly contiguous CDS pieces, warning on genuine overlaps.
///
/// Comparison happens in **transcript order**: on the minus strand the next piece in the CDS has
/// *lower* genomic coordinates, so testing `next.start <= cur.end + 1` in genomic space would treat
/// every minus-strand intron as a contiguous overlap and silently erase it. Coordinates are mapped
/// through a monotone-increasing transform, merged, and mapped back.
fn merge_contiguous(
    segs: &mut Vec<AnnSeg>,
    strand: Strand,
    diags: &mut Vec<Diagnostic>,
    subject: &str,
) {
    // Transcript-space bounds: ascending along the mRNA for either strand.
    let to_tx = |s: Seg| -> (i64, i64) {
        match strand {
            Strand::Minus => (-(s.end as i64), -(s.start as i64)),
            _ => (s.start as i64, s.end as i64),
        }
    };
    let from_tx = |a: i64, b: i64, frame: Option<u8>| -> Seg {
        match strand {
            Strand::Minus => Seg {
                start: (-b) as u64,
                end: (-a) as u64,
                frame,
            },
            _ => Seg {
                start: a as u64,
                end: b as u64,
                frame,
            },
        }
    };

    let mut out: Vec<(i64, i64, Option<u8>)> = Vec::with_capacity(segs.len());
    for s in segs.iter() {
        let (a, b) = to_tx(s.seg);
        match out.last_mut() {
            Some(last) if a <= last.1.saturating_add(1) => {
                if a <= last.1 {
                    note(
                        diags,
                        Level::Warn,
                        "overlapping_cds",
                        subject,
                        format!(
                            "CDS {}-{} overlaps {}-{}, merged",
                            s.seg.start,
                            s.seg.end,
                            from_tx(last.0, last.1, None).start,
                            from_tx(last.0, last.1, None).end
                        ),
                    );
                }
                last.1 = last.1.max(b);
                last.2 = last.2.or(s.orig_frame);
            }
            _ => out.push((a, b, s.orig_frame)),
        }
    }
    *segs = out
        .into_iter()
        .map(|(a, b, f)| AnnSeg {
            seg: from_tx(a, b, f),
            orig_frame: f,
        })
        .collect();
}

/// Genomic (low, high) bounds of the intron between two consecutive coding segments given in
/// transcript order.
fn gap_between(strand: Strand, cur: Seg, next: Seg) -> (u64, u64) {
    match strand {
        Strand::Minus => (next.end + 1, cur.start.saturating_sub(1)),
        // Plus and Unknown both assume ascending coordinates.
        _ => (cur.end + 1, next.start.saturating_sub(1)),
    }
}

/// Introns implied by the exon structure, in transcript order.
fn exon_gaps(
    exons: &[(u64, u64)],
    strand: Strand,
    diags: &mut Vec<Diagnostic>,
    subject: &str,
) -> Vec<(u64, u64)> {
    let mut e = exons.to_vec();
    e.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(e.len());
    for (s, t) in e {
        match merged.last_mut() {
            Some(last) if s <= last.1.saturating_add(1) => last.1 = last.1.max(t),
            _ => merged.push((s, t)),
        }
    }
    if merged.len() < exons.len() {
        note(
            diags,
            Level::Warn,
            "overlapping_exons",
            subject,
            format!(
                "{} overlapping/contiguous exon feature(s) were merged",
                exons.len() - merged.len()
            ),
        );
    }
    if strand == Strand::Minus {
        merged.reverse();
    }
    let mut gaps = Vec::with_capacity(merged.len().saturating_sub(1));
    for w in merged.windows(2) {
        let (lo, hi) = if strand == Strand::Minus {
            (w[1].1 + 1, w[0].0.saturating_sub(1))
        } else {
            (w[0].1 + 1, w[1].0.saturating_sub(1))
        };
        if hi >= lo {
            gaps.push((lo, hi));
        }
    }
    gaps
}

/// The GFF3 `phase` (GTF `frame`) a CDS feature must carry when `coding_before` coding bases precede
/// it in transcript order.
///
/// The attribute is *not* the number of bases already emitted: it is how many bases to remove from
/// the start of the feature to reach the first base of the next complete codon. After 1 leftover
/// base the answer is 2, after 2 leftover bases it is 1.
pub(crate) fn expected_frame(coding_before: i64) -> u8 {
    ((3 - coding_before.rem_euclid(3)) % 3) as u8
}

/// Bases of the first CDS feature, in transcript order, that precede the first complete codon (its
/// annotated `phase`; 0 when absent or out of range).
fn leading_skip(segs: &[AnnSeg]) -> u64 {
    segs.first()
        .and_then(|s| s.orig_frame)
        .filter(|f| *f <= 2)
        .map_or(0, u64::from)
}

/// Best guess at strand when the annotation omits it: whichever direction makes more CDS frames
/// consistent, falling back to ascending coordinates.
fn infer_strand(acc: &TxAcc) -> Strand {
    if acc.cds.len() < 2 {
        return Strand::Plus;
    }
    let frames_ok = |strand: Strand| -> u64 {
        let segs = order_segments(&acc.cds, strand);
        let skip = leading_skip(&segs) as i64;
        segs.iter()
            .scan(0i64, |cum, s| {
                let ok = s
                    .orig_frame
                    .is_some_and(|f| f == expected_frame(*cum - skip))
                    as u64;
                *cum += s.seg.len() as i64;
                Some(ok)
            })
            .sum()
    };
    if frames_ok(Strand::Minus) > frames_ok(Strand::Plus) {
        Strand::Minus
    } else {
        Strand::Plus
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_with(body: &str, tweak: impl FnOnce(&mut GffOptions)) -> AnnotationSet {
        // A per-call directory: fixtures otherwise overwrite each other's GFF when two bodies happen
        // to be the same length.
        let dir = std::env::temp_dir().join(format!(
            "ip_test_{}_{}",
            std::process::id(),
            ANCHOR.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("in.gff3");
        std::fs::write(&path, body).unwrap();
        let mut src = GffSource::path(path);
        tweak(&mut src.opts);
        src.load().unwrap()
    }

    static ANCHOR: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

    fn load(body: &str) -> AnnotationSet {
        load_with(body, |_| {})
    }

    #[test]
    fn derives_coding_offsets_on_plus_strand() {
        let s = load(
            "##gff-version 3\n\
chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t1;protein_id=P1\n\
chr1\tt\tCDS\t200\t221\t.\t+\t0\tID=c2;Parent=t1;protein_id=P1\n\
chr1\tt\tCDS\t300\t320\t.\t+\t2\tID=c3;Parent=t1;protein_id=P1\n",
        );
        assert_eq!(s.transcripts.len(), 1);
        let tx = &s.transcripts[0];
        assert_eq!(tx.protein_id.as_deref(), Some("P1"));
        assert_eq!(tx.cds_len, 15 + 22 + 21);
        assert_eq!(tx.cds_segments, 3);
        assert_eq!(tx.introns.len(), 2);
        assert_eq!((tx.introns[0].start, tx.introns[0].end), (115, 199));
        assert_eq!(
            tx.introns[0].cds_offset,
            Some(15),
            "15 coding bases precede the first intron"
        );
        assert_eq!(tx.introns[1].cds_offset, Some(37));
        assert_eq!(tx.introns[0].length(), 85);
        assert_eq!(
            tx.frame_mismatches,
            0,
            "consistent frames are not flagged: {:?}",
            s.diagnostics.retained()
        );
        let sites: Vec<_> = tx.protein_sites().map(|p| (p.codon(), p.phase)).collect();
        assert_eq!(sites, vec![(5, 0), (13, 1)]);
    }

    #[test]
    fn minus_strand_orders_from_high_to_low() {
        let s = load(
            "chr1\tt\tCDS\t300\t320\t.\t-\t0\tID=c1;Parent=t2;protein_id=P2\n\
chr1\tt\tCDS\t200\t214\t.\t-\t0\tID=c2;Parent=t2;protein_id=P2\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.strand, Strand::Minus);
        assert_eq!(tx.introns.len(), 1);
        assert_eq!((tx.introns[0].start, tx.introns[0].end), (215, 299));
        assert_eq!(
            tx.introns[0].cds_offset,
            Some(21),
            "300-320 comes first in transcript order"
        );
        assert_eq!(tx.frame_mismatches, 0);
    }

    #[test]
    fn utr_introns_are_kept_but_not_projectable() {
        // 5'UTR exon 50-89 (so 90-99 is a UTR intron), CDS exons 100-114 + 200-214, 3'UTR 300-330.
        let s = load(
            "chr1\tt\texon\t50\t89\t.\t+\t.\tID=e1;Parent=t3\n\
chr1\tt\texon\t100\t114\t.\t+\t.\tID=e2;Parent=t3\n\
chr1\tt\texon\t200\t214\t.\t+\t.\tID=e3;Parent=t3\n\
chr1\tt\texon\t300\t330\t.\t+\t.\tID=e4;Parent=t3\n\
chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t3;protein_id=P3\n\
chr1\tt\tCDS\t200\t214\t.\t+\t0\tID=c2;Parent=t3;protein_id=P3\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.introns.len(), 3, "UTR introns are still counted");
        assert_eq!(
            tx.introns[0].cds_offset, None,
            "5'UTR intron has no protein position"
        );
        assert_eq!(tx.introns[1].cds_offset, Some(15));
        assert_eq!(
            tx.introns[2].cds_offset, None,
            "3'UTR intron has no protein position"
        );
        assert_eq!(tx.n_coding_introns(), 1);
        assert_eq!(tx.introns[1].rank, 1);
    }

    #[test]
    fn flags_frame_inconsistency() {
        // c3 starts after 30 coding bases, so its frame must be 0; it claims 2.
        let s = load(
            "chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t4;protein_id=P4\n\
chr1\tt\tCDS\t200\t214\t.\t+\t0\tID=c2;Parent=t4;protein_id=P4\n\
chr1\tt\tCDS\t300\t314\t.\t+\t2\tID=c3;Parent=t4;protein_id=P4\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.frame_mismatches, 1);
        assert!(
            s.diagnostics.count_of("frame_mismatch") >= 1,
            "{:?}",
            s.diagnostics.tally()
        );
    }

    #[test]
    fn contiguous_cds_pieces_collapse_without_an_intron() {
        let s = load(
            "chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t5;protein_id=P5\n\
chr1\tt\tCDS\t115\t130\t.\t+\t0\tID=c2;Parent=t5;protein_id=P5\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(
            tx.introns.len(),
            0,
            "contiguous pieces are one coding segment"
        );
        assert_eq!(tx.cds_len, 31);
        assert_eq!(tx.cds_segments, 2, "raw feature count is preserved");
    }

    #[test]
    fn overlapping_cds_is_flagged_and_merged() {
        let s = load(
            "chr1\tt\tCDS\t100\t120\t.\t+\t0\tID=c1;Parent=ov;protein_id=PO\n\
chr1\tt\tCDS\t115\t132\t.\t+\t0\tID=c2;Parent=ov;protein_id=PO\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.cds_len, 33, "union 100-132");
        assert_eq!(tx.introns.len(), 0);
        assert!(s.diagnostics.count_of("overlapping_cds") >= 1);
    }

    #[test]
    fn gtf_quoted_attributes_and_explicit_introns() {
        let src = |body: &str| load_with(body, |o| o.keep_noncoding = true);
        let s = src(
            "chr1\thavana\texon\t100\t114\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n\
chr1\thavana\texon\t200\t214\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n",
        );
        assert_eq!(s.transcripts.len(), 1);
        let tx = &s.transcripts[0];
        assert_eq!(tx.gene_id.as_deref(), Some("G"));
        assert_eq!(tx.introns.len(), 1);
        assert_eq!(
            tx.n_coding_introns(),
            0,
            "no CDS features, so nothing is projectable"
        );
    }

    #[test]
    fn explicit_intron_features_are_used_when_exons_are_absent() {
        let s = load(
            "chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=ti;protein_id=PI\n\
chr1\tt\tCDS\t200\t214\t.\t+\t0\tID=c2;Parent=ti;protein_id=PI\n\
chr1\tt\tintron\t115\t199\t.\t+\t.\tID=i1;Parent=ti\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.introns.len(), 1);
        assert_eq!(tx.introns[0].cds_offset, Some(15));
        assert_eq!(
            s.diagnostics.count_of("intron_feature_mismatch"),
            0,
            "explicit intron matches"
        );
    }

    #[test]
    fn transcripts_without_cds_are_dropped_by_default() {
        let s = load("chr1\tt\texon\t100\t114\t.\t+\t.\tParent=txN\nchr1\tt\texon\t200\t214\t.\t+\t.\tParent=txN\n");
        assert!(s.transcripts.is_empty());
        assert!(s.diagnostics.count_of("no_cds") >= 1);
    }

    #[test]
    fn derive_cds_from_exons_rescues_protein2genome_style() {
        let s = load_with(
            "chr1\texonerate\texon\t100\t114\t.\t+\t.\tID=e1;Parent=p2g1\n\
             chr1\texonerate\texon\t200\t214\t.\t+\t.\tID=e2;Parent=p2g1\n",
            |o| o.derive_cds_from_exons = true,
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.cds_len, 30);
        assert_eq!(tx.introns[0].cds_offset, Some(15));
        assert!(s.diagnostics.count_of("cds_from_exons") >= 1);
    }

    #[test]
    fn namespace_stripping_and_version_control() {
        let s = load("chr1\tens\tCDS\t100\t114\t.\t+\t0\tID=cds:ENSC1;Parent=rna:ENST1.3;protein_id=ENSP9.1\n");
        let tx = &s.transcripts[0];
        assert_eq!(
            tx.transcript_id, "ENST1.3",
            "rna: prefix stripped, version kept by default"
        );
        assert_eq!(tx.protein_id.as_deref(), Some("ENSP9.1"));

        let s = load_with(
            "chr1\tens\tCDS\t100\t114\t.\t+\t0\tID=cds:ENSC1;Parent=rna:ENST1;protein_id=ENSP9.1\n",
            |o| o.strip_version = true,
        );
        assert_eq!(s.transcripts[0].protein_id.as_deref(), Some("ENSP9"));
    }

    #[test]
    fn id_map_rewrites_join_keys() {
        let mut opts = GffOptions::default();
        opts.id_map.insert("OLD".into(), "NEW".into());
        let s = load_with(
            "chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t6;protein_id=OLD\n",
            |o| *o = opts,
        );
        assert_eq!(s.transcripts[0].protein_id.as_deref(), Some("NEW"));
    }

    #[test]
    fn protein_child_feature_resolves_through_cds_id() {
        let s = load(
            "chr1\tmaker\tCDS\t100\t114\t.\t+\t0\tID=cds:g1.t1;Parent=g1.t1\n\
chr1\tmaker\texon\t100\t114\t.\t+\t.\tID=exon:g1.t1.1;Parent=g1.t1\n\
chr1\tmaker\tprotein\t100\t114\t.\t+\t.\tID=g1.p1;Parent=cds:g1.t1\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(
            tx.protein_id.as_deref(),
            Some("g1.p1"),
            "protein feature's Parent points at the CDS id"
        );
        assert_eq!(
            tx.gene_id.as_deref(),
            None,
            "no gene attribute was invented from the CDS id"
        );
    }

    #[test]
    fn transcript_fallback_when_no_protein_id_exists() {
        let s = load("chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t7\n");
        assert_eq!(s.transcripts[0].protein_id.as_deref(), Some("t7"));
        assert!(s.diagnostics.count_of("protein_id_fallback_transcript") >= 1);
    }

    #[test]
    fn missing_transcript_key_is_reported_not_fatal() {
        let s = load("chr1\tt\tCDS\t100\t114\t.\t+\t0\tID=c1;note=orphan\n");
        assert!(s.transcripts.is_empty());
        assert_eq!(s.diagnostics.count_of("no_transcript_key"), 1);
    }

    #[test]
    fn two_files_do_not_silently_merge_distinct_transcripts() {
        let mk = |body: &str| {
            let dir = std::env::temp_dir().join(format!(
                "ip_multi_{}_{}",
                std::process::id(),
                ANCHOR.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            ));
            std::fs::create_dir_all(&dir).unwrap();
            let p = dir.join(format!("a-{}.gff3", std::process::id()));
            std::fs::write(&p, body).unwrap();
            p
        };
        let mut src = GffSource::new(
            vec![
                mk("c1\tt\tCDS\t100\t114\t.\t+\t0\tID=x;Parent=dup;protein_id=PA\n"),
                mk("c1\tt\tCDS\t200\t214\t.\t+\t0\tID=y;Parent=dup;protein_id=PB\n"),
            ],
            GffOptions::default(),
        );
        let s = src.load().unwrap();
        // Files are independent annotations (one species each): the same id in two of them is two
        // transcripts, never one chimera built from both, and the repeat is reported.
        assert_eq!(s.transcripts.len(), 2);
        assert!(s
            .transcripts
            .iter()
            .all(|t| t.cds_len == 15 && t.introns.is_empty()));
        let prots: Vec<_> = s
            .transcripts
            .iter()
            .map(|t| t.protein_id.clone().unwrap())
            .collect();
        assert_eq!(
            prots,
            vec!["PA".to_string(), "PB".to_string()],
            "file order is kept"
        );
        assert_eq!(
            s.diagnostics.count_of("transcript_id_collision"),
            1,
            "{:?}",
            s.diagnostics.tally()
        );
    }

    #[test]
    fn mrna_row_supplies_gene_and_strand_to_child_cds() {
        let s = load(
            "chr1\tt\tgene\t50\t330\t.\t+\t.\tID=g1;gene_id=GENE-A\n\
chr1\tt\tmRNA\t50\t330\t.\t+\t.\tID=t8;Parent=g1;protein_id=P8\n\
chr1\tt\tCDS\t100\t114\t.\t+\t0\tParent=t8\n\
chr1\tt\tCDS\t200\t214\t.\t+\t0\tParent=t8\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(
            tx.gene_id.as_deref(),
            Some("g1"),
            "mRNA Parent is the gene id"
        );
        assert_eq!(
            tx.protein_id.as_deref(),
            Some("P8"),
            "mRNA protein_id reaches the CDS rows"
        );
        assert_eq!(tx.strand, Strand::Plus);
        assert_eq!(tx.introns[0].cds_offset, Some(15));
    }

    #[test]
    fn five_prime_partial_cds_offsets_start_at_first_full_codon() {
        // First CDS (10 nt) has phase 1: its first base belongs to a codon cut off at the 5' end, so
        // 9 coding bases (3 full codons) precede the intron -> offset 9, phase 0. The next feature
        // then has phase 0 and must not be flagged.
        let s = load(
            "c\tt\tCDS\t1\t10\t.\t+\t1\tParent=pp;protein_id=PP\n\
c\tt\tCDS\t101\t112\t.\t+\t0\tParent=pp;protein_id=PP\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.introns[0].cds_offset, Some(9));
        assert_eq!(tx.frame_mismatches, 0, "{:?}", s.diagnostics.retained());
        assert_eq!(tx.cds_len, 21, "22 CDS bases minus 1 skipped");
    }

    #[test]
    fn five_prime_partial_cds_with_phase_two_on_minus_strand() {
        // Transcript order is 112..101 (12 nt, phase 2 = skip 2) then 10..1: offset 10, and the next
        // feature must carry phase (3 - 10 % 3) % 3 = 2.
        let s = load(
            "c\tt\tCDS\t101\t112\t.\t-\t2\tParent=pm;protein_id=PM\n\
c\tt\tCDS\t1\t10\t.\t-\t2\tParent=pm;protein_id=PM\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.introns[0].cds_offset, Some(10));
        assert_eq!(tx.frame_mismatches, 0, "{:?}", s.diagnostics.retained());
    }

    #[test]
    fn multi_parent_feature_joins_every_parent_and_makes_no_phantom_transcript() {
        let s = load_with(
            "c\tt\tmRNA\t1\t300\t.\t+\t.\tID=t1\n\
c\tt\tmRNA\t1\t300\t.\t+\t.\tID=t2\n\
c\tt\tCDS\t1\t12\t.\t+\t0\tParent=t1,t2;protein_id=P\n\
c\tt\tCDS\t101\t112\t.\t+\t0\tParent=t1;protein_id=P\n\
c\tt\tCDS\t101\t115\t.\t+\t0\tParent=t2;protein_id=P\n",
            |o| o.keep_noncoding = true,
        );
        let ids: Vec<&str> = s
            .transcripts
            .iter()
            .map(|t| t.transcript_id.as_str())
            .collect();
        assert_eq!(ids, vec!["t1", "t2"]);
        assert!(s.transcripts.iter().all(|t| t.introns.len() == 1));
    }

    #[test]
    fn ncbi_exon_rows_with_transcript_id_join_cds_rows_with_parent_only() {
        let s = load(
            "c\tR\tmRNA\t1\t400\t.\t+\t.\tID=rna-NM_1.1;Parent=gene-A;transcript_id=NM_1.1\n\
c\tR\texon\t1\t50\t.\t+\t.\tID=e1;Parent=rna-NM_1.1;transcript_id=NM_1.1\n\
c\tR\texon\t101\t200\t.\t+\t.\tID=e2;Parent=rna-NM_1.1;transcript_id=NM_1.1\n\
c\tR\texon\t301\t400\t.\t+\t.\tID=e3;Parent=rna-NM_1.1;transcript_id=NM_1.1\n\
c\tR\tCDS\t110\t200\t.\t+\t0\tID=cds-NP_1.1;Parent=rna-NM_1.1;protein_id=NP_1.1\n\
c\tR\tCDS\t301\t320\t.\t+\t2\tID=cds-NP_1.1;Parent=rna-NM_1.1;protein_id=NP_1.1\n",
        );
        assert_eq!(s.transcripts.len(), 1, "{:?}", s.diagnostics.tally());
        let tx = &s.transcripts[0];
        assert_eq!(
            tx.introns.len(),
            2,
            "one 5' UTR intron plus one coding intron"
        );
        assert_eq!(tx.n_coding_introns(), 1);
    }

    #[test]
    fn diagnostic_counts_are_not_truncated() {
        let mut body = String::new();
        for i in 0..10 {
            body.push_str(&format!(
                "c\tt\tCDS\t{}\t{}\t.\t+\t0\tID=x{i}\n",
                1 + i * 20,
                12 + i * 20
            ));
        }
        let s = load(&body);
        assert_eq!(s.diagnostics.count_of("no_transcript_key"), 10);
    }

    #[test]
    fn latin1_bytes_do_not_abort_the_parse() {
        let mut body = b"c\tt\tCDS\t1\t12\t.\t+\t0\tParent=t1;protein_id=P;Note=caf".to_vec();
        body.push(0xe9);
        body.extend_from_slice(b"\nc\tt\tCDS\t101\t112\t.\t+\t0\tParent=t1;protein_id=P\n");
        let dir = std::env::temp_dir().join(format!("ip_latin1_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("in.gff3");
        std::fs::write(&p, &body).unwrap();
        let s = GffSource::path(p).load().expect("lossy read, not an error");
        assert_eq!(s.transcripts[0].introns.len(), 1);
        assert_eq!(s.diagnostics.count_of("invalid_utf8"), 1);
    }

    #[test]
    fn quoted_gff3_values_are_unquoted() {
        let s = load_with(
            "c\tt\tCDS\t1\t12\t.\t+\t0\tParent=t1;protein_id=\"P1\"\n",
            |o| o.keep_noncoding = true,
        );
        assert_eq!(s.transcripts[0].protein_id.as_deref(), Some("P1"));
    }

    #[test]
    fn strip_version_is_applied_consistently_to_annotation_and_sequences() {
        let s = load_with(
            "c\tt\tCDS\t1\t12\t.\t+\t0\tParent=t1;protein_id=A.1.2\n",
            |o| {
                o.strip_version = true;
                o.keep_noncoding = true;
            },
        );
        let pid = s.transcripts[0].protein_id.clone().unwrap();
        let mut ps = crate::seqio::ProteinSet::new();
        ps.read_one(
            std::io::Cursor::new(b">A.1.2\nMKVL\n"),
            std::path::Path::new("m"),
        )
        .unwrap();
        assert_eq!(ps.lookup(&pid), Some(0));
    }

    #[test]
    fn expected_frame_is_bases_to_skip_not_bases_emitted() {
        assert_eq!(
            (0..7i64).map(expected_frame).collect::<Vec<_>>(),
            vec![0, 2, 1, 0, 2, 1, 0]
        );
    }

    #[test]
    fn spec_phases_are_accepted_when_exons_are_not_multiples_of_three() {
        // exon 1 = 13 nt leaves 1 base of a codon, so exon 2 has phase 2; 25 coding bases then
        // precede exon 3, and (3 - 25 % 3) % 3 = 2.
        let s = load(
            "chr1\tt\tCDS\t100\t112\t.\t+\t0\tID=c1;Parent=np;protein_id=NP1\n\
chr1\tt\tCDS\t200\t211\t.\t+\t2\tID=c2;Parent=np;protein_id=NP1\n\
chr1\tt\tCDS\t300\t314\t.\t+\t2\tID=c3;Parent=np;protein_id=NP1\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.frame_mismatches, 0, "{:?}", s.diagnostics.retained());
        assert_eq!(tx.introns[0].cds_offset, Some(13));
        assert_eq!(tx.introns[1].cds_offset, Some(25));
    }

    #[test]
    fn emitted_bases_mod_three_is_flagged_as_a_mismatch() {
        // The old (wrong) reading would call phase 1 correct after 13 coding bases.
        let s = load(
            "chr1\tt\tCDS\t100\t112\t.\t+\t0\tID=c1;Parent=nq;protein_id=NQ1\n\
chr1\tt\tCDS\t200\t211\t.\t+\t1\tID=c2;Parent=nq;protein_id=NQ1\n",
        );
        assert_eq!(s.transcripts[0].frame_mismatches, 1);
    }

    #[test]
    fn strand_inference_uses_spec_phases() {
        // Minus-strand transcript, no strand column. Transcript order is 320-300 (21 nt... use 22 nt:
        // 299-320) then 200-211, so the second piece follows 22 coding bases -> phase 2.
        let s = load(
            "chr1\tt\tCDS\t299\t320\t.\t.\t0\tID=c1;Parent=ms;protein_id=MS\n\
chr1\tt\tCDS\t200\t211\t.\t.\t2\tID=c2;Parent=ms;protein_id=MS\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(tx.introns[0].cds_offset, Some(22));
        assert_eq!(tx.frame_mismatches, 0);
    }

    #[test]
    fn strand_is_inferred_from_frames_when_absent() {
        // Two CDS pieces, no usable strand, frames consistent only with descending order.
        let s = load(
            "chr1\tt\tCDS\t300\t320\t.\t.\t0\tID=c1;Parent=ts;protein_id=PS\n\
chr1\tt\tCDS\t200\t209\t.\t.\t.\tID=c2;Parent=ts;protein_id=PS\n",
        );
        let tx = &s.transcripts[0];
        assert_eq!(
            tx.introns[0].cds_offset,
            Some(21),
            "300-320 first implies offset 21"
        );
        assert!(s.diagnostics.count_of("strand_unknown") >= 1);
    }

    #[test]
    fn many_transcripts_project_independently() {
        let mut body = String::new();
        for i in 0..200u32 {
            body.push_str(&format!(
                "c1\tt\tCDS\t100\t114\t.\t+\t0\tID=c{i}a;Parent=t{i};protein_id=P{i}\n"
            ));
            body.push_str(&format!(
                "c1\tt\tCDS\t200\t{}\t.\t+\t0\tID=c{i}b;Parent=t{i};protein_id=P{i}\n",
                200 + i
            ));
        }
        let s = load(&body);
        assert_eq!(s.transcripts.len(), 200);
        let t150 = s
            .transcripts
            .iter()
            .find(|t| t.protein_id.as_deref() == Some("P150"))
            .unwrap();
        assert_eq!(t150.cds_len, 15 + 151);
        assert_eq!(t150.introns[0].cds_offset, Some(15));
        assert_eq!(t150.introns[0].end, 199);
    }

    // ---- proposed: random CDS layouts vs a naive per-base model, both strands ----
    #[test]
    fn random_cds_layouts_match_a_naive_model_on_both_strands() {
        let mut s: u64 = 0x2545_F491_4F6C_DD1D;
        let mut rnd = move |n: u64| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s % n
        };
        for case in 0..300 {
            let n_seg = 1 + rnd(5) as usize;
            let minus = rnd(2) == 1;
            let skip = rnd(3); // 5'-partial: bases before the first full codon
            let mut segs = Vec::new();
            let mut pos = 1 + rnd(50);
            for _ in 0..n_seg {
                let len = 3 + rnd(20); // >= 3 so the first segment always covers the skip
                segs.push((pos, pos + len - 1));
                pos += len + 1 + rnd(30);
            }
            let tx_order: Vec<(u64, u64)> = if minus {
                segs.iter().rev().copied().collect()
            } else {
                segs.clone()
            };
            let total: u64 = segs.iter().map(|(a, b)| b - a + 1).sum();
            let mut body = String::new();
            let mut before = 0i64;
            for (k, &(a, b)) in tx_order.iter().enumerate() {
                // GFF3 phase = bases to remove to reach the next codon start.
                let phase = if k == 0 {
                    skip as i64
                } else {
                    (3 - (before - skip as i64).rem_euclid(3)) % 3
                };
                let st = if minus { '-' } else { '+' };
                body += &format!(
                    "c\tt\tCDS\t{a}\t{b}\t.\t{st}\t{phase}\tParent=r{case};protein_id=P{case}\n"
                );
                before += (b - a + 1) as i64;
            }
            let set = load(&body);
            let tx = &set.transcripts[0];
            let mut want = Vec::new();
            let mut cum = 0u64;
            for w in tx_order.windows(2) {
                cum += w[0].1 - w[0].0 + 1;
                let (lo, hi) = if minus {
                    (w[1].1 + 1, w[0].0 - 1)
                } else {
                    (w[0].1 + 1, w[1].0 - 1)
                };
                want.push((lo, hi, Some(cum - skip)));
            }
            let got: Vec<_> = tx
                .introns
                .iter()
                .map(|i| (i.start, i.end, i.cds_offset.map(u64::from)))
                .collect();
            assert_eq!(got, want, "case {case}:\n{body}");
            assert_eq!(tx.cds_len, total - skip, "case {case}");
            assert_eq!(tx.frame_mismatches, 0, "case {case}:\n{body}");
            assert_eq!(tx.strand, if minus { Strand::Minus } else { Strand::Plus });
        }
    }

    #[test]
    fn many_files_read_in_parallel_equal_one_concatenated_file() {
        let dir = std::env::temp_dir().join(format!("ip_par_files_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut paths = Vec::new();
        let mut all = String::new();
        for f in 0..9u32 {
            let mut body = String::new();
            for g in 0..25u32 {
                let id = format!("f{f}g{g}");
                let strand = if (f + g) % 2 == 0 { '+' } else { '-' };
                // Two CDS pieces of 15 + (g % 7 + 1) * 3 nt, so offsets differ between genes.
                let len2 = 3 * (g % 7 + 1);
                let (a1, b1, a2, b2) = (100, 114, 200, 199 + len2);
                for (a, b, ph) in [(a1, b1, 0), (a2, b2, 0)] {
                    body += &format!("c{f}\tt\tCDS\t{a}\t{b}\t.\t{strand}\t{ph}\tID={id}.c;Parent={id};protein_id=P{id}\n");
                }
            }
            all.push_str(&body);
            let p = dir.join(format!("part{f}.gff3"));
            std::fs::write(&p, body).unwrap();
            paths.push(p);
        }
        let one = dir.join("all.gff3");
        std::fs::write(&one, all).unwrap();
        let multi = GffSource::new(paths, GffOptions::default()).load().unwrap();
        let single = GffSource::path(one).load().unwrap();
        let key = |t: &TranscriptAnnotation| {
            (
                t.transcript_id.clone(),
                t.protein_id.clone(),
                t.strand.to_string(),
                t.cds_len,
                t.introns
                    .iter()
                    .map(|i| (i.start, i.end, i.cds_offset))
                    .collect::<Vec<_>>(),
            )
        };
        let mut a: Vec<_> = multi.transcripts.iter().map(key).collect();
        let mut b: Vec<_> = single.transcripts.iter().map(key).collect();
        a.sort();
        b.sort();
        assert_eq!(a.len(), 9 * 25);
        assert_eq!(
            a, b,
            "parallel per-file reading must not change any transcript"
        );
        assert_eq!(multi.diagnostics.count_of("transcript_id_collision"), 0);
        // Deterministic: file order, then transcript id within a file.
        let ids: Vec<&str> = multi
            .transcripts
            .iter()
            .map(|t| t.transcript_id.as_str())
            .collect();
        let again = GffSource::new(
            (0..9).map(|f| dir.join(format!("part{f}.gff3"))).collect(),
            GffOptions::default(),
        )
        .load()
        .unwrap();
        let ids2: Vec<&str> = again
            .transcripts
            .iter()
            .map(|t| t.transcript_id.as_str())
            .collect();
        assert_eq!(ids, ids2);
    }
}
