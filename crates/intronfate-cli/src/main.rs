//! `intronfate` — project annotated intron positions into protein alignments and score them.
//!
//! Three composable steps, so an expensive projection runs once and cheap re-scoring can sweep
//! parameters or topologies:
//!
//! ```text
//! intronfate project --gff a.gff.gz -p a.faa                    # introns -> protein space
//! intronfate matrix  --gff a.gff.gz -p a.faa --msa family.faa   # protein space -> columns
//! intronfate score   --matrix out/family/matrix.tsv -t tree.nwk # columns -> Dollo / MP
//! intronfate analyze --gff a.gff.gz -p a.faa --msa family.faa -t tree.nwk
//! ```

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use intronfate_core::annotation::{AnnotationSet, AnnotationSource, GffOptions, GffSource};
use intronfate_core::matrix::SiteMatrix;
use intronfate_core::msa::Alignment;
use intronfate_core::project::{
    join_annotations, project_into_alignment, AlignmentProjection, JoinOptions, JoiningSink,
    PhaseMode, ProjectOptions, Projection, TailPolicy,
};
use intronfate_core::report::Diagnostics;
use intronfate_core::score::{self, Method, PrepareReport, ScoringData};
use intronfate_core::seqio::ProteinSet;
use intronfate_core::tables;
use intronfate_core::tree::Tree;

#[derive(Parser, Debug)]
#[command(
    name = "intronfate",
    version,
    about = "Project intron positions from annotation into protein alignments, then score how they are shared.",
    long_about = "Project intron positions from GFF3/GTF annotation into protein space, map them onto \
                  MSA columns for a set of orthologs, and score which positions are shared using Dollo \
                  parsimony, maximum parsimony, or plain counts.\n\n\
                  A position is reported as cds_offset (coding bases 5' of the intron), codon, and \
                  phase (0/1/2), plus the alignment column it lands in. Two introns count as the same \
                  site when they share a column, and — unless --ignore-phase is given — a phase too."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Print timing and per-stage record counts.
    #[arg(short = 'v', long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Write protein-space intron positions for every annotated transcript.
    Project(ProjectArgs),
    /// Project intron positions onto alignments and write site x taxon matrices.
    Matrix(MatrixArgs),
    /// Score an existing matrix with Dollo parsimony, maximum parsimony and/or counts.
    Score(ScoreArgs),
    /// Run project + matrix + score in one pass.
    Analyze(Box<AnalyzeArgs>),
    /// Convert between this tool's matrix.tsv files and MALIN intron tables.
    Convert(ConvertArgs),
    /// Run the whole pipeline on a built-in toy dataset and print the tables.
    Demo(DemoArgs),
}

#[derive(Args, Debug)]
struct AnnotationArgs {
    /// GFF3/GTF file(s). Repeatable, e.g. one per species.
    #[arg(long, value_name = "FILE")]
    gff: Vec<PathBuf>,

    /// Text file listing GFF3/GTF files, one per line (`#` comments; relative paths are relative to the list). Repeatable.
    #[arg(long, value_name = "LIST")]
    gff_list: Vec<PathBuf>,

    /// Folder of GFF3/GTF files (matched by `--gff-ext`, default gff3, gff, gtf; `.gz` allowed). Repeatable.
    #[arg(long, value_name = "DIR")]
    gff_dir: Vec<PathBuf>,

    /// Extension(s) selecting files in `--gff-dir`, e.g. `gff3` or `annot.gff`.
    #[arg(long, value_delimiter = ',', value_name = "EXT")]
    gff_ext: Vec<String>,

    /// Protein FASTA file(s) corresponding to the annotation. Repeatable.
    #[arg(long, short = 'p', value_name = "FILE")]
    proteins: Vec<PathBuf>,

    /// Text file listing protein FASTA files, one per line. Repeatable.
    #[arg(long, value_name = "LIST")]
    proteins_list: Vec<PathBuf>,

    /// Folder of protein FASTA files (matched by `--proteins-ext`, default faa, fa, fasta, pep; `.gz` allowed). Repeatable.
    #[arg(long, value_name = "DIR")]
    proteins_dir: Vec<PathBuf>,

    /// Extension(s) selecting files in `--proteins-dir`.
    #[arg(long, value_delimiter = ',', value_name = "EXT")]
    proteins_ext: Vec<String>,

    /// Attribute(s) naming the transcript a CDS/exon belongs to, tried in order.
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "Parent,transcript_id,transcript,rna_id"
    )]
    transcript_attrs: Vec<String>,

    /// Attribute(s) naming the protein a transcript encodes, tried in order.
    #[arg(long, value_delimiter = ',', default_value = "protein_id,protein")]
    protein_attrs: Vec<String>,

    /// Attribute(s) naming the parent gene, tried in order.
    #[arg(long, value_delimiter = ',', default_value = "gene_id,gene,locus_tag")]
    gene_attrs: Vec<String>,

    /// Feature type(s) to treat as CDS.
    #[arg(long, value_delimiter = ',', default_value = "cds")]
    cds_types: Vec<String>,

    /// Two-column `annotation_id -> protein_id` mapping file (TSV or CSV), for stubborn joins.
    #[arg(long, value_name = "FILE")]
    id_map: Option<PathBuf>,

    /// Ignore trailing `.N` versions when matching ids.
    #[arg(long)]
    strip_version: bool,

    /// Treat `exon` features as coding where a transcript has no CDS (protein2genome-style GFF).
    #[arg(long)]
    derive_cds_from_exons: bool,

    /// Keep transcripts with no CDS in the report; their introns are never projectable.
    #[arg(long)]
    keep_noncoding: bool,

    /// Abort on a malformed annotation line instead of counting and skipping it.
    #[arg(long)]
    strict: bool,
}

#[derive(Args, Debug)]
struct ProjectArgs {
    #[command(flatten)]
    ann: AnnotationArgs,

    /// Directory for output tables.
    #[arg(long, short = 'o', default_value = "intron_out")]
    out_dir: PathBuf,

    /// Tolerated difference between protein length and CDS/3, in residues. `none` disables the check.
    #[arg(long, default_value = "3")]
    length_slack: String,
}

#[derive(Args, Debug)]
struct MatrixArgs {
    #[command(flatten)]
    ann: AnnotationArgs,

    /// Protein alignment(s) of orthologs: aligned FASTA or Stockholm. Repeatable.
    #[arg(long, value_name = "FILE")]
    msa: Vec<PathBuf>,

    /// Text file listing alignments, one per line (`#` comments; relative paths are relative to the list). Repeatable.
    #[arg(long, value_name = "LIST")]
    msa_list: Vec<PathBuf>,

    /// Folder of alignments (matched by `--msa-ext`, default aln, aln.fasta, aln.fa, aln.faa, afa, msa, sto, stockholm; `.gz` allowed; plain `.fasta` is not a default because proteins use it too). Repeatable.
    #[arg(long, value_name = "DIR")]
    msa_dir: Vec<PathBuf>,

    /// Extension(s) selecting files in `--msa-dir`, e.g. `aln.fasta` (a name matches when it is, or ends in, `.ext`).
    #[arg(long, value_delimiter = ',', value_name = "EXT")]
    msa_ext: Vec<String>,

    #[arg(long, short = 'o', default_value = "intron_out")]
    out_dir: PathBuf,

    /// Also write one matrix that stacks the sites of every alignment (`<out-dir>/combined/matrix.tsv`, with `family` and
    /// `family_column` columns); `analyze` scores it too, so the totals over all families need no summing afterwards.
    #[arg(long)]
    combine: bool,

    /// With `--combine`: write only `<out-dir>/combined/` (plus the run-level `introns.tsv`, `unmatched.tsv` and diagnostics), not one
    /// folder per alignment. For runs over thousands of families, where the per-family folders are most of the output.
    #[arg(long, requires = "combine")]
    combined_only: bool,

    /// Group sites by alignment column only, ignoring codon phase.
    #[arg(long)]
    ignore_phase: bool,

    /// Place sites whose codon falls outside an alignment row on that row's last column.
    #[arg(long)]
    clamp_tail: bool,

    /// Drift tolerance in columns for the pairwise sharing statistics.
    #[arg(long, default_value = "0")]
    tolerance: u32,

    /// Two-column `alignment_row_id -> taxon` file (TSV or CSV). Needed when one alignment per gene
    /// family holds a differently named protein for each species, so the tree's tips (species) cannot
    /// equal the row ids. Rows not listed keep their own id; two rows of one alignment may not map to
    /// the same taxon.
    #[arg(long, value_name = "FILE")]
    taxon_map: Option<PathBuf>,

    /// Residues a taxon must have on each side of a site before its absence counts as an observation
    /// (otherwise the cell is `?`, missing data, and costs no loss). 0 = only the intron's own
    /// residue. Phase-0 introns always need one residue on each side.
    #[arg(long, default_value = "1")]
    flank: u32,

    #[arg(long, default_value = "3")]
    length_slack: String,
}

#[derive(Args, Debug)]
struct AnalyzeArgs {
    #[command(flatten)]
    matrix: MatrixArgs,

    /// Newick tree whose tips are matrix taxon names.
    #[arg(long, short = 't', value_name = "FILE")]
    tree: Option<PathBuf>,

    /// Scoring method(s): dollo, mp, counts. Defaults to all three.
    #[arg(long, short = 'm', value_delimiter = ',')]
    methods: Vec<String>,

    /// Outgroup tip name(s) used to root the tree. Repeatable.
    #[arg(long = "outgroup", value_name = "NAME")]
    outgroups: Vec<String>,

    /// Ignore sites carried by fewer than this many taxa.
    #[arg(long, default_value = "1")]
    min_taxa: usize,

    /// Ignore sites observed (present or absent, not `?`) in fewer than this many taxa.
    #[arg(long, default_value = "0")]
    min_observed: usize,

    /// Write per-node ancestral states.
    #[arg(long)]
    reconstructions: bool,
}

#[derive(Args, Debug)]
struct ScoreArgs {
    /// `matrix.tsv` written by the matrix or analyze command. One file scores into `--out-dir` itself; with several
    /// (or any of the list/folder forms) each is scored into `<out-dir>/<name>/`. Repeatable.
    #[arg(
        long,
        value_name = "FILE",
        required_unless_present_any = ["matrix_list", "matrix_dir"]
    )]
    matrix: Vec<PathBuf>,

    /// Text file listing matrix files, one per line. Repeatable.
    #[arg(long, value_name = "LIST")]
    matrix_list: Vec<PathBuf>,

    /// Folder of matrices: the folder and its immediate sub-folders are searched, so the output folder of
    /// `analyze`/`matrix` works as is (`<dir>/<family>/matrix.tsv`). Repeatable.
    #[arg(long, value_name = "DIR")]
    matrix_dir: Vec<PathBuf>,

    /// Extension(s) selecting files in `--matrix-dir` (default `matrix.tsv`, which matches `matrix.tsv` and `X.matrix.tsv`).
    #[arg(long, value_delimiter = ',', value_name = "EXT")]
    matrix_ext: Vec<String>,

    /// Stack all the given matrices into one combined matrix and score that once, writing everything into `-o`
    /// (sites are labelled `<alignment>:<column>:<phase>`; taxa missing from an alignment are `?` there).
    #[arg(long)]
    combine: bool,

    /// Newick tree whose tips are matrix taxon names.
    #[arg(long, short = 't', value_name = "FILE")]
    tree: Option<PathBuf>,

    /// Scoring method(s): dollo, mp, counts. Defaults to all three.
    #[arg(long, short = 'm', value_delimiter = ',')]
    methods: Vec<String>,

    /// Outgroup tip name(s) used to root the tree. Repeatable.
    #[arg(long = "outgroup", value_name = "NAME")]
    outgroups: Vec<String>,

    /// Ignore sites carried by fewer than this many taxa.
    #[arg(long, default_value = "1")]
    min_taxa: usize,

    /// Ignore sites observed (present or absent, not `?`) in fewer than this many taxa.
    #[arg(long, default_value = "0")]
    min_observed: usize,

    /// Drift tolerance in columns for the pairwise sharing statistics.
    #[arg(long, default_value = "0")]
    tolerance: u32,

    /// Write per-node ancestral states.
    #[arg(long)]
    reconstructions: bool,

    #[arg(long, short = 'o', default_value = "intron_out")]
    out_dir: PathBuf,
}

#[derive(Args, Debug)]
struct ConvertArgs {
    /// Import a MALIN table. Writes one `<alignment>/matrix.tsv` per alignment named in its `#MAP`
    /// lines (or a single `<table-name>/matrix.tsv` for a 0/1 table without them) under `--out`.
    #[arg(long, value_name = "TABLE", conflicts_with = "to_malin")]
    from_malin: Option<PathBuf>,

    /// Export matrix.tsv file(s) as one MALIN table with `#MAP` lines. All matrices must share the
    /// same taxa, in the same order, and must not have been built with --ignore-phase.
    #[arg(long, value_name = "MATRIX", num_args = 1.., conflicts_with = "from_malin")]
    to_malin: Vec<PathBuf>,

    /// Export: text file listing matrix files, one per line (`#` comments; relative paths are relative to the list). Repeatable.
    #[arg(long, value_name = "LIST", conflicts_with = "from_malin")]
    matrix_list: Vec<PathBuf>,

    /// Export: folder of matrices; the folder and its immediate sub-folders are searched (so `analyze`'s output folder works; a
    /// sub-folder named `combined` is skipped). Repeatable.
    #[arg(long, value_name = "DIR", conflicts_with = "from_malin")]
    matrix_dir: Vec<PathBuf>,

    /// Export: extension(s) selecting files in `--matrix-dir` (default `matrix.tsv`).
    #[arg(
        long,
        value_delimiter = ',',
        value_name = "EXT",
        conflicts_with = "from_malin"
    )]
    matrix_ext: Vec<String>,

    /// Output directory (import) or output file (export).
    #[arg(long, short = 'o')]
    out: PathBuf,
}

#[derive(Args, Debug)]
struct DemoArgs {
    #[arg(long, short = 'o', default_value = "intron_demo")]
    out_dir: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let pool = rayon::ThreadPoolBuilder::new()
        .build()
        .context("could not start the worker pool")?;
    let started = Instant::now();
    let verbose = cli.verbose;
    let result = pool.install(|| run(cli));
    match result {
        Ok(()) => {
            if verbose > 0 {
                eprintln!("done in {:.2}s", started.elapsed().as_secs_f64());
            }
            Ok(())
        }
        Err(e) => {
            eprintln!("error: {e:#}");
            // 1 = the run failed; clap itself exits 2 for a usage error, so the two stay distinguishable.
            std::process::exit(1);
        }
    }
}

fn run(cli: Cli) -> Result<()> {
    let verbose = cli.verbose;
    match cli.command {
        Command::Project(mut a) => {
            a.ann.resolve_inputs()?;
            cmd_project(&a, verbose)
        }
        Command::Matrix(mut a) => {
            a.resolve_inputs()?;
            cmd_matrix(&a, verbose)
        }
        Command::Analyze(mut a) => {
            a.matrix.resolve_inputs()?;
            cmd_analyze(&a, verbose)
        }
        Command::Score(mut a) => {
            a.resolve_inputs()?;
            cmd_score(&a, verbose)
        }
        Command::Demo(a) => cmd_demo(&a, verbose),
        Command::Convert(mut a) => {
            a.to_malin = expand_inputs(
                "matrix",
                (&a.to_malin, &a.matrix_list, &a.matrix_dir),
                &a.matrix_ext,
                &["matrix.tsv"],
                (1, &["combined"]),
            )?;
            cmd_convert(&a)
        }
    }
}

/// Where an input's files come from: `--x`, `--x-list`, `--x-dir`.
type Sources<'a> = (&'a [PathBuf], &'a [PathBuf], &'a [PathBuf]);

/// Fold `--x`, `--x-list` and `--x-dir/--x-ext` into one ordered file list (see `intronfate_core::inputs`).
fn expand_inputs(
    what: &str,
    (files, lists, dirs): Sources,
    exts: &[String],
    default_exts: &[&str],
    depth_and_skip: (usize, &[&str]),
) -> Result<Vec<PathBuf>> {
    Ok(intronfate_core::inputs::expand(
        &intronfate_core::inputs::Inputs {
            what,
            files,
            lists,
            dirs,
            exts,
            default_exts,
            depth: depth_and_skip.0,
            skip_dirs: depth_and_skip.1,
        },
    )?)
}

impl AnnotationArgs {
    fn resolve_inputs(&mut self) -> Result<()> {
        self.gff = expand_inputs(
            "GFF",
            (&self.gff, &self.gff_list, &self.gff_dir),
            &self.gff_ext,
            &["gff3", "gff", "gtf"],
            (0, &[]),
        )?;
        self.proteins = expand_inputs(
            "protein",
            (&self.proteins, &self.proteins_list, &self.proteins_dir),
            &self.proteins_ext,
            &["faa", "fa", "fasta", "pep"],
            (0, &[]),
        )?;
        Ok(())
    }
}

impl MatrixArgs {
    fn resolve_inputs(&mut self) -> Result<()> {
        self.ann.resolve_inputs()?;
        self.msa = expand_inputs(
            "alignment",
            (&self.msa, &self.msa_list, &self.msa_dir),
            &self.msa_ext,
            &[
                "aln",
                "aln.fasta",
                "aln.fa",
                "aln.faa",
                "afa",
                "msa",
                "sto",
                "stockholm",
            ],
            (0, &[]),
        )?;
        Ok(())
    }
}

impl ScoreArgs {
    /// `true` when the matrices came from anything but a single `--matrix FILE`: each then gets its own output folder.
    fn many(&self) -> bool {
        self.matrix.len() > 1 || !self.matrix_list.is_empty() || !self.matrix_dir.is_empty()
    }

    fn resolve_inputs(&mut self) -> Result<()> {
        self.matrix = expand_inputs(
            "matrix",
            (&self.matrix, &self.matrix_list, &self.matrix_dir),
            &self.matrix_ext,
            &["matrix.tsv"],
            (1, &["combined"]),
        )?;
        Ok(())
    }
}

/// Output-folder name for a matrix file: the family folder for `<family>/matrix.tsv`, else the file name without
/// `.matrix.tsv`/`.tsv`.
fn matrix_stem(p: &Path) -> String {
    let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("matrix");
    let name = name.strip_suffix(".gz").unwrap_or(name);
    if name == "matrix.tsv" {
        if let Some(parent) = p
            .parent()
            .and_then(|d| d.file_name())
            .and_then(|n| n.to_str())
        {
            return parent.to_string();
        }
    }
    name.strip_suffix(".matrix.tsv")
        .or_else(|| name.strip_suffix(".tsv"))
        .unwrap_or(name)
        .to_string()
}

/// Everything `AnnotationArgs::load_projection` returns.
struct Loaded {
    proj: Projection,
    proteins: ProteinSet,
    /// Diagnostics raised while reading and grouping the annotation.
    annotation_diagnostics: Diagnostics,
    n_transcripts: usize,
}

impl AnnotationArgs {
    fn gff_options(&self) -> Result<GffOptions> {
        if self.gff.is_empty() {
            bail!("--gff is required");
        }
        if self.proteins.is_empty() {
            bail!("--proteins is required");
        }
        let mut o = GffOptions {
            cds_types: self.cds_types.clone(),
            transcript_attrs: self.transcript_attrs.clone(),
            protein_attrs: self.protein_attrs.clone(),
            gene_attrs: self.gene_attrs.clone(),
            strip_version: self.strip_version,
            keep_noncoding: self.keep_noncoding,
            derive_cds_from_exons: self.derive_cds_from_exons,
            strict: self.strict,
            ..Default::default()
        };
        if let Some(path) = &self.id_map {
            o.id_map = GffSource::load_id_map(path)
                .with_context(|| format!("could not read id map {}", path.display()))?;
        }
        Ok(o)
    }

    /// Read the proteins, then stream the annotation straight into the join: each transcript is joined
    /// to its protein as it is read and then dropped, so memory follows what the join keeps, not the
    /// size of the annotation. `keep_residues = false` reads protein ids and lengths only (enough for
    /// `project`). Annotation diagnostics (frame mismatches, dropped transcripts, id collisions, ...) come
    /// back in `Loaded::annotation_diagnostics`.
    fn load_projection(
        &self,
        verbose: u8,
        keep_residues: bool,
        opts: &JoinOptions,
    ) -> Result<Loaded> {
        // Validate the options first so a missing --gff is reported before any file is read.
        let gff_options = self.gff_options()?;
        let proteins = timed(verbose, "proteins_read", || {
            if keep_residues {
                ProteinSet::read(&self.proteins)
            } else {
                ProteinSet::read_lengths(&self.proteins)
            }
        })
        .context("reading protein sequences")?;
        if verbose > 0 {
            eprintln!(
                "sequences: {} proteins, {} residues",
                proteins.len(),
                proteins.total_residues
            );
            let bytes: u64 = self
                .gff
                .iter()
                .filter_map(|p| std::fs::metadata(p).ok())
                .map(|m| m.len())
                .sum();
            eprintln!("[size] gff_bytes_on_disk\t{bytes}");
        }
        for (level, msg) in proteins.report() {
            eprintln!("{level}: {msg}");
        }
        let mut src = GffSource::new(self.gff.clone(), gff_options);
        let mut sink = JoiningSink::new(&proteins, opts);
        timed(verbose, "gff_parse_and_join", || src.stream_into(&mut sink))
            .context("reading annotation")?;
        let (proj, annotation_diagnostics, n_transcripts) =
            timed(verbose, "join_finish", || sink.finish())
                .context("joining annotation to proteins")?;
        if verbose > 0 {
            eprintln!("annotation: {n_transcripts} transcripts");
        }
        Ok(Loaded {
            proj,
            proteins,
            annotation_diagnostics,
            n_transcripts,
        })
    }
}

fn parse_slack(s: &str) -> Result<u64> {
    match s.to_ascii_lowercase().as_str() {
        "none" | "off" | "no" | "disable" => Ok(u64::MAX),
        other => other
            .parse::<u64>()
            .with_context(|| format!("--length-slack must be a number or 'none', got {s:?}")),
    }
}

fn project_options(ignore_phase: bool, clamp_tail: bool, flank: u32) -> ProjectOptions {
    ProjectOptions {
        phase_mode: if ignore_phase {
            PhaseMode::Ignore
        } else {
            PhaseMode::Strict
        },
        tail: if clamp_tail {
            TailPolicy::Clamp
        } else {
            TailPolicy::Drop
        },
        skip_incoherent: false,
        flank,
        taxon_names: Default::default(),
    }
}

/// Print the quality tally and keep the individual rows in a table.
fn report_diagnostics(d: &Diagnostics, out_dir: &Path, verbose: u8) -> Result<()> {
    // Informational records alone (e.g. the transcripts-retained summary) are not worth a file.
    if d.warnings == 0 && d.errors == 0 {
        // A previous run's report must not survive next to fresh results.
        for stale in ["diagnostics.tsv", "diagnostics_detail.tsv"] {
            let _ = std::fs::remove_file(out_dir.join(stale));
        }
        return Ok(());
    }
    let tally = tables::diagnostics_tally_table(d);
    tally.write(out_dir.join("diagnostics.tsv"))?;
    tables::diagnostics_table(d).write(out_dir.join("diagnostics_detail.tsv"))?;
    eprintln!(
        "annotation quality: {} kinds flagged ({} errors, {} warnings) — details in diagnostics.tsv",
        tally.len(),
        d.errors,
        d.warnings
    );
    if verbose > 0 {
        for x in d.retained().iter().take(25) {
            eprintln!("  {} {} {}: {}", x.level, x.kind, x.subject, x.message);
        }
    } else {
        for row in tally.rows.iter().take(5) {
            eprintln!("  {} x{}", row[0], row[1]);
        }
        eprintln!("  (rerun with -v for individual records)");
    }
    Ok(())
}

/// Run `f`, and with `-v` print `[time] <label>\t<seconds>` to stderr (one machine-readable line per stage).
fn timed<T>(verbose: u8, label: &str, f: impl FnOnce() -> T) -> T {
    let t0 = Instant::now();
    let out = f();
    if verbose > 0 {
        eprintln!("[time] {label}\t{:.3}", t0.elapsed().as_secs_f64());
        if let Some((rss, hwm)) = process_memory_kb() {
            // Resident set now, and its high-water mark so far, right after the stage.
            eprintln!("[mem] {label}\t{rss}\t{hwm}");
        }
    }
    out
}

/// `(VmRSS, VmHWM)` in kB from `/proc/self/status`; `None` where that file does not exist (macOS).
fn process_memory_kb() -> Option<(u64, u64)> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let field = |name: &str| -> Option<u64> {
        status
            .lines()
            .find_map(|l| l.strip_prefix(name))?
            .trim()
            .trim_end_matches("kB")
            .trim()
            .parse()
            .ok()
    };
    Some((field("VmRSS:")?, field("VmHWM:")?))
}

fn cmd_project(a: &ProjectArgs, verbose: u8) -> Result<()> {
    let opts = JoinOptions {
        length_slack: parse_slack(&a.length_slack)?,
        ..Default::default()
    };
    let Loaded {
        proj,
        proteins,
        annotation_diagnostics,
        n_transcripts,
    } = a.ann.load_projection(verbose, false, &opts)?;
    std::fs::create_dir_all(&a.out_dir)
        .with_context(|| format!("could not create output directory {}", a.out_dir.display()))?;

    let (n, nu) = timed(verbose, "write_introns_and_unmatched", || -> Result<_> {
        let n = tables::write_introns_tsv(&proj, a.out_dir.join("introns.tsv"))?;
        let nu = tables::unmatched_table(&proj).write(a.out_dir.join("unmatched.tsv"))?;
        Ok((n, nu))
    })?;
    let mut d = annotation_diagnostics;
    d.merge(&proj.diagnostics);
    d.merge(&proteins_diagnostics(&proteins));

    tables::summary_table(&[
        ("stage", "project".to_string()),
        ("transcripts_read", n_transcripts.to_string()),
        ("proteins_in_fasta", proteins.len().to_string()),
        ("proteins_with_introns", proj.joined_proteins.to_string()),
        ("intron_records", proj.rows.len().to_string()),
        ("introns_total", proj.total_introns.to_string()),
        ("introns_in_utr", proj.utr_introns.to_string()),
        ("unmatched_records", nu.to_string()),
        ("out_dir", a.out_dir.display().to_string()),
    ])
    .write(a.out_dir.join("summary_project.tsv"))?;
    report_diagnostics(&d, &a.out_dir, verbose)?;

    println!(
        "wrote {} intron positions for {} protein(s) to {}/introns.tsv",
        n,
        proj.joined_proteins,
        a.out_dir.display()
    );
    Ok(())
}

fn proteins_diagnostics(proteins: &ProteinSet) -> Diagnostics {
    let mut d = Diagnostics::new(64);
    for (level, msg) in proteins.report() {
        d.push(level, "sequence_set", "proteins", msg);
    }
    d
}

/// Project every alignment, writing one directory per alignment.
/// Project every alignment, write its matrix tables, and hand it to `per_alignment` — then drop it.
///
/// Alignments are processed one at a time and not kept: a run over thousands of families would
/// otherwise hold every site matrix (and its per-column residue sets) in memory until the end.
fn build_matrices(
    a: &MatrixArgs,
    verbose: u8,
    mut per_alignment: impl FnMut(&str, &SiteMatrix, &AlignmentProjection) -> Result<()>,
) -> Result<usize> {
    if a.msa.is_empty() {
        bail!("--msa is required: give at least one protein alignment of orthologs");
    }
    let opts = JoinOptions {
        length_slack: parse_slack(&a.length_slack)?,
        ..Default::default()
    };
    let Loaded {
        proj,
        proteins,
        annotation_diagnostics,
        ..
    } = a.ann.load_projection(verbose, true, &opts)?;

    std::fs::create_dir_all(&a.out_dir)
        .with_context(|| format!("could not create output directory {}", a.out_dir.display()))?;
    timed(verbose, "write_introns_and_unmatched", || -> Result<()> {
        tables::write_introns_tsv(&proj, a.out_dir.join("introns.tsv"))?;
        tables::unmatched_table(&proj).write(a.out_dir.join("unmatched.tsv"))?;
        Ok(())
    })?;

    let mut popts = project_options(a.ignore_phase, a.clamp_tail, a.flank);
    if let Some(path) = &a.taxon_map {
        popts.taxon_names = GffSource::load_id_map(path)
            .with_context(|| format!("could not read --taxon-map {}", path.display()))?;
    }
    let mut n_done = 0usize;
    let mut diag = annotation_diagnostics;
    diag.merge(&proj.diagnostics);
    diag.merge(&proteins_diagnostics(&proteins));
    let mut used_stems = std::collections::HashSet::new();
    // Per-stage totals over all alignments, printed with -v.
    let (mut t_read, mut t_project, mut t_build, mut t_tables, mut t_callback) = (
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
    );

    // The diagnostics report is written even when an alignment fails part-way through.
    let outcome: Result<()> = (|| {
        for path in &a.msa {
            let t0 = Instant::now();
            let aln = Alignment::read(path)
                .with_context(|| format!("could not read alignment {}", path.display()))?;
            t_read += t0.elapsed();
            if aln.n_rows() < 2 {
                bail!(
                "alignment {} has {} row(s); at least two sequences are needed to compare intron positions",
                path.display(),
                aln.n_rows()
            );
            }
            let t0 = Instant::now();
            let ap = project_into_alignment(&aln, &proj, &proteins, &popts);
            t_project += t0.elapsed();
            {
                let mut seen: std::collections::HashMap<&str, &str> =
                    std::collections::HashMap::new();
                for t in &ap.taxa {
                    if let Some(first) = seen.insert(&t.taxon_name, &t.row_id) {
                        bail!(
                        "alignment {}: rows {first:?} and {:?} both map to taxon {:?}; a taxon may have one \
                         row per alignment (check --taxon-map, or drop paralogs from the alignment)",
                        aln.name,
                        t.row_id,
                        t.taxon_name
                    );
                    }
                }
            }
            let t0 = Instant::now();
            let m = SiteMatrix::build(&ap);
            t_build += t0.elapsed();

            let t0 = Instant::now();
            let stem = claim_stem(&mut used_stems, &aln.name)?;
            let dir = a.out_dir.join(&stem);
            if !a.combined_only {
                std::fs::create_dir_all(&dir)?;
                tables::sites_table(&ap, &proteins).write(dir.join("sites.tsv"))?;
                tables::taxa_table(&ap, &proteins, &m).write(dir.join("taxa.tsv"))?;
                tables::write_matrix_tsv(&m, dir.join("matrix.tsv"))?;
                tables::presence_table(&m).write(dir.join("presence.tsv"))?;
                let pairs = m.pairwise(a.tolerance, !a.ignore_phase);
                tables::pairwise_table(&pairs).write(dir.join("pairwise.tsv"))?;
                tables::sharing_matrix_table_from(&m, &pairs)
                    .write(dir.join("sharing_matrix.tsv"))?;
                tables::spectrum_table(&m).write(dir.join("spectrum.tsv"))?;
                tables::counts_table(&score::score_counts(&m))
                    .write(dir.join("counts_sites.tsv"))?;
                std::fs::write(
                    dir.join("characters.nex"),
                    intronfate_core::matrix::characters_nexus(&m),
                )
                .with_context(|| {
                    format!("could not write {}", dir.join("characters.nex").display())
                })?;
            }

            t_tables += t0.elapsed();
            for taxon in &ap.taxa {
                if !taxon.has_annotation {
                    diag.warn(
                    taxon.unannotated_reason.unwrap_or("unannotated_taxon"),
                    &taxon.row_id,
                    "alignment row has no annotation, so its absences are scored as missing data (?), not losses",
                );
                }
                if taxon.report.tail_sites > 0 {
                    diag.warn(
                        "tail_sites_dropped",
                        &taxon.row_id,
                        format!(
                        "{} intron position(s) anchor outside the {} residues this row covers; {}",
                        taxon.report.tail_sites,
                        taxon.residues_aligned,
                        if a.clamp_tail {
                            format!(
                                "{} placed on the last column (--clamp-tail), the rest dropped",
                                taxon.report.clamped
                            )
                        } else {
                            "dropped".to_string()
                        }
                    ),
                    );
                }
            }
            if m.masked_present > 0 {
                diag.warn(
                "present_not_observable",
                &stem,
                format!(
                    "{} annotated intron(s) sit where the taxon is not aligned around the site (flank {}); they are reported as `?`, not as presence",
                    m.masked_present, m.flank
                ),
            );
            }
            if ap.coverage.n_absent_from_alignment > 0 {
                // One record per alignment, not one per protein. In a multi-family run this is true of
                // nearly every protein in nearly every alignment, so it is informational.
                diag.info(
                    "proteins_not_in_alignment",
                    &stem,
                    format!(
                        "{} annotated protein(s) have no row in this alignment (first: {})",
                        ap.coverage.n_absent_from_alignment,
                        ap.coverage.absent_examples.join(", ")
                    ),
                );
            }
            if verbose > 0 {
                eprintln!("{stem}: {}", m.describe());
            }
            let t0 = Instant::now();
            per_alignment(&stem, &m, &ap)?;
            t_callback += t0.elapsed();
            n_done += 1;
        }
        Ok(())
    })();
    if verbose > 0 {
        for (label, t) in [
            ("alignments_read", t_read),
            ("alignments_project", t_project),
            ("matrix_build", t_build),
            ("matrix_tables_written", t_tables),
            ("per_alignment_callback(score)", t_callback),
        ] {
            eprintln!("[time] {label}\t{:.3}", t.as_secs_f64());
        }
        eprintln!("[size] alignments\t{n_done}");
    }
    report_diagnostics(&diag, &a.out_dir, verbose)?;
    outcome?;
    Ok(n_done)
}

fn cmd_matrix(a: &MatrixArgs, verbose: u8) -> Result<()> {
    let mut parts: Vec<SiteMatrix> = Vec::new();
    build_matrices(a, verbose, |stem, m, _| {
        if !a.combined_only {
            println!(
                "{stem}: {} sites x {} taxa -> {}/{stem}/matrix.tsv",
                m.n_sites(),
                m.n_taxa(),
                a.out_dir.display()
            );
        }
        if a.combine {
            parts.push(m.clone());
        }
        Ok(())
    })?;
    if a.combine {
        let c = write_combined(&parts, &a.out_dir.join("combined"), a.tolerance)?;
        println!(
            "combined: {} sites from {} alignments x {} taxa -> {}/combined/matrix.tsv",
            c.n_sites(),
            parts.len(),
            c.n_taxa(),
            a.out_dir.display()
        );
    }
    Ok(())
}

/// Stack `parts` into one matrix and write it with its matrix-level tables (not the per-alignment projection tables) into `dir`.
fn write_combined(parts: &[SiteMatrix], dir: &Path, tolerance: u32) -> Result<SiteMatrix> {
    let c = SiteMatrix::combine(parts, "combined")?;
    if tolerance >= c.family_pad {
        bail!(
            "--tolerance {tolerance} is not below the {} columns that separate alignments in a combined matrix",
            c.family_pad
        );
    }
    std::fs::create_dir_all(dir)?;
    tables::write_matrix_tsv(&c, dir.join("matrix.tsv"))?;
    tables::presence_table(&c).write(dir.join("presence.tsv"))?;
    let strict = matches!(c.phase_mode, intronfate_core::project::PhaseMode::Strict);
    let pairs = c.pairwise(tolerance, strict);
    tables::pairwise_table(&pairs).write(dir.join("pairwise.tsv"))?;
    tables::sharing_matrix_table_from(&c, &pairs).write(dir.join("sharing_matrix.tsv"))?;
    tables::spectrum_table(&c).write(dir.join("spectrum.tsv"))?;
    tables::counts_table(&score::score_counts(&c)).write(dir.join("counts_sites.tsv"))?;
    Ok(c)
}

fn cmd_convert(a: &ConvertArgs) -> Result<()> {
    if let Some(table) = &a.from_malin {
        let text = std::fs::read_to_string(table)
            .with_context(|| format!("could not read MALIN table {}", table.display()))?;
        let stem = table
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("malin")
            .to_string();
        let matrices = intronfate_core::malin::from_malin(&text, &stem)
            .with_context(|| format!("could not parse {}", table.display()))?;
        let mut used = std::collections::HashSet::new();
        for m in &matrices {
            let dir = a.out.join(claim_stem(&mut used, &m.alignment)?);
            tables::matrix_table(m).write(dir.join("matrix.tsv"))?;
        }
        println!(
            "{} taxa, {} alignment(s), {} sites -> {}",
            matrices[0].n_taxa(),
            matrices.len(),
            matrices.iter().map(|m| m.n_sites()).sum::<usize>(),
            a.out.display()
        );
        return Ok(());
    }
    if a.to_malin.is_empty() {
        bail!("give --from-malin <table>, or --to-malin <matrix.tsv>... (or --matrix-list / --matrix-dir) to export");
    }
    let matrices = a
        .to_malin
        .iter()
        .map(|p| {
            SiteMatrix::read_tsv(p)
                .with_context(|| format!("could not read matrix {}", p.display()))
        })
        .collect::<Result<Vec<_>>>()?;
    let refs: Vec<&SiteMatrix> = matrices.iter().collect();
    let text = intronfate_core::malin::to_malin(&refs)?;
    intronfate_core::util::write_atomic(&a.out, text.as_bytes())?;
    println!(
        "{} alignment(s), {} sites -> {}",
        matrices.len(),
        matrices.iter().map(|m| m.n_sites()).sum::<usize>(),
        a.out.display()
    );
    Ok(())
}

fn cmd_analyze(a: &AnalyzeArgs, verbose: u8) -> Result<()> {
    let methods = resolve_methods(&a.methods)?;
    let mut parts: Vec<SiteMatrix> = Vec::new();
    build_matrices(&a.matrix, verbose, |stem, m, _| {
        if a.matrix.combine {
            parts.push(m.clone());
        }
        if a.matrix.combined_only {
            return Ok(());
        }
        let dir = a.matrix.out_dir.join(stem);
        let summary = score_matrix(
            m,
            &dir,
            &ScoreOpts {
                methods: methods.clone(),
                tree: a.tree.clone(),
                outgroups: a.outgroups.clone(),
                min_taxa: a.min_taxa,
                min_observed: a.min_observed,
                reconstructions: a.reconstructions,
                matrix_label: dir.join("matrix.tsv").display().to_string(),
            },
            verbose,
        )
        .with_context(|| format!("scoring alignment {stem}"))?;
        if let Some(line) = summary {
            println!("{stem}: {line}");
        }
        Ok(())
    })?;
    if a.matrix.combine {
        let dir = a.matrix.out_dir.join("combined");
        let c = write_combined(&parts, &dir, a.matrix.tolerance)?;
        let summary = score_matrix(
            &c,
            &dir,
            &ScoreOpts {
                methods,
                tree: a.tree.clone(),
                outgroups: a.outgroups.clone(),
                min_taxa: a.min_taxa,
                min_observed: a.min_observed,
                reconstructions: a.reconstructions,
                matrix_label: dir.join("matrix.tsv").display().to_string(),
            },
            verbose,
        )
        .context("scoring the combined matrix")?;
        if let Some(line) = summary {
            println!("combined: {line}");
        }
        println!(
            "combined: {} sites from {} alignments x {} taxa -> {}",
            c.n_sites(),
            parts.len(),
            c.n_taxa(),
            dir.display()
        );
    }
    Ok(())
}

fn cmd_score(a: &ScoreArgs, verbose: u8) -> Result<()> {
    if a.matrix.is_empty() {
        bail!("--matrix is required: give a matrix.tsv, or use --matrix-list / --matrix-dir");
    }
    if a.combine {
        let parts = a
            .matrix
            .iter()
            .map(|p| {
                SiteMatrix::read_tsv(p)
                    .with_context(|| format!("could not read matrix {}", p.display()))
            })
            .collect::<Result<Vec<_>>>()?;
        let c = write_combined(&parts, &a.out_dir, a.tolerance)?;
        let summary = score_matrix(
            &c,
            &a.out_dir,
            &ScoreOpts {
                methods: resolve_methods(&a.methods)?,
                tree: a.tree.clone(),
                outgroups: a.outgroups.clone(),
                min_taxa: a.min_taxa,
                min_observed: a.min_observed,
                reconstructions: a.reconstructions,
                matrix_label: format!("{} matrices combined", parts.len()),
            },
            verbose,
        )?;
        if let Some(line) = summary {
            println!("{line}");
        }
        println!(
            "combined {} matrices: {} sites x {} taxa -> {}",
            parts.len(),
            c.n_sites(),
            c.n_taxa(),
            a.out_dir.display()
        );
        return Ok(());
    }
    if !a.many() {
        return score_one(a, &a.matrix[0], &a.out_dir, verbose);
    }
    let mut used = std::collections::HashSet::new();
    for path in &a.matrix {
        let stem = claim_stem(&mut used, &matrix_stem(path))?;
        let out = a.out_dir.join(&stem);
        println!("{stem}:");
        score_one(a, path, &out, verbose).with_context(|| format!("scoring {}", path.display()))?;
    }
    println!(
        "scored {} matrices into {}",
        a.matrix.len(),
        a.out_dir.display()
    );
    Ok(())
}

fn score_one(a: &ScoreArgs, matrix: &Path, out_dir: &Path, verbose: u8) -> Result<()> {
    let m = SiteMatrix::read_tsv(matrix)
        .with_context(|| format!("could not read matrix {}", matrix.display()))?;
    std::fs::create_dir_all(out_dir)?;
    // Re-emit pairwise statistics at the requested tolerance without re-reading the annotation.
    let pairs = m.pairwise(a.tolerance, true);
    tables::pairwise_table(&pairs).write(out_dir.join("pairwise.tsv"))?;
    tables::sharing_matrix_table_from(&m, &pairs).write(out_dir.join("sharing_matrix.tsv"))?;

    let methods = resolve_methods(&a.methods)?;
    let summary = score_matrix(
        &m,
        out_dir,
        &ScoreOpts {
            methods,
            tree: a.tree.clone(),
            outgroups: a.outgroups.clone(),
            min_taxa: a.min_taxa,
            min_observed: a.min_observed,
            reconstructions: a.reconstructions,
            matrix_label: matrix.display().to_string(),
        },
        verbose,
    )?;
    if let Some(line) = summary {
        println!("{}", line);
    }
    println!(
        "scored {} sites x {} taxa from {}",
        m.n_sites(),
        m.n_taxa(),
        matrix.display()
    );
    Ok(())
}

struct ScoreOpts {
    methods: Vec<Method>,
    tree: Option<PathBuf>,
    outgroups: Vec<String>,
    min_taxa: usize,
    min_observed: usize,
    reconstructions: bool,
    matrix_label: String,
}

/// Score one matrix; returns a one-line summary of what was computed.
fn score_matrix(m: &SiteMatrix, dir: &Path, o: &ScoreOpts, verbose: u8) -> Result<Option<String>> {
    let mut summary: Vec<(&str, String)> = vec![
        ("stage", "score".into()),
        ("matrix", o.matrix_label.clone()),
        ("alignment", m.alignment.clone()),
        ("sites", m.n_sites().to_string()),
        ("taxa", m.n_taxa().to_string()),
        ("annotated_taxa", m.n_annotated.to_string()),
        (
            "methods",
            o.methods
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(","),
        ),
    ];
    let mut notes = Vec::new();

    if o.methods.contains(&Method::Counts) {
        let c = score::score_counts(m);
        let n = tables::counts_table(&c).write(dir.join("counts_sites.tsv"))?;
        summary.push(("counts_sites", n.to_string()));
        if verbose > 0 {
            eprintln!("counts: {n} positions");
        }
    }

    let needs_tree = o
        .methods
        .iter()
        .any(|x| matches!(x, Method::Dollo | Method::MaxParsimony));
    if !needs_tree {
        tables::summary_table(&summary).write(dir.join("summary_score.tsv"))?;
        return Ok(Some(format!("counts over {} sites", m.n_sites())));
    }
    let Some(tree_path) = &o.tree else {
        bail!("--tree is required to score with dollo or mp (or select only --methods counts)");
    };
    let tree = Tree::read(tree_path)
        .with_context(|| format!("could not read tree {}", tree_path.display()))?;
    if verbose > 0 {
        eprintln!(
            "tree: {} tips, {}",
            tree.n_tips(),
            if tree.looks_unrooted() {
                "unrooted"
            } else {
                "rooted"
            }
        );
    }
    let (data, report) = score::prepare_with(
        m,
        &tree,
        &o.outgroups,
        &score::PrepareOptions {
            min_taxa: o.min_taxa,
            min_observed: o.min_observed,
        },
    )
    .context("matching matrix taxa to tree tips")?;
    if !report.taxa_not_in_tree.is_empty() {
        eprintln!(
            "warn: {} matrix taxon(s) are not tips of the tree and were not scored: {}",
            report.taxa_not_in_tree.len(),
            intronfate_core::util::ellipsis(&report.taxa_not_in_tree.join(", "), 200)
        );
    }
    if !report.tree_tip_without_data.is_empty() {
        eprintln!(
            "warn: {} tree tip(s) have no matrix row and are scored as missing (?): {}",
            report.tree_tip_without_data.len(),
            intronfate_core::util::ellipsis(&report.tree_tip_without_data.join(", "), 200)
        );
    }
    if data.tree.looks_unrooted() && o.outgroups.is_empty() {
        eprintln!(
            "warn: the tree looks unrooted ({} children at the root) and no --outgroup was given, so \
             gain/loss placement follows the Newick as written",
            data.tree.children(data.tree.root()).len()
        );
        notes.push("unrooted".to_string());
    }

    summary.extend(prepare_summary(&data, &report, &o.outgroups));
    let mut recon_rows: Vec<(&str, &str, Option<Vec<u8>>)> = Vec::new();

    let dollo_out = o
        .methods
        .contains(&Method::Dollo)
        .then(|| score::score_dollo(&data, o.reconstructions));
    if let Some(d) = &dollo_out {
        let n = tables::dollo_table(d).write(dir.join("dollo_sites.tsv"))?;
        tables::dollo_branch_table(d).write(dir.join("dollo_branches.tsv"))?;
        std::fs::write(
            dir.join("tree_gains.nwk"),
            tables::tree_with_gain_counts(&data, d),
        )?;
        summary.push(("dollo_sites", n.to_string()));
        summary.push(("dollo_gains", d.total_gains.to_string()));
        summary.push(("dollo_losses", d.total_losses.to_string()));
        summary.push(("dollo_compatible_sites", d.compatible.to_string()));
        summary.push(("dollo_needing_losses", d.incompatible.to_string()));
        summary.push(("dollo_monophyletic_sites", d.monophyletic.to_string()));
        if d.incompatible > 0 {
            eprintln!(
                "note: {} position(s) cannot be explained by a single origin without {} independent \
                 loss event(s) — candidates for convergent insertion, or for the topology",
                d.incompatible, d.total_losses
            );
        }
        for s in &d.sites {
            recon_rows.push(("dollo", s.label.as_str(), s.reconstruction.clone()));
        }
        if verbose > 0 {
            eprintln!(
                "dollo: {} gains, {} losses over {} sites",
                d.total_gains, d.total_losses, n
            );
        }
        notes.push(format!("dollo {}g/{}l", d.total_gains, d.total_losses));
    }

    let mp_out = o
        .methods
        .contains(&Method::MaxParsimony)
        .then(|| score::score_mp(&data, o.reconstructions));
    if let Some(mp) = &mp_out {
        let n = tables::mp_table(mp).write(dir.join("mp_sites.tsv"))?;
        summary.push(("mp_sites", n.to_string()));
        summary.push(("mp_steps", mp.total_steps.to_string()));
        summary.push(("mp_informative_sites", mp.informative.to_string()));
        summary.push(("mp_ambiguous_sites", mp.ambiguous.to_string()));
        summary.push(("mp_consistency_index", format!("{:.6}", mp.ci)));
        for s in &mp.sites {
            recon_rows.push(("mp", s.label.as_str(), s.reconstruction.clone()));
        }
        if verbose > 0 {
            eprintln!(
                "mp: {} steps over {} informative sites",
                mp.total_steps, mp.informative
            );
        }
        notes.push(format!("mp {} steps", mp.total_steps));
    }

    if o.reconstructions && !recon_rows.is_empty() {
        let t = tables::reconstruction_table(&data, &recon_rows);
        summary.push(("reconstructed_nodes", t.len().to_string()));
        t.write(dir.join("reconstructions.tsv"))?;
    }
    tables::summary_table(&summary).write(dir.join("summary_score.tsv"))?;
    Ok(Some(notes.join(", ")))
}

fn prepare_summary(
    data: &ScoringData,
    report: &PrepareReport,
    outgroups: &[String],
) -> Vec<(&'static str, String)> {
    vec![
        ("tree_tips", data.tree.n_tips().to_string()),
        ("tips_scored", report.used.len().to_string()),
        (
            "tips_without_data",
            report.tree_tip_without_data.len().to_string(),
        ),
        (
            "matrix_taxa_not_in_tree",
            report.taxa_not_in_tree.len().to_string(),
        ),
        ("characters_scored", data.characters.len().to_string()),
        ("characters_skipped", report.uninformative.to_string()),
        ("characters_too_sparse", report.too_sparse.to_string()),
        (
            "rooting",
            if !outgroups.is_empty() {
                format!("outgroup:{}", outgroups.join("+"))
            } else if data.tree.looks_unrooted() {
                "as-given (unrooted)".to_string()
            } else {
                "as-given".to_string()
            },
        ),
    ]
}

fn resolve_methods(requested: &[String]) -> Result<Vec<Method>> {
    if requested.is_empty() {
        return Ok(vec![Method::Counts, Method::Dollo, Method::MaxParsimony]);
    }
    let mut out = Vec::new();
    for name in requested {
        let m = Method::parse(name).ok_or_else(|| {
            anyhow::anyhow!("unknown method {name:?}: choose dollo, mp or counts")
        })?;
        if !out.contains(&m) {
            out.push(m);
        }
    }
    Ok(out)
}

fn safe_name(s: &str) -> String {
    let out: String = s
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '.' || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    // `.` and `..` are legal characters above but would name the current or parent directory.
    if out.is_empty() || out.chars().all(|c| c == '.') {
        "alignment".to_string()
    } else {
        out
    }
}

/// Reserve an output folder name for an alignment; names are compared case-insensitively because
/// macOS and Windows filesystems would otherwise write two alignments into one folder.
fn claim_stem(used: &mut std::collections::HashSet<String>, name: &str) -> Result<String> {
    let stem = safe_name(name);
    if !used.insert(stem.to_lowercase()) {
        bail!(
            "two alignments map to the output folder {stem:?} (e.g. {name:?}); rename one so results are \
             not silently overwritten"
        );
    }
    Ok(stem)
}

// ---------------------------------------------------------------------------------------------
// Demo
// ---------------------------------------------------------------------------------------------

/// Built-in example: five orthologs whose intron positions are shared, lost, and newly gained.
///
/// Coding layout (all plus strand, all introns at codon boundaries so the story is about position
/// rather than phase):
///
/// | taxon      | CDS segments (nt)      | protein-space introns | residues |
/// |------------|------------------------|-----------------------|----------|
/// | species1   | 12 + 15 + 12 = 39      | cds 12, 27            | 13       |
/// | species2/3 | 12 + 6 + 24 = 42       | cds 12, 18            | 14       |
/// | species4   | 27 + 12 = 39           | cds 27                | 13       |
/// | outgroup   | 12 + 15 + 12 = 39      | cds 12, 27            | 13       |
///
/// Projected onto the alignment, cds 12 lands on column 4 for everyone, so: column 4 is ancestral
/// (lost in species4), column 9 is ancestral (lost in species2 and 3 together), and column 6 is a
/// derived insertion on the species2/species3 branch.
fn cmd_demo(a: &DemoArgs, verbose: u8) -> Result<()> {
    let dir = &a.out_dir;
    std::fs::create_dir_all(dir)
        .with_context(|| format!("could not create demo output directory {}", dir.display()))?;
    let (gff, faa, msa, nwk) = (demo_gff(), demo_faa(), demo_msa(), demo_tree());
    std::fs::write(dir.join("demo.gff3"), &gff)?;
    std::fs::write(dir.join("demo.faa"), &faa)?;
    std::fs::write(dir.join("demo.aln.fasta"), &msa)?;
    std::fs::write(dir.join("demo.nwk"), &nwk)?;

    let mut set = AnnotationSet::default();
    GffSource::new(vec![dir.join("demo.gff3")], GffOptions::default())
        .stream_into(&mut set)
        .context("demo annotation")?;
    let proteins = ProteinSet::read([dir.join("demo.faa")])?;
    let proj = join_annotations(&set.transcripts, &proteins, &JoinOptions::default())?;
    let aln = Alignment::from_fasta_str(&msa, "demo")?;
    let ap = project_into_alignment(&aln, &proj, &proteins, &ProjectOptions::default());
    let m = SiteMatrix::build(&ap);

    let tree = Tree::from_newick(&nwk)?;
    let (data, report) = score::prepare(&m, &tree, &["outgroup".into()], 1)?;
    let d = score::score_dollo(&data, true);
    let mp = score::score_mp(&data, true);

    tables::introns_table(&proj).write(dir.join("introns.tsv"))?;
    tables::sites_table(&ap, &proteins).write(dir.join("sites.tsv"))?;
    tables::matrix_table(&m).write(dir.join("matrix.tsv"))?;
    tables::spectrum_table(&m).write(dir.join("spectrum.tsv"))?;
    tables::pairwise_table(&m.pairwise(0, true)).write(dir.join("pairwise.tsv"))?;
    tables::dollo_table(&d).write(dir.join("dollo_sites.tsv"))?;
    tables::mp_table(&mp).write(dir.join("mp_sites.tsv"))?;
    tables::reconstruction_table(
        &data,
        &d.sites
            .iter()
            .map(|s| ("dollo", s.label.as_str(), s.reconstruction.clone()))
            .collect::<Vec<_>>(),
    )
    .write(dir.join("reconstructions.tsv"))?;
    std::fs::write(
        dir.join("tree_gains.nwk"),
        tables::tree_with_gain_counts(&data, &d),
    )?;
    tables::summary_table(&[
        ("stage", "demo".into()),
        ("sites", m.n_sites().to_string()),
        ("taxa", m.n_taxa().to_string()),
        ("dollo_gains", d.total_gains.to_string()),
        ("dollo_losses", d.total_losses.to_string()),
        ("mp_steps", mp.total_steps.to_string()),
        (
            "tips_without_data",
            report.tree_tip_without_data.len().to_string(),
        ),
    ])
    .write(dir.join("summary.tsv"))?;

    print!("{}", tables::matrix_table(&m).to_tsv());
    println!();
    print!("{}", tables::dollo_table(&d).to_tsv());
    println!();
    print!("{}", tables::mp_table(&mp).to_tsv());
    if verbose > 0 {
        eprintln!("rooted on the outgroup: {}", data.tree.to_newick());
    }
    println!("\ndemo inputs and tables written to {}", dir.display());
    println!(
        "expect 3 sites, 3 Dollo gains, 2 Dollo losses, 3 mp steps: got {} sites, {} gains, {} losses, {} steps",
        m.n_sites(),
        d.total_gains,
        d.total_losses,
        mp.total_steps
    );
    Ok(())
}

/// CDS layout for the demo; `frame` is always 0 because every segment is a multiple of 3 nt.
fn demo_gff() -> String {
    let mut s = String::from("##gff-version 3\n#!demo — intron projection example\n");
    let layout: [(&str, &[(u64, u64)]); 5] = [
        ("species1", &[(1, 12), (100, 114), (200, 211)]),
        ("species2", &[(1, 12), (100, 105), (200, 223)]),
        ("species3", &[(1, 12), (100, 105), (200, 223)]),
        ("species4", &[(1, 27), (100, 111)]),
        ("outgroup", &[(1, 12), (100, 114), (200, 211)]),
    ];
    for (taxon, segs) in layout {
        for (i, (start, end)) in segs.iter().enumerate() {
            s.push_str(&format!(
                "chr1\tdemo\tCDS\t{start}\t{end}\t.\t+\t0\tID={taxon}.c{};Parent={taxon}.t1;protein_id={taxon};gene={taxon}.g\n",
                i + 1
            ));
        }
    }
    s
}

fn demo_faa() -> String {
    let mut s = String::new();
    for (taxon, seq) in [
        ("species1", "MKVLGAWETSIGL"),
        ("species2", "MKVLGAWETSIGLQ"),
        ("species3", "MKVLGAWETSIGLQ"),
        ("species4", "MKVLGAWETSIGL"),
        ("outgroup", "MKVLGAWETSIGL"),
    ] {
        s.push_str(&format!(">{taxon} demo ortholog\n{seq}\n"));
    }
    s
}

fn demo_msa() -> String {
    let mut s = String::new();
    for (taxon, seq) in [
        ("species1", "MKVLGAWETSIGL--"),
        ("species2", "MKVLGAWETSIGLQ-"),
        ("species3", "MKVLGAWETSIGLQ-"),
        ("species4", "MKVLGAWETSIGL--"),
        ("outgroup", "MKVLGAWETSIGL--"),
    ] {
        s.push_str(&format!(">{taxon}\n{seq}\n"));
    }
    s
}

fn demo_tree() -> String {
    "((species1,(species2,species3))A,species4,outgroup);".to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_accepts_the_documented_invocations() {
        Cli::try_parse_from(["ip", "project", "--gff", "a.gff3", "-p", "a.faa"]).unwrap();
        Cli::try_parse_from([
            "ip", "matrix", "--gff", "a.gff3", "-p", "a.faa", "--msa", "x.faa",
        ])
        .unwrap();
        Cli::try_parse_from([
            "ip", "score", "--matrix", "m.tsv", "-t", "t.nwk", "-m", "dollo,mp",
        ])
        .unwrap();
        Cli::try_parse_from([
            "ip", "analyze", "--gff", "a.gff3", "-p", "a.faa", "--msa", "x.faa", "-t", "t.nwk",
        ])
        .unwrap();
        // A missing --msa is a run-time error with a clear message, not a panic.
        let c = Cli::try_parse_from(["ip", "matrix", "--gff", "a.gff3", "-p", "a.faa"]).unwrap();
        match c.command {
            Command::Matrix(a) => {
                let e = build_matrices(&a, 0, |_, _, _| Ok(())).unwrap_err();
                assert!(e.to_string().contains("--msa is required"), "{e}");
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn slack_and_methods_parse() {
        assert_eq!(parse_slack("3").unwrap(), 3);
        assert_eq!(parse_slack("none").unwrap(), u64::MAX);
        assert!(parse_slack("nope").is_err());
        assert_eq!(resolve_methods(&[]).unwrap().len(), 3);
        assert_eq!(
            resolve_methods(&["dollo".into(), "dollo".into()]).unwrap(),
            vec![Method::Dollo]
        );
        assert!(resolve_methods(&["nope".into()]).is_err());
    }

    #[test]
    fn safe_name_replaces_unexpected_characters() {
        assert_eq!(safe_name("fam/one.fasta"), "fam_one.fasta");
        assert_eq!(safe_name(""), "alignment");
        assert_eq!(safe_name("a b"), "a_b");
    }

    #[test]
    fn demo_inputs_are_internally_consistent() {
        // CDS length must imply the protein length, or the length check would fire on our own demo.
        let dir = std::env::temp_dir().join(format!("ip_demochk_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gff = dir.join("demo.gff3");
        std::fs::write(&gff, demo_gff()).unwrap();
        let mut set = AnnotationSet::default();
        GffSource::new(vec![gff], GffOptions::default())
            .stream_into(&mut set)
            .unwrap();
        assert_eq!(
            set.diagnostics.count_of("frame_mismatch"),
            0,
            "{:?}",
            set.diagnostics.retained()
        );
        let proteins = ProteinSet::read([demo_faa_fixture(&dir)]).unwrap();
        let proj = join_annotations(&set.transcripts, &proteins, &JoinOptions::default()).unwrap();
        assert_eq!(
            proj.diagnostics.count_of("length_mismatch"),
            0,
            "{:?}",
            proj.diagnostics.retained()
        );
        assert_eq!(proj.diagnostics.count_of("site_beyond_protein_end"), 0);
        assert_eq!(proj.joined_proteins, 5);

        // And the projection must produce the three intended columns.
        let aln = Alignment::from_fasta_str(&demo_msa(), "demo").unwrap();
        let ap = project_into_alignment(&aln, &proj, &proteins, &ProjectOptions::default());
        let m = SiteMatrix::build(&ap);
        let cols: Vec<(u32, u8, u32)> = m
            .sites
            .iter()
            .map(|s| (s.col, s.phase, s.n_taxa()))
            .collect();
        assert_eq!(
            cols,
            vec![(4, 0, 4), (6, 0, 2), (9, 0, 3)],
            "ancestral 4:0, gained 6:0, ancestral 9:0"
        );
    }

    fn demo_faa_fixture(dir: &std::path::Path) -> PathBuf {
        let p = dir.join("demo.faa");
        std::fs::write(&p, demo_faa()).unwrap();
        p
    }

    #[test]
    fn demo_end_to_end_scores_as_documented() {
        let dir = std::env::temp_dir().join(format!("ip_demo_{}", std::process::id()));
        let a = DemoArgs {
            out_dir: dir.clone(),
        };
        cmd_demo(&a, 0).expect("demo runs");

        let m = SiteMatrix::read_tsv(dir.join("matrix.tsv")).unwrap();
        assert_eq!(m.n_sites(), 3);
        assert_eq!(m.n_taxa(), 5);
        assert_eq!(
            m.alen, 15,
            "alignment length round-trips through the metadata line"
        );
        assert_eq!(
            m.taxa,
            vec!["species1", "species2", "species3", "species4", "outgroup"]
        );

        let tree = Tree::from_newick(&demo_tree()).unwrap();
        let (data, report) = score::prepare(&m, &tree, &["outgroup".into()], 1).unwrap();
        assert!(
            report.taxa_not_in_tree.is_empty(),
            "{:?}",
            report.taxa_not_in_tree
        );
        let d = score::score_dollo(&data, false);
        // One origin per position (3 gains). Losses: column 4 was lost in species4, column 9 was
        // lost once on the branch leading to the (species2,species3) clade, column 6 nowhere.
        assert_eq!(
            d.total_gains,
            3,
            "{:?}",
            d.sites
                .iter()
                .map(|s| (&s.label, s.gains, s.losses))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            d.total_losses,
            2,
            "{:?}",
            d.sites
                .iter()
                .map(|s| (&s.label, s.losses))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            d.compatible, 1,
            "only the species2/3-specific position needs no loss"
        );
        let mp = score::score_mp(&data, false);
        // Unrestricted MP is cheaper on this fixture: one change per position, because in both
        // ancestral cases the losing taxa are sisters (or a single taxon).
        assert_eq!(
            mp.total_steps,
            3,
            "{:?}",
            mp.sites
                .iter()
                .map(|s| (&s.label, s.steps))
                .collect::<Vec<_>>()
        );
        assert!(dir.join("tree_gains.nwk").exists());
        assert!(dir.join("reconstructions.tsv").exists());
    }

    #[test]
    fn matrix_files_survive_a_round_trip_and_re_score() {
        let dir = std::env::temp_dir().join(format!("ip_rt_{}", std::process::id()));
        cmd_demo(
            &DemoArgs {
                out_dir: dir.clone(),
            },
            0,
        )
        .unwrap();
        // Score from the file alone, with no annotation in sight.
        let a = ScoreArgs {
            matrix: vec![dir.join("matrix.tsv")],
            matrix_list: vec![],
            matrix_dir: vec![],
            matrix_ext: vec![],
            combine: false,
            tree: Some(dir.join("demo.nwk")),
            methods: vec!["dollo".into()],
            outgroups: vec!["outgroup".into()],
            min_taxa: 1,
            min_observed: 0,
            tolerance: 0,
            reconstructions: false,
            out_dir: dir.join("scored"),
        };
        cmd_score(&a, 0).expect("scoring a matrix file works standalone");
        assert!(dir.join("scored/dollo_sites.tsv").exists());
        assert!(dir.join("scored/summary_score.tsv").exists());
        let text = std::fs::read_to_string(dir.join("scored/dollo_sites.tsv")).unwrap();
        assert_eq!(text.lines().count(), 4, "header + 3 sites");
    }

    #[test]
    fn score_without_a_tree_still_produces_counts() {
        let dir = std::env::temp_dir().join(format!("ip_nocount_{}", std::process::id()));
        cmd_demo(
            &DemoArgs {
                out_dir: dir.clone(),
            },
            0,
        )
        .unwrap();
        let a = ScoreArgs {
            matrix: vec![dir.join("matrix.tsv")],
            matrix_list: vec![],
            matrix_dir: vec![],
            matrix_ext: vec![],
            combine: false,
            tree: None,
            methods: vec!["counts".into()],
            outgroups: vec![],
            min_taxa: 1,
            min_observed: 0,
            tolerance: 0,
            reconstructions: false,
            out_dir: dir.join("counts"),
        };
        cmd_score(&a, 0).unwrap();
        assert!(dir.join("counts/counts_sites.tsv").exists());
        assert!(!dir.join("counts/dollo_sites.tsv").exists());
    }

    #[test]
    fn tree_is_required_when_a_tree_based_method_is_asked_for() {
        let dir = std::env::temp_dir().join(format!("ip_needtree_{}", std::process::id()));
        cmd_demo(
            &DemoArgs {
                out_dir: dir.clone(),
            },
            0,
        )
        .unwrap();
        let a = ScoreArgs {
            matrix: vec![dir.join("matrix.tsv")],
            matrix_list: vec![],
            matrix_dir: vec![],
            matrix_ext: vec![],
            combine: false,
            tree: None,
            methods: vec!["dollo".into()],
            outgroups: vec![],
            min_taxa: 1,
            min_observed: 0,
            tolerance: 0,
            reconstructions: false,
            out_dir: dir.join("x"),
        };
        let e = cmd_score(&a, 0).unwrap_err();
        assert!(e.to_string().contains("--tree is required"), "{e}");
    }

    #[test]
    fn analyze_runs_the_whole_pipeline() {
        let dir = std::env::temp_dir().join(format!("ip_analyze_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gff = dir.join("a.gff3");
        let faa = dir.join("a.faa");
        let aln = dir.join("a.fasta");
        let nwk = dir.join("a.nwk");
        std::fs::write(&gff, demo_gff()).unwrap();
        std::fs::write(&faa, demo_faa()).unwrap();
        std::fs::write(&aln, demo_msa()).unwrap();
        std::fs::write(&nwk, demo_tree()).unwrap();
        let args = AnalyzeArgs {
            matrix: MatrixArgs {
                ann: AnnotationArgs {
                    gff: vec![gff],
                    gff_list: vec![],
                    gff_dir: vec![],
                    gff_ext: vec![],
                    proteins: vec![faa],
                    proteins_list: vec![],
                    proteins_dir: vec![],
                    proteins_ext: vec![],
                    transcript_attrs: GffOptions::default().transcript_attrs,
                    protein_attrs: GffOptions::default().protein_attrs,
                    gene_attrs: GffOptions::default().gene_attrs,
                    cds_types: vec!["cds".into()],
                    id_map: None,
                    strip_version: false,
                    derive_cds_from_exons: false,
                    keep_noncoding: false,
                    strict: false,
                },
                msa: vec![aln],
                msa_list: vec![],
                msa_dir: vec![],
                msa_ext: vec![],
                combine: false,
                combined_only: false,
                out_dir: dir.join("out"),
                ignore_phase: false,
                clamp_tail: false,
                flank: 1,
                taxon_map: None,
                tolerance: 0,
                length_slack: "3".into(),
            },
            tree: Some(nwk),
            methods: vec![],
            outgroups: vec!["outgroup".into()],
            min_taxa: 1,
            min_observed: 0,
            reconstructions: true,
        };
        cmd_analyze(&args, 1).expect("analyze runs");
        assert!(dir.join("out/a/matrix.tsv").exists());
        assert!(dir.join("out/a/dollo_sites.tsv").exists());
        assert!(dir.join("out/a/reconstructions.tsv").exists());
        assert!(dir.join("out/introns.tsv").exists());
    }

    #[test]
    fn ignore_phase_merges_positions_that_differ_only_in_phase() {
        let dir = std::env::temp_dir().join(format!("ip_phase_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gff = dir.join("p.gff3");
        // Two species with introns one nucleotide apart inside the same codon: 12 (phase 0) and 13 (phase 1).
        std::fs::write(
            &gff,
            "c1\tx\tCDS\t1\t12\t.\t+\t0\tID=a1;Parent=t1;protein_id=S1\n\
             c1\tx\tCDS\t100\t111\t.\t+\t0\tID=a2;Parent=t1;protein_id=S1\n\
             c1\tx\tCDS\t1\t13\t.\t+\t0\tID=b1;Parent=t2;protein_id=S2\n\
             c1\tx\tCDS\t100\t111\t.\t+\t1\tID=b2;Parent=t2;protein_id=S2\n",
        )
        .unwrap();
        let faa = dir.join("p.faa");
        std::fs::write(&faa, ">S1\nMKVLGAWETSIGL\n>S2\nMKVLGAWETSIGL\n").unwrap();
        let aln = dir.join("p.fasta");
        std::fs::write(&aln, ">S1\nMKVLGAWETSIGL\n>S2\nMKVLGAWETSIGL\n").unwrap();

        let strict = SiteMatrix::build(&run_projection(&gff, &faa, &aln, false));
        assert_eq!(strict.n_sites(), 2, "12:0 and 13:1 are different positions");
        let loose = SiteMatrix::build(&run_projection(&gff, &faa, &aln, true));
        assert_eq!(
            loose.n_sites(),
            1,
            "--ignore-phase collapses them onto one column"
        );
    }

    fn run_projection(
        gff: &Path,
        faa: &Path,
        aln: &Path,
        ignore_phase: bool,
    ) -> AlignmentProjection {
        let mut set = AnnotationSet::default();
        GffSource::new(vec![gff.to_path_buf()], GffOptions::default())
            .stream_into(&mut set)
            .unwrap();
        let proteins = ProteinSet::read([faa]).unwrap();
        let proj = join_annotations(&set.transcripts, &proteins, &JoinOptions::default()).unwrap();
        let a = Alignment::read(aln).unwrap();
        project_into_alignment(
            &a,
            &proj,
            &proteins,
            &project_options(ignore_phase, false, 1),
        )
    }
}
