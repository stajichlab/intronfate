//! # intronfate-core
//!
//! Project intron positions from genome annotation into **protein space**, then into
//! **multiple-sequence-alignment (MSA) columns**, and score how those positions are shared
//! across homologs with parsimony frameworks (Dollo, unrestricted maximum parsimony, counts).
//!
//! ## The coordinate model
//!
//! An intron is a gap between two consecutive *coding* segments of a transcript. Walking the
//! CDS in transcript order (5'→3' on the transcribed strand) and counting coding nucleotides
//! gives each intron a unique `cds_offset`: how many coding bases lie 5' of it. That single
//! integer is lossless and, because protein sequence is a `3`-to-`1` read of coding sequence,
//! it splits into the classical comparative-genomics pair *(codon, phase)*:
//!
//! ```text
//! codon_index0 = cds_offset / 3     // 0-based codon the intron anchors to
//! phase        = cds_offset % 3     // coding bases of the interrupted codon that precede it
//! ```
//!
//! * `phase 0` — intron sits exactly at a codon boundary (between the complete codons
//!   `cds_offset/3` and `cds_offset/3 + 1`).
//! * `phase 1` / `phase 2` — intron splits codon `cds_offset/3 + 1` after 1 or 2 bases.
//!
//! [`ProteinSite`] carries all of these views and every table emits them as separate columns
//! (`cds_offset`, `phase`, `residue_before`, `residue_after`, `codon`, `site_key`) so downstream
//! code can pick its own convention instead of re-deriving it.
//!
//! Two introns in different species are *homologous* when their protein-space positions land in
//! the **same MSA column** — that is the whole point of projecting into the alignment: residue
//! numbering drifts with indels, alignment columns do not.
//!
//! ## Pipeline
//!
//! ```text
//! GFF3/GTF ──► AnnotationSource ──► TranscriptAnnotation (genomic introns + cds_offset)
//!                                        │
//! protein FASTA ──► ProteinSet ──────────┤  join on protein id
//!                                        ▼
//!                              ProteinSite per intron
//!                                        │  MSA row: residue index ──► column
//!                                        ▼
//!                              SiteMatrix (column × taxon, 0/1)
//!                                        │
//!                 ┌──────────────────────┼───────────────────────┐
//!                 ▼                      ▼                       ▼
//!        counts / pairwise      Dollo parsimony (Newick)   maximum parsimony (Fitch)
//! ```
//!
//! ## Feature areas
//!
//! * [`gff`] / [`annotation`] — streaming, format-tolerant GFF3 ingestion behind the
//!   [`AnnotationSource`](annotation::AnnotationSource) trait, so a SQL/Chado/Ensembl reader can
//!   be dropped in without touching anything downstream.
//! * [`msa`] — aligned FASTA / Stockholm with residue↔column maps.
//! * [`project`] — protein sites to alignment columns, with tail and phase policies.
//! * [`matrix`] — dense bitset `column × taxon` matrix, pairwise sharing, Jaccard.
//! * [`tree`] — Newick reading and re-rooting.
//! * [`score`] — [`Dollo`](score::Dollo) and [`MaxParsimony`](score::MaxParsimony) scoring plus
//!   ancestral state reconstruction.
//! * [`tables`] — TSV/NEXUS writers.

#![warn(missing_debug_implementations)]
#![warn(clippy::all)]

pub mod annotation;
pub mod error;
pub mod gff;
pub mod inputs;
pub mod likelihood;
pub mod malin;
pub mod matrix;
pub mod model;
pub mod msa;
pub mod project;
pub mod report;
pub mod score;
pub mod seqio;
pub mod tables;
pub mod tree;
pub mod util;

pub use error::{Error, Result};
pub use model::{AnnotatedIntron, Phase, ProteinSite, Strand, TranscriptAnnotation};
