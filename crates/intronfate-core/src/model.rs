//! Core data model: strands, introns, and intron positions in protein space.

use std::fmt;

/// Which strand a feature was annotated on.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash, Default)]
pub enum Strand {
    Plus,
    Minus,
    /// GFF writes `.` when the strand is unknown.
    #[default]
    Unknown,
}

impl Strand {
    pub fn parse(b: u8) -> Option<Strand> {
        match b {
            b'+' => Some(Strand::Plus),
            b'-' => Some(Strand::Minus),
            b'.' | b'?' => Some(Strand::Unknown),
            _ => None,
        }
    }

    /// `1` for plus, `-1` for minus, `0` for unknown.
    pub fn sign(self) -> i64 {
        match self {
            Strand::Plus => 1,
            Strand::Minus => -1,
            Strand::Unknown => 0,
        }
    }

    pub fn is_known(self) -> bool {
        !matches!(self, Strand::Unknown)
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Unknown => ".",
        })
    }
}

/// Codon phase of an intron: how many coding bases of the interrupted codon precede it.
pub type Phase = u8;

/// An intron position expressed in protein space.
///
/// `cds_offset` is the number of coding nucleotides 5' of the intron within the CDS (0-based),
/// which is enough to recover both the codon and the phase; the two are stored together so the
/// value is self-describing in tables and APIs.
///
/// ```
/// use intronfate_core::model::ProteinSite;
/// // 15 coding bases before the intron -> exactly after complete codon 5, boundary intron.
/// let s = ProteinSite::from_cds_offset(15);
/// assert_eq!(s.phase, 0);
/// assert_eq!(s.codon(), 5);
/// assert_eq!((s.residue_before(), s.residue_after()), (5, 6));
/// assert_eq!(s.to_string(), "5:0");
///
/// // 37 -> inside codon 13, after its first base.
/// let s = ProteinSite::from_cds_offset(37);
/// assert_eq!(s.phase, 1);
/// assert_eq!(s.codon(), 13);
/// assert_eq!((s.residue_before(), s.residue_after()), (13, 13));
/// assert_eq!(s.to_string(), "13:1");
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct ProteinSite {
    /// 0-based count of coding nucleotides lying 5' of the intron.
    pub cds_offset: u32,
    /// `cds_offset % 3` — 0 for boundary introns, 1 or 2 for introns inside a codon.
    pub phase: Phase,
}

impl ProteinSite {
    /// Build a site from a CDS nucleotide offset, deriving the phase.
    pub const fn from_cds_offset(cds_offset: u32) -> Self {
        ProteinSite {
            cds_offset,
            phase: (cds_offset % 3) as Phase,
        }
    }

    pub const fn new(cds_offset: u32, phase: Phase) -> Self {
        ProteinSite { cds_offset, phase }
    }

    /// `true` when `cds_offset % 3 == phase`, i.e. the two views agree.
    pub const fn is_coherent(self) -> bool {
        self.phase as u32 == self.cds_offset % 3
    }

    /// An intron exactly at a codon boundary.
    pub const fn is_boundary(self) -> bool {
        self.phase == 0
    }

    /// 0-based index of the codon this site anchors to in the protein.
    ///
    /// Uniform across phases: for `phase > 0` it is the interrupted codon, for `phase == 0` it
    /// is the codon *following* the intron. That is also the index used to look the site up in an
    /// alignment row, so two homologous sites collide only when they share a phase as well.
    pub const fn codon_index0(self) -> u32 {
        self.cds_offset / 3
    }

    /// 1-based number of the last complete codon 5' of the intron (`0` if the intron precedes
    /// the start codon).
    pub const fn full_codons_before(self) -> u32 {
        self.cds_offset / 3
    }

    /// 1-based residue immediately 5' of the intron, or the interrupted residue for `phase > 0`.
    pub const fn residue_before(self) -> u32 {
        if self.phase == 0 {
            self.cds_offset / 3
        } else {
            self.cds_offset / 3 + 1
        }
    }

    /// 1-based residue immediately 3' of the intron (the interrupted residue for `phase > 0`).
    pub const fn residue_after(self) -> u32 {
        self.cds_offset / 3 + 1
    }

    /// 1-based codon number used in the `<codon>:<phase>` label.
    pub fn codon(self) -> u32 {
        self.residue_before().max(1)
    }

    /// Stable machine key, `cds_offset` as text.
    pub fn site_key(self) -> String {
        self.cds_offset.to_string()
    }

    /// ROPE-style `<codon>:<phase>` label, e.g. `13:1`.
    pub fn label(self) -> String {
        format!("{}:{}", self.codon(), self.phase)
    }
}

impl fmt::Display for ProteinSite {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.codon(), self.phase)
    }
}

/// One intron of a transcript, in genomic coordinates and (when it is a coding intron) in
/// protein space.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AnnotatedIntron {
    /// 1-based inclusive genomic start (lower coordinate).
    pub start: u64,
    /// 1-based inclusive genomic end (higher coordinate).
    pub end: u64,
    /// 0-based rank of this intron among all introns of the transcript, in transcript order.
    pub rank: u32,
    /// Coding offset, or `None` when the intron lies in UTR / non-coding sequence and therefore
    /// has no protein-space position.
    pub cds_offset: Option<u32>,
}

impl AnnotatedIntron {
    pub const fn length(self) -> u64 {
        self.end - self.start + 1
    }

    pub const fn is_coding(self) -> bool {
        self.cds_offset.is_some()
    }

    /// The protein-space position, `None` for non-coding introns.
    pub const fn site(self) -> Option<ProteinSite> {
        match self.cds_offset {
            Some(o) => Some(ProteinSite::from_cds_offset(o)),
            None => None,
        }
    }
}

/// A transcript's splicing structure, as produced by any [`AnnotationSource`](
/// crate::annotation::AnnotationSource) — GFF3 today, a database reader later.
#[derive(Clone, Debug)]
pub struct TranscriptAnnotation {
    pub transcript_id: String,
    pub gene_id: Option<String>,
    /// Identifier to join against the protein FASTA.
    pub protein_id: Option<String>,
    pub seqid: String,
    pub strand: Strand,
    /// Total coding length in bases, over all merged CDS segments.
    pub cds_len: u64,
    /// Number of CDS segments before merging contiguous neighbours.
    pub cds_segments: u32,
    /// Introns in transcript order.
    pub introns: Vec<AnnotatedIntron>,
    /// CDS `phase`/`frame` attributes that disagreed with the derived reading frame.
    pub frame_mismatches: u32,
    /// Where this record came from (file path, table name, …) for diagnostics.
    pub source: String,
}

impl TranscriptAnnotation {
    pub fn coding_introns(&self) -> impl Iterator<Item = &AnnotatedIntron> {
        self.introns.iter().filter(|i| i.is_coding())
    }

    pub fn n_coding_introns(&self) -> usize {
        self.coding_introns().count()
    }

    /// Protein-space positions of every coding intron, in transcript order.
    pub fn protein_sites(&self) -> impl Iterator<Item = ProteinSite> + '_ {
        self.introns.iter().filter_map(|i| i.site())
    }

    /// Expected protein length if `cds_len` includes a terminal stop codon.
    pub fn expected_protein_len_min(&self) -> u64 {
        self.cds_len.saturating_sub(3) / 3
    }

    pub fn expected_protein_len_max(&self) -> u64 {
        self.cds_len.div_ceil(3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_boundaries() {
        for off in 0..12u32 {
            let s = ProteinSite::from_cds_offset(off);
            assert!(s.is_coherent(), "{off}");
            assert_eq!(s.cds_offset, s.codon_index0() * 3 + s.phase as u32);
        }
        let s = ProteinSite::from_cds_offset(0);
        assert!(s.is_boundary());
        assert_eq!(s.residue_before(), 0);
        assert_eq!(s.residue_after(), 1);
        assert_eq!(s.label(), "1:0");

        let s = ProteinSite::from_cds_offset(2);
        assert_eq!((s.codon(), s.phase), (1, 2));
        assert_eq!((s.residue_before(), s.residue_after()), (1, 1));
    }

    #[test]
    fn coding_and_utr_introns() {
        let tx = TranscriptAnnotation {
            transcript_id: "t1".into(),
            gene_id: None,
            protein_id: None,
            seqid: "I".into(),
            strand: Strand::Plus,
            cds_len: 0,
            cds_segments: 0,
            introns: vec![
                AnnotatedIntron {
                    start: 10,
                    end: 20,
                    rank: 0,
                    cds_offset: Some(15),
                },
                AnnotatedIntron {
                    start: 30,
                    end: 40,
                    rank: 1,
                    cds_offset: None,
                },
            ],
            frame_mismatches: 0,
            source: "test".into(),
        };
        assert_eq!(tx.n_coding_introns(), 1);
        assert_eq!(
            tx.protein_sites().collect::<Vec<_>>(),
            vec![ProteinSite::from_cds_offset(15)]
        );
    }
}
