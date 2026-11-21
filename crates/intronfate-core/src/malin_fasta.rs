//! MALIN's own input: one aligned protein FASTA per gene family, with the organism and the intron
//! positions in each header (SPEC §6c).
//!
//! ```text
//! >gene1 /organism=Hs {i 23,48 i}
//! MKV-LLA...
//! ```
//!
//! `p` in `{i p1,p2 i}` is the number of coding bases 5' of the intron (1 = after the first nucleotide),
//! which is this tool's `cds_offset`. The tags may come in either order. A file becomes the three inputs
//! the rest of the pipeline already takes: an [`Alignment`], a [`ProteinSet`] of the ungapped rows, and one
//! [`TranscriptAnnotation`] per row, plus a row-id → organism map for the taxon names.

use rustc_hash::FxHashMap;

use crate::error::{Error, Result};
use crate::model::{AnnotatedIntron, Strand, TranscriptAnnotation};
use crate::msa::Alignment;
use crate::seqio::{is_gap, ProteinSet};

/// One alignment row with its header tags.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MalinRow {
    pub id: String,
    pub organism: String,
    /// Intron positions, sorted, without repeats.
    pub introns: Vec<u32>,
    /// The aligned sequence, gaps included.
    pub aligned: String,
}

impl MalinRow {
    fn ungapped(&self) -> String {
        self.aligned
            .bytes()
            .filter(|&b| !is_gap(b))
            .map(char::from)
            .collect()
    }
}

/// A parsed MALIN FASTA file.
#[derive(Clone, Debug)]
pub struct MalinFasta {
    pub name: String,
    pub rows: Vec<MalinRow>,
}

/// `/organism=NAME` from the header remainder, `None` when absent.
fn organism_tag(rest: &str) -> Option<&str> {
    let at = rest.find("/organism=")? + "/organism=".len();
    rest[at..]
        .split_whitespace()
        .next()
        .filter(|s| !s.is_empty())
}

/// The positions inside `{i ... i}`. `Ok(None)` when the header has no such tag.
fn intron_tag(rest: &str, row: &str) -> Result<Option<Vec<u32>>> {
    let Some(open) = rest.find("{i") else {
        return Ok(None);
    };
    let body = &rest[open + 2..];
    let Some(close) = body.find("i}") else {
        return Err(Error::parse(format!(
            "row {row:?}: `{{i` has no closing `i}}` in the header"
        )));
    };
    let mut out = Vec::new();
    for tok in body[..close].split(|c: char| c == ',' || c.is_whitespace()) {
        if tok.is_empty() {
            continue;
        }
        let p: u32 = tok.parse().map_err(|_| {
            Error::parse(format!(
                "row {row:?}: intron position {tok:?} is not a whole number"
            ))
        })?;
        if p == 0 {
            return Err(Error::parse(format!(
                "row {row:?}: intron position 0 is not valid (positions count the coding bases before the intron, so the first is 1)"
            )));
        }
        out.push(p);
    }
    out.sort_unstable();
    out.dedup();
    Ok(Some(out))
}

impl MalinFasta {
    /// Parse one aligned FASTA. `name` becomes the alignment name.
    pub fn parse(text: &str, name: impl Into<String>) -> Result<MalinFasta> {
        let name = name.into();
        let mut rows: Vec<MalinRow> = Vec::new();
        let mut cur: Option<MalinRow> = None;
        for line in text.lines() {
            let line = line.trim_end_matches('\r');
            if line.trim().is_empty() {
                continue;
            }
            if let Some(h) = line.strip_prefix('>') {
                rows.extend(cur.take());
                let h = h.trim();
                let id = h.split_whitespace().next().unwrap_or("");
                if id.is_empty() || id.starts_with('/') || id.starts_with("{i") {
                    return Err(Error::parse(format!(
                        "alignment {name:?}: header {line:?} must start with the row identifier, before the /organism= and {{i … i}} tags"
                    )));
                }
                let rest = &h[id.len()..];
                let organism = organism_tag(rest).ok_or_else(|| {
                    Error::parse(format!(
                        "alignment {name:?}: row {id:?} has no /organism= tag"
                    ))
                })?;
                let introns = intron_tag(rest, id)?.ok_or_else(|| {
                    Error::parse(format!(
                        "alignment {name:?}: row {id:?} has no {{i … i}} intron tag (write `{{i i}}` for a gene with no introns)"
                    ))
                })?;
                cur = Some(MalinRow {
                    id: id.to_string(),
                    organism: organism.to_string(),
                    introns,
                    aligned: String::new(),
                });
                continue;
            }
            match cur.as_mut() {
                Some(r) => r
                    .aligned
                    .extend(line.chars().filter(|c| !c.is_whitespace())),
                None => {
                    return Err(Error::parse(format!(
                        "alignment {name:?}: sequence data before any header"
                    )))
                }
            }
        }
        rows.extend(cur);
        if rows.is_empty() {
            return Err(Error::parse(format!(
                "alignment {name:?}: no FASTA records"
            )));
        }
        Ok(MalinFasta { name, rows })
    }

    /// The protein alignment (row ids are the header identifiers).
    pub fn alignment(&self) -> Result<Alignment> {
        let mut text = String::new();
        for r in &self.rows {
            text.push('>');
            text.push_str(&r.id);
            text.push('\n');
            text.push_str(&r.aligned);
            text.push('\n');
        }
        Alignment::from_fasta_str(&text, self.name.clone())
    }

    /// The ungapped rows as a protein set.
    pub fn proteins(&self) -> Result<ProteinSet> {
        let mut text = String::new();
        for r in &self.rows {
            text.push('>');
            text.push_str(&r.id);
            text.push('\n');
            text.push_str(&r.ungapped());
            text.push('\n');
        }
        let mut set = ProteinSet::new();
        set.read_one(text.as_bytes(), std::path::Path::new(&self.name))
            .map_err(|e| Error::parse(format!("alignment {:?}: {e}", self.name)))?;
        Ok(set)
    }

    /// One annotation per row. There are no genomic coordinates in this source: each intron carries
    /// start = end = 0 and the sequence name is `.`.
    pub fn transcripts(&self) -> Vec<TranscriptAnnotation> {
        self.rows
            .iter()
            .map(|r| {
                let n = r.ungapped().len() as u64;
                TranscriptAnnotation {
                    transcript_id: r.id.clone(),
                    gene_id: None,
                    protein_id: Some(r.id.clone()),
                    seqid: ".".to_string(),
                    strand: Strand::Unknown,
                    cds_len: 3 * n,
                    cds_segments: 1,
                    introns: r
                        .introns
                        .iter()
                        .enumerate()
                        .map(|(rank, &p)| AnnotatedIntron {
                            start: 0,
                            end: 0,
                            rank: rank as u32,
                            cds_offset: Some(p),
                        })
                        .collect(),
                    frame_mismatches: 0,
                    source: self.name.clone(),
                }
            })
            .collect()
    }

    /// Row id → organism, the taxon names of the matrix.
    pub fn taxon_names(&self) -> FxHashMap<String, String> {
        self.rows
            .iter()
            .map(|r| (r.id.clone(), r.organism.clone()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GUIDE: &str =
        ">a1 /organism=Hs {i 3,9 i}\nMKV-LLA\nGT\n>b1 {i 6 i} /organism=Dm\nMKVALLAGT\n";

    #[test]
    fn parses_tags_in_either_order() {
        let f = MalinFasta::parse(GUIDE, "fam").unwrap();
        assert_eq!(f.rows.len(), 2);
        assert_eq!(f.rows[0].organism, "Hs");
        assert_eq!(f.rows[0].introns, vec![3, 9]);
        assert_eq!(f.rows[0].aligned, "MKV-LLAGT");
        assert_eq!(f.rows[1].organism, "Dm");
        assert_eq!(f.rows[1].introns, vec![6]);
    }

    #[test]
    fn empty_tag_spaces_crlf_and_repeats() {
        let t = ">x /organism=A {i i}\r\nMKV\r\n>y /organism=B {i 9,  3 , 3 i}\r\nMKVLLA\r\n";
        let f = MalinFasta::parse(t, "f").unwrap();
        assert!(f.rows[0].introns.is_empty());
        assert_eq!(f.rows[0].aligned, "MKV");
        assert_eq!(f.rows[1].introns, vec![3, 9]);
    }

    #[test]
    fn refusals_name_the_row() {
        let e = MalinFasta::parse(">x {i 3 i}\nMKV\n", "f")
            .unwrap_err()
            .to_string();
        assert!(e.contains("/organism="), "{e}");
        let e = MalinFasta::parse(">x /organism=A\nMKV\n", "f")
            .unwrap_err()
            .to_string();
        assert!(e.contains("intron tag"), "{e}");
        let e = MalinFasta::parse(">x /organism=A {i 0 i}\nMKV\n", "f")
            .unwrap_err()
            .to_string();
        assert!(e.contains("position 0"), "{e}");
        let e = MalinFasta::parse(">x /organism=A {i a i}\nMKV\n", "f")
            .unwrap_err()
            .to_string();
        assert!(e.contains("whole number"), "{e}");
        let e = MalinFasta::parse(">x /organism=A {i 3\nMKV\n", "f")
            .unwrap_err()
            .to_string();
        assert!(e.contains("closing"), "{e}");
        let e = MalinFasta::parse(">/organism=A {i 3 i}\nMKV\n", "f")
            .unwrap_err()
            .to_string();
        assert!(e.contains("identifier"), "{e}");
        assert!(MalinFasta::parse("MKV\n", "f").is_err());
        assert!(MalinFasta::parse("", "f").is_err());
    }

    #[test]
    fn builds_the_three_inputs() {
        let f = MalinFasta::parse(GUIDE, "fam").unwrap();
        let aln = f.alignment().unwrap();
        assert_eq!((aln.n_rows(), aln.len), (2, 9));
        let prot = f.proteins().unwrap();
        assert_eq!(prot.len(), 2);
        assert_eq!(prot.residues(0).unwrap(), b"MKVLLAGT");
        let tx = f.transcripts();
        assert_eq!(tx[0].cds_len, 24);
        assert_eq!(tx[0].introns[1].cds_offset, Some(9));
        assert_eq!(f.taxon_names()["b1"], "Dm");
    }
}
