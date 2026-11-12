//! Protein alignments and the residue → column maps that make projection possible.
//!
//! An MSA column, not a residue number, is the homology unit: residue 42 of one ortholog and
//! residue 39 of another can be the *same* position, and only the alignment says so. Each row
//! therefore precomputes `res2col`, the alignment column of its n-th residue, which turns
//! projection into an array lookup.
//!
//! Stockholm insert columns matter here. In a Pfam-style alignment the row may carry extra
//! lowercase insert columns that are **not** part of the source sequence, and `#=GS <id> SF/SE`
//! says which source range the row covers. Ignoring that silently shifts every downstream residue
//! number, so it is honoured when present and reported when absent.

use std::io::BufRead;
use std::path::PathBuf;

use rustc_hash::FxHashMap;

use crate::error::{Error, Result};
use crate::seqio::{is_gap, is_residue};
use crate::util;

/// Formats accepted by [`Alignment::read`].
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AlignmentFormat {
    Fasta,
    Stockholm,
}

impl std::fmt::Display for AlignmentFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            AlignmentFormat::Fasta => "fasta",
            AlignmentFormat::Stockholm => "stockholm",
        })
    }
}

/// One aligned sequence.
#[derive(Clone, Debug)]
pub struct AlignmentRow {
    pub id: String,
    /// Aligned characters, gaps included.
    pub aligned: Vec<u8>,
    /// Alignment column of each residue present in the row, in order.
    pub res2col: Vec<u32>,
    /// Index in the *source protein* of `res2col[0]`. Stockholm `SF` omits the N-terminal part of a
    /// sequence from the alignment, so protein residue numbers start here rather than at 0.
    pub residue_offset: usize,
    /// How many entries of `res2col` belong to the source range (`SF`..`SE`).
    pub numbered: usize,
}

impl AlignmentRow {
    /// Residues of the source protein that this row actually covers.
    pub fn n_residues(&self) -> usize {
        self.numbered
    }

    /// Column holding `residue` (0-based in the source protein), `None` when the row does not
    /// cover it — an out-of-window position is unknown, never silently clamped.
    pub fn column_of(&self, residue: usize) -> Option<usize> {
        let q = residue.checked_sub(self.residue_offset)?;
        if q >= self.numbered {
            return None;
        }
        self.res2col.get(q).map(|&c| c as usize)
    }

    /// Column of the last covered residue, used by `TailPolicy::Clamp`.
    ///
    /// `res2col` is indexed by *aligned* residue, while `residue_offset` is a position in the
    /// *source*, so the last covered residue is simply entry `numbered - 1`.
    pub fn last_column(&self) -> Option<usize> {
        self.res2col
            .get(self.numbered.saturating_sub(1))
            .map(|&c| c as usize)
    }

    /// Last column at or before `col` holding a residue of this row.
    pub fn prev_residue_col(&self, col: usize) -> Option<usize> {
        self.res2col
            .binary_search_by(|c| (*c as usize).cmp(&col))
            .map_or_else(
                |i| i.checked_sub(1).map(|i| self.res2col[i] as usize),
                |i| Some(self.res2col[i] as usize),
            )
    }

    /// First column at or after `col` holding a residue of this row.
    pub fn next_residue_col(&self, col: usize) -> Option<usize> {
        let i = self
            .res2col
            .binary_search_by(|c| (*c as usize).cmp(&col))
            .unwrap_or_else(|i| i);
        self.res2col.get(i).map(|&c| c as usize)
    }

    pub fn n_gaps(&self) -> usize {
        self.aligned.len() - self.res2col.len()
    }
}

/// A multiple sequence alignment.
#[derive(Clone, Debug)]
pub struct Alignment {
    pub rows: Vec<AlignmentRow>,
    /// Alignment length in columns.
    pub len: usize,
    /// Label used in tables; the file stem unless overridden.
    pub name: String,
    pub format: AlignmentFormat,
    pub source: Option<PathBuf>,
    index: FxHashMap<String, usize>,
    /// Rows that were all gaps, or had no source residues.
    pub gap_only_rows: Vec<String>,
    /// `#=GF` (ID, value) pairs from a Stockholm file.
    pub gf: Vec<(String, String)>,
}

impl Alignment {
    /// Read an alignment, sniffing FASTA vs Stockholm from the first line.
    pub fn read(path: impl Into<PathBuf>) -> Result<Alignment> {
        let path = path.into();
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.strip_suffix(".gz").unwrap_or(s).to_string())
            .unwrap_or_else(|| path.display().to_string());
        let reader = util::open_buffered(&path, 1 << 18)?;
        let mut aln = Self::read_reader(reader, name)?;
        aln.source = Some(path);
        Ok(aln)
    }

    pub fn read_reader<R: BufRead>(reader: R, name: String) -> Result<Alignment> {
        let mut lines: Vec<String> = Vec::new();
        let mut br = reader;
        let mut scratch: Vec<u8> = Vec::new();
        loop {
            let mut line = String::new();
            if util::read_line_lossy(&mut br, &mut line, &mut scratch)?.0 == 0 {
                break;
            }
            lines.push(line);
        }
        let is_stockholm = lines
            .iter()
            .find(|l| !l.trim().is_empty())
            .map(|l| l.trim_start().starts_with("# STOCKHOLM"))
            .unwrap_or(false);
        if is_stockholm {
            Self::from_stockholm_lines(&lines, name)
        } else {
            Self::from_fasta_lines(&lines, name)
        }
    }

    pub fn from_fasta_str(text: &str, name: impl Into<String>) -> Result<Alignment> {
        Self::from_fasta_lines(
            &text.lines().map(|l| l.to_string()).collect::<Vec<_>>(),
            name.into(),
        )
    }

    pub fn from_fasta_lines(lines: &[String], name: String) -> Result<Alignment> {
        let mut raw: Vec<(String, Vec<u8>)> = Vec::new();
        let mut cur: Option<(String, Vec<u8>)> = None;
        for l in lines {
            let l = l.trim_end_matches(['\n', '\r']);
            if l.is_empty() {
                continue;
            }
            if let Some(h) = l.strip_prefix('>') {
                if let Some(c) = cur.take() {
                    raw.push(c);
                }
                let id = h.split_whitespace().next().unwrap_or("").to_string();
                if id.is_empty() {
                    return Err(Error::parse("FASTA header has no identifier"));
                }
                cur = Some((id, Vec::new()));
                continue;
            }
            match cur.as_mut() {
                Some((_, s)) => s.extend(l.bytes().filter(|b| !b.is_ascii_whitespace())),
                None => {
                    return Err(Error::parse(format!(
                        "sequence data before any header in alignment {name:?}: {l:?}"
                    )))
                }
            }
        }
        if let Some(c) = cur {
            raw.push(c);
        }
        Self::assemble(raw, name, AlignmentFormat::Fasta, &Default::default())
    }

    pub fn from_stockholm_str(text: &str, name: impl Into<String>) -> Result<Alignment> {
        Self::from_stockholm_lines(
            &text.lines().map(|l| l.to_string()).collect::<Vec<_>>(),
            name.into(),
        )
    }

    /// Parse Stockholm, using `#=GS <id> SF/SE` to exclude insert columns from residue numbering.
    pub fn from_stockholm_lines(lines: &[String], name: String) -> Result<Alignment> {
        // A Stockholm file may interleave blocks, so accumulate per id in first-seen order.
        let mut chunks: FxHashMap<String, Vec<u8>> = Default::default();
        let mut order: Vec<String> = Vec::new();
        let mut gf = Vec::new();
        let mut bounds: FxHashMap<String, (Option<u64>, Option<u64>)> = Default::default();
        for l in lines {
            let l = l.trim_end_matches(['\n', '\r']);
            let t = l.trim_start();
            if t.is_empty() {
                continue;
            }
            if t == "//" {
                break;
            }
            if let Some(rest) = t.strip_prefix("#=GS") {
                // `#=GS <id> <fc> <value>`
                let mut it = rest.split_whitespace();
                let (Some(id), Some(fc)) = (it.next(), it.next()) else {
                    continue;
                };
                let value = it.next();
                let entry = bounds.entry(id.to_string()).or_default();
                match fc.to_ascii_uppercase().as_str() {
                    "SF" => entry.0 = value.and_then(|v| v.parse().ok()),
                    "SE" => entry.1 = value.and_then(|v| v.parse().ok()),
                    _ => {}
                }
                continue;
            }
            if t.starts_with("#=GF") {
                let mut it = t.splitn(3, char::is_whitespace);
                let (Some(_), Some(k)) = (it.next(), it.next()) else {
                    continue;
                };
                gf.push((k.to_string(), it.next().unwrap_or("").to_string()));
                continue;
            }
            if t.starts_with('#') || t.starts_with("//") {
                continue;
            }
            let Some((id, seq)) = t.split_once(char::is_whitespace) else {
                return Err(Error::parse(format!(
                    "malformed Stockholm line in {name:?}: {t:?}"
                )));
            };
            let buf = chunks.entry(id.to_string()).or_default();
            if buf.is_empty() {
                order.push(id.to_string());
            }
            buf.extend(seq.bytes().filter(|b| !b.is_ascii_whitespace()));
        }
        let raw: Vec<(String, Vec<u8>)> = order
            .into_iter()
            .map(|id| {
                let seq = chunks.remove(&id).unwrap_or_default();
                (id, seq)
            })
            .collect();
        let mut aln = Self::assemble(raw, name, AlignmentFormat::Stockholm, &bounds)?;
        aln.gf = gf;
        Ok(aln)
    }

    fn assemble(
        raw: Vec<(String, Vec<u8>)>,
        name: String,
        format: AlignmentFormat,
        bounds: &FxHashMap<String, (Option<u64>, Option<u64>)>,
    ) -> Result<Alignment> {
        if raw.is_empty() {
            return Err(Error::parse(format!(
                "alignment {name:?} contains no sequences"
            )));
        }
        let len = raw[0].1.len();
        if len == 0 {
            return Err(Error::parse(format!(
                "alignment {name:?} has zero-length rows"
            )));
        }
        let mut bad: Vec<String> = Vec::new();
        for (id, s) in &raw {
            if s.len() != len {
                bad.push(format!("{id}: {}", s.len()));
            }
        }
        if !bad.is_empty() {
            return Err(Error::parse(format!(
                "alignment {name:?} rows must all be {len} columns, but {} row(s) differ [{}]",
                bad.len(),
                util::ellipsis(&bad.join(", "), 300)
            )));
        }
        let mut ids: FxHashMap<String, usize> = Default::default();
        let mut rows = Vec::with_capacity(raw.len());
        let mut gap_only_rows = Vec::new();
        for (id, aligned) in raw {
            if ids.contains_key(&id) {
                return Err(Error::parse(format!(
                    "alignment {name:?}: duplicate sequence id {id:?}"
                )));
            }
            let all: Vec<u32> = {
                let mut v = Vec::with_capacity(aligned.len());
                for (c, &b) in aligned.iter().enumerate() {
                    if is_residue(b) {
                        v.push(c as u32);
                    }
                }
                v
            };
            // SF/SE give the source range covered by this row (1-based, inclusive). Residues of the
            // protein outside it are absent from the alignment, so they must have no column.
            let (lead, keep) = match bounds.get(&id) {
                Some((sf, se)) => {
                    // Aligned residue k is source residue `SF + k`; the row covers SF..SE, so the
                    // window length is SE-SF+1 (bounded by what the row actually holds).
                    let lead = sf.map(|v| (v as usize).saturating_sub(1)).unwrap_or(0);
                    // SE is the source position of the last aligned residue, so the covered span is
                    // SE-SF+1 = SE-lead. The only cap is how many residues the row actually holds:
                    // `lead` is a position in the *source*, not something to subtract from it.
                    let span = se
                        .map(|se| (se as usize).saturating_sub(lead))
                        .unwrap_or(all.len());
                    (lead, span.min(all.len()))
                }
                None => (0, all.len()),
            };
            let res2col = all;
            let numbered = keep;
            if numbered == 0 {
                gap_only_rows.push(id.clone());
            }
            ids.insert(id.clone(), rows.len());
            rows.push(AlignmentRow {
                id,
                aligned,
                res2col,
                residue_offset: lead,
                numbered,
            });
        }
        Ok(Alignment {
            rows,
            len,
            name,
            format,
            source: None,
            index: ids,
            gap_only_rows,
            gf: Vec::new(),
        })
    }

    pub fn n_rows(&self) -> usize {
        self.rows.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = &str> {
        self.rows.iter().map(|r| r.id.as_str())
    }

    pub fn taxon_index(&self, id: &str) -> Option<usize> {
        self.index.get(id).copied()
    }

    pub fn row(&self, idx: usize) -> Option<&AlignmentRow> {
        self.rows.get(idx)
    }

    /// Characters of column `col`, one per row.
    pub fn column(&self, col: usize) -> impl Iterator<Item = u8> + '_ {
        self.rows
            .iter()
            .map(move |r| r.aligned.get(col).copied().unwrap_or(b'-'))
    }

    /// The consensus residue of a column (most frequent non-gap), `None` when fully gapped.
    pub fn consensus(&self, col: usize) -> Option<u8> {
        let mut counts = [0u32; 256];
        let mut best = 0u32;
        let mut letter = None;
        for b in self.column(col) {
            if is_gap(b) || !is_residue(b) {
                continue;
            }
            let u = b.to_ascii_uppercase() as usize;
            counts[u] += 1;
            if counts[u] > best {
                best = counts[u];
                letter = Some(b.to_ascii_uppercase());
            }
        }
        letter
    }

    /// `true` when no row has a residue in this column.
    pub fn column_is_empty(&self, col: usize) -> bool {
        self.column(col).all(|b| is_gap(b) || !is_residue(b))
    }

    /// Source-protein residue count per row.
    pub fn residues_aligned(&self) -> Vec<usize> {
        self.rows.iter().map(|r| r.n_residues()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_residue_to_column_maps() {
        let a = Alignment::from_fasta_str(">A\nMK-VL--A\n>B\nM-VLKA--\n", "t").unwrap();
        assert_eq!(a.len, 8);
        assert_eq!(a.n_rows(), 2);
        let a0 = a.row(0).unwrap();
        assert_eq!(a0.n_residues(), 5);
        assert_eq!(a0.column_of(0), Some(0));
        assert_eq!(a0.column_of(1), Some(1));
        assert_eq!(a0.column_of(2), Some(3), "the gap at column 2 is skipped");
        assert_eq!(a0.column_of(4), Some(7));
        assert_eq!(a0.column_of(5), None, "past the end of the sequence");
        assert_eq!(a.row(1).unwrap().column_of(0), Some(0));
        assert_eq!(a.row(1).unwrap().column_of(1), Some(2));
    }

    #[test]
    fn neighbour_columns_around_a_boundary() {
        let a = Alignment::from_fasta_str(">A\nMK--VL\n", "t").unwrap();
        let r = a.row(0).unwrap();
        assert_eq!(
            r.prev_residue_col(2),
            Some(1),
            "intron before a gap block anchors left"
        );
        assert_eq!(r.next_residue_col(2), Some(4));
        assert_eq!(r.prev_residue_col(0), Some(0));
        assert_eq!(r.next_residue_col(5), Some(5));
        assert_eq!(r.next_residue_col(6), None);
        assert_eq!(r.prev_residue_col(9), Some(5));
    }

    #[test]
    fn ragged_rows_are_rejected() {
        let e = Alignment::from_fasta_str(">A\nMKVL\n>B\nMKV\n", "t").unwrap_err();
        assert!(e.to_string().contains("must all be 4 columns"), "{e}");
        assert!(e.to_string().contains("B: 3"), "{e}");
    }

    #[test]
    fn duplicate_ids_rejected() {
        let e = Alignment::from_fasta_str(">A\nMK\n>A\nVL\n", "t").unwrap_err();
        assert!(e.to_string().contains("duplicate"), "{e}");
    }

    #[test]
    fn empty_alignment_rejected() {
        assert!(Alignment::from_fasta_str("", "t").is_err());
        assert!(Alignment::from_fasta_str(">A\n", "t").is_err());
        assert!(Alignment::from_fasta_str("SEQUENCEWITHOUTHEADER\n", "t").is_err());
    }

    #[test]
    fn stockholm_sf_se_defines_where_protein_numbering_starts() {
        // `SF 3 / SE 6` means this row carries source residues 3..6: the protein's first two
        // residues are absent from the alignment, so its aligned residue 0 is protein residue 3.
        // Getting this wrong shifts every projected intron position in the row by two.
        let sto = "# STOCKHOLM 1.0\n\
                   #=GS A SF 3\n\
                   #=GS A SE 6\n\
                   A      MKVL--\n\
                   B      -MKVL-\n\
                   #=GC RF  xxxxxx\n\
                   //\n";
        let a = Alignment::from_stockholm_str(sto, "s").unwrap();
        assert_eq!(a.format, AlignmentFormat::Stockholm);
        assert_eq!(a.len, 6);
        let r = a.row(0).unwrap();
        assert_eq!((r.residue_offset, r.n_residues()), (2, 4));
        assert_eq!(
            r.column_of(0),
            None,
            "source residue 1 is not in the alignment"
        );
        assert_eq!(r.column_of(1), None, "nor is residue 2");
        assert_eq!(
            r.column_of(2),
            Some(0),
            "residue 3 is the first aligned residue"
        );
        assert_eq!(r.column_of(5), Some(3));
        assert_eq!(
            r.column_of(6),
            None,
            "past SE there is nothing to anchor to"
        );
        assert_eq!(r.last_column(), Some(3));
        // Row B carries its whole sequence.
        assert_eq!(
            (
                a.row(1).unwrap().residue_offset,
                a.row(1).unwrap().n_residues()
            ),
            (0, 4)
        );
        assert_eq!(a.row(1).unwrap().column_of(0), Some(1));
    }

    #[test]
    fn stockholm_without_sf_covers_every_residue() {
        let sto = "# STOCKHOLM 1.0\nA  MKVL\nB  MK-L\n//\n";
        let a = Alignment::from_stockholm_str(sto, "s").unwrap();
        assert_eq!(a.row(0).unwrap().n_residues(), 4);
        assert_eq!(a.row(1).unwrap().n_residues(), 3);
    }

    #[test]
    fn stockholm_block_rows_are_concatenated() {
        let sto = "# STOCKHOLM 1.0\nA  MK\nB  VL\n\nA  VL\nB  MK\n//\n";
        let a = Alignment::from_stockholm_str(sto, "s").unwrap();
        assert_eq!(a.row(0).unwrap().aligned, b"MKVL");
        assert_eq!(a.row(1).unwrap().aligned, b"VLMK");
    }

    #[test]
    fn format_sniffing_from_content() {
        let sto = "# STOCKHOLM 1.0\nA  MK\n//\n";
        let a = Alignment::read_reader(sto.as_bytes(), "x".into()).unwrap();
        assert_eq!(a.format, AlignmentFormat::Stockholm);
        let fa = ">A\nMK\n";
        let a = Alignment::read_reader(fa.as_bytes(), "x".into()).unwrap();
        assert_eq!(a.format, AlignmentFormat::Fasta);
    }

    #[test]
    fn column_queries_and_gap_only_rows() {
        let a = Alignment::from_fasta_str(">A\nMK--\n>B\n-V--\n>C\n----\n", "t").unwrap();
        assert_eq!(a.column(0).collect::<Vec<_>>(), vec![b'M', b'-', b'-']);
        assert_eq!(a.consensus(0), Some(b'M'));
        assert!(a.column_is_empty(2));
        assert_eq!(a.gap_only_rows, vec!["C"]);
        assert_eq!(a.taxon_index("B"), Some(1));
        assert_eq!(a.taxon_index("Z"), None);
    }
}
