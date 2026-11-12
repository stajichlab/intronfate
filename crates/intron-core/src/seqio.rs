//! Protein FASTA reading and identifier resolution.
//!
//! The join between annotation and sequence is the fragile edge of this whole pipeline: a GFF may
//! say `protein_id=ENSP000002.3` while the FASTA header says `ENSP000002.3 cytochrome b [species]`,
//! or the two may disagree on versions. [`ProteinSet`] therefore keeps two indexes — exact and
//! version-stripped — and reports every sequence and every annotation record it could not join
//! instead of quietly dropping it.

use std::io::BufRead;
use std::path::{Path, PathBuf};

use rustc_hash::FxHashMap;

use crate::error::{Error, Result};
use crate::report::Level;
use crate::util;

/// A residue character in an aligned or unaligned protein.
pub fn is_residue(b: u8) -> bool {
    b.is_ascii_alphabetic() || b == b'*'
}

/// Gap characters used by the alignment formats we accept.
pub fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'.' | b'~')
}

/// How an identifier matched the sequence set.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Resolution {
    /// The id is present verbatim.
    Exact,
    /// Only the version-stripped id matched.
    VersionInsensitive,
    /// The version-stripped id matches several sequences, so no join is made.
    Ambiguous,
    /// No match.
    Missing,
}

/// One unaligned sequence.
#[derive(Clone, Debug)]
pub struct Sequence {
    pub id: String,
    /// Always `None`: header descriptions are not kept.
    pub description: Option<String>,
    /// Empty when the set was read with [`ProteinSet::read_lengths`]; use `len` for the length.
    pub residues: Vec<u8>,
    /// Number of residues, whether or not they were kept.
    pub len: usize,
}

/// A set of protein sequences with fast, forgiving lookup.
#[derive(Clone, Debug, Default)]
pub struct ProteinSet {
    pub seqs: Vec<Sequence>,
    exact: FxHashMap<String, usize>,
    /// Version-stripped id → index, only consulted when the exact lookup misses.
    /// `None` marks a stripped key shared by several sequences: refusing to guess beats joining the
    /// wrong paralog.
    ///
    /// Built on the first lookup that misses exactly (never, when every annotation id matches a header
    /// verbatim), because it holds a second copy of every id and a proteome set can be ten million long.
    stripped: std::sync::OnceLock<FxHashMap<String, Option<usize>>>,
    /// `(wanted_id, first_id_that_occupied_it)` for every duplicate header.
    pub duplicates: Vec<(String, String)>,
    pub malformed: u64,
    /// Lines that held invalid UTF-8 and were read lossily.
    pub lossy_lines: u64,
    pub empty_records: u64,
    pub total_residues: u64,
    pub sources: Vec<PathBuf>,
    strip_versions: bool,
    /// `false` for a lengths-only set: residues are counted, not stored.
    keep_residues: bool,
}

impl ProteinSet {
    /// Empty set with the forgiving (version-insensitive) join enabled.
    pub fn new() -> Self {
        ProteinSet {
            strip_versions: true,
            keep_residues: true,
            ..Default::default()
        }
    }

    /// Read protein FASTA files, keeping every residue.
    pub fn read(paths: impl IntoIterator<Item = impl AsRef<Path>>) -> Result<Self> {
        Self::read_with(paths, true)
    }

    /// Read ids and lengths only. A proteome is mostly sequence, so this is what `project` uses: it needs
    /// each protein's length for the CDS cross-check and nothing else.
    pub fn read_lengths(paths: impl IntoIterator<Item = impl AsRef<Path>>) -> Result<Self> {
        Self::read_with(paths, false)
    }

    /// `true` when residues are stored (not a lengths-only set).
    pub fn keeps_residues(&self) -> bool {
        self.keep_residues
    }

    /// Files are read in parallel and merged in the order given, so the result (ids, indices, duplicate
    /// reports) does not depend on thread timing.
    fn read_with(
        paths: impl IntoIterator<Item = impl AsRef<Path>>,
        keep_residues: bool,
    ) -> Result<Self> {
        use rayon::prelude::*;
        let paths: Vec<PathBuf> = paths
            .into_iter()
            .map(|p| p.as_ref().to_path_buf())
            .collect();
        let parts: Vec<Result<ProteinSet>> = paths
            .par_iter()
            .map(|path| {
                let mut part = ProteinSet {
                    strip_versions: true,
                    keep_residues,
                    ..Default::default()
                };
                let reader = util::open_buffered(path, 1 << 18)?;
                part.read_one(reader, path)?;
                Ok(part)
            })
            .collect();
        let mut set = ProteinSet {
            sources: paths.clone(),
            strip_versions: true,
            keep_residues,
            ..Default::default()
        };
        for part in parts {
            set.merge_from(part?);
        }
        if set.seqs.is_empty() {
            return Err(Error::EmptyInput {
                thing: "protein sequences",
                checked: paths.len(),
            });
        }
        Ok(set)
    }

    /// Append another set's sequences (moving them), re-indexing and reporting duplicates across the two.
    fn merge_from(&mut self, other: ProteinSet) {
        self.malformed += other.malformed;
        self.lossy_lines += other.lossy_lines;
        self.empty_records += other.empty_records;
        self.duplicates.extend(other.duplicates);
        for seq in other.seqs {
            self.push(seq);
        }
    }

    /// Read one stream, appending to this set.
    pub fn read_one<R: BufRead>(&mut self, mut reader: R, path: &Path) -> Result<()> {
        let mut line = String::new();
        let mut scratch: Vec<u8> = Vec::new();
        let mut current: Option<Sequence> = None;
        loop {
            let (n, lossy) = util::read_line_lossy(&mut reader, &mut line, &mut scratch)?;
            self.lossy_lines += u64::from(lossy);
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if trimmed.is_empty() {
                continue;
            }
            if let Some(header) = trimmed.strip_prefix('>') {
                if let Some(seq) = current.take() {
                    self.push(seq);
                }
                let id = header
                    .find(char::is_whitespace)
                    .map_or(header, |i| &header[..i]);
                if id.is_empty() {
                    self.malformed += 1;
                    current = None;
                    continue;
                }
                current = Some(Sequence {
                    id: id.to_string(),
                    // Not stored: nothing reads it, and at ten million proteins it is ~0.5 GB.
                    description: None,
                    residues: Vec::new(),
                    len: 0,
                });
                continue;
            }
            let Some(seq) = current.as_mut() else {
                // Sequence bytes before any header: not FASTA, or a truncated record.
                self.malformed += 1;
                continue;
            };
            for b in trimmed.bytes() {
                match b {
                    b' ' | b'\t' => continue,
                    b if is_residue(b) => {
                        seq.len += 1;
                        if self.keep_residues {
                            seq.residues.push(b.to_ascii_uppercase());
                        }
                    }
                    b if is_gap(b) => continue, // some exports keep gaps in unaligned files
                    _ => {
                        // X, Z, O and stray symbols: keep the residue-like ones, count the rest.
                        if b.is_ascii_punctuation() || b.is_ascii_digit() {
                            self.malformed += 1;
                        } else {
                            seq.len += 1;
                            if self.keep_residues {
                                seq.residues.push(b.to_ascii_uppercase());
                            }
                        }
                    }
                }
            }
        }
        if let Some(seq) = current.take() {
            self.push(seq);
        }
        let _ = path;
        Ok(())
    }

    fn push(&mut self, seq: Sequence) {
        if seq.len == 0 {
            self.empty_records += 1;
            return;
        }
        self.total_residues += seq.len as u64;
        let idx = self.seqs.len();
        match self.exact.get(&seq.id) {
            Some(&first) => {
                self.duplicates
                    .push((seq.id.clone(), self.seqs[first].id.clone()));
                return;
            }
            None => {
                self.exact.insert(seq.id.clone(), idx);
            }
        }
        // A sequence added after the stripped index was built invalidates it.
        if self.stripped.get().is_some() {
            self.stripped = Default::default();
        }
        self.seqs.push(seq);
    }

    /// Version-stripped id → index (`None` when several sequences share the stripped key).
    fn stripped_map(&self) -> &FxHashMap<String, Option<usize>> {
        self.stripped.get_or_init(|| {
            let mut m: FxHashMap<String, Option<usize>> = FxHashMap::default();
            for (idx, seq) in self.seqs.iter().enumerate() {
                m.entry(util::strip_version(&seq.id).to_string())
                    .and_modify(|v| *v = None)
                    .or_insert(Some(idx));
            }
            m
        })
    }

    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    /// Resolve an identifier, exact first and then version-insensitive.
    ///
    /// Both sides are normalised through the same idempotent [`util::strip_version`], so an
    /// unversioned annotation id finds a versioned sequence and vice versa.
    pub fn lookup(&self, id: &str) -> Option<usize> {
        if let Some(&i) = self.exact.get(id) {
            return Some(i);
        }
        if self.strip_versions {
            return self
                .stripped_map()
                .get(util::strip_version(id))
                .copied()
                .flatten();
        }
        None
    }

    /// How `id` was resolved: exactly, through version stripping, or not at all.
    pub fn resolve(&self, id: &str) -> Resolution {
        if self.exact.contains_key(id) {
            return Resolution::Exact;
        }
        if self.strip_versions {
            match self.stripped_map().get(util::strip_version(id)) {
                Some(Some(_)) => return Resolution::VersionInsensitive,
                Some(None) => return Resolution::Ambiguous,
                None => {}
            }
        }
        Resolution::Missing
    }

    /// `true` when `id` was stored exactly (not through a version-insensitive fallback).
    pub fn has_exact(&self, id: &str) -> bool {
        self.exact.contains_key(id)
    }

    pub fn get(&self, idx: usize) -> Option<&Sequence> {
        self.seqs.get(idx)
    }

    pub fn residues(&self, idx: usize) -> Option<&[u8]> {
        self.seqs.get(idx).map(|s| s.residues.as_slice())
    }

    pub fn id(&self, idx: usize) -> Option<&str> {
        self.seqs.get(idx).map(|s| s.id.as_str())
    }

    /// Ids in file order.
    pub fn ids(&self) -> impl Iterator<Item = &str> {
        self.seqs.iter().map(|s| s.id.as_str())
    }

    /// Warn/error-level lines describing join quality.
    pub fn report(&self) -> Vec<(Level, String)> {
        let mut out = Vec::new();
        if self.malformed > 0 {
            out.push((
                Level::Warn,
                format!("{} malformed sequence line(s) were skipped", self.malformed),
            ));
        }
        if self.empty_records > 0 {
            out.push((
                Level::Warn,
                format!("{} record(s) had no residues", self.empty_records),
            ));
        }
        if !self.duplicates.is_empty() {
            out.push((
                Level::Warn,
                format!(
                    "{} duplicate sequence id(s), first occurrence kept (e.g. {})",
                    self.duplicates.len(),
                    self.duplicates[0].0
                ),
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn ambiguous_version_stripped_ids_are_refused() {
        let mut ps = ProteinSet::new();
        ps.read_one(
            std::io::Cursor::new(b">T.1\nMKVL\n>T.2\nMMMM\n>U.1\nMKKK\n"),
            Path::new("m"),
        )
        .unwrap();
        assert_eq!(ps.lookup("T.1"), Some(0), "exact ids still resolve");
        assert_eq!(
            ps.lookup("T.9"),
            None,
            "T.9 could be either isoform: refuse"
        );
        assert_eq!(ps.resolve("T.9"), Resolution::Ambiguous);
        assert_eq!(
            ps.lookup("U.7"),
            Some(2),
            "a unique stripped key still resolves"
        );
        assert_eq!(ps.resolve("U.7"), Resolution::VersionInsensitive);
    }

    #[test]
    fn invalid_utf8_is_read_lossily_and_counted() {
        let mut ps = ProteinSet::new();
        let mut data = b">P1 caf".to_vec();
        data.push(0xe9);
        data.extend_from_slice(b"\nMKVL\n");
        ps.read_one(std::io::Cursor::new(data), Path::new("m"))
            .unwrap();
        assert_eq!(ps.len(), 1);
        assert_eq!(ps.lossy_lines, 1);
    }

    use super::*;
    use std::io::Cursor;

    fn set(data: &[u8]) -> ProteinSet {
        let mut s = ProteinSet::new();
        s.read_one(Cursor::new(data), Path::new("mem")).unwrap();
        s
    }

    #[test]
    fn reads_wrapped_records_and_headers() {
        let s = set(b">P1 first protein\nMKVL\nAAATA\n>P2\nSEQUENCE\n");
        assert_eq!(s.len(), 2);
        assert_eq!(s.id(0), Some("P1"));
        assert_eq!(s.seqs[0].description, None, "descriptions are not kept");
        assert_eq!(s.residues(0), Some(&b"MKVLAAATA"[..]));
        assert_eq!(s.total_residues, 9 + 8, "SEQUENCE is 8 residues");
        assert_eq!(s.malformed, 0);
    }

    #[test]
    fn header_parsing_edge_cases() {
        let s = set(b">P1  spaced  description here\nAAA\n>\nBBB\n>P3 *XU*\nMMM*X\n");
        assert_eq!(s.ids().collect::<Vec<_>>(), vec!["P1", "P3"]);
        assert_eq!(
            s.malformed, 2,
            "the id-less header and its orphaned sequence line are both bad"
        );
        assert_eq!(
            s.seqs[1].residues, b"MMM*X",
            "stop and ambiguous codes are residues"
        );
    }

    #[test]
    fn empty_and_gap_only_records() {
        let s = set(b">A\n\n>B\n---\n>C\nACDE\n");
        assert_eq!(s.len(), 1, "records without residues are dropped");
        assert_eq!(s.empty_records, 2);
    }

    #[test]
    fn duplicate_ids_keep_the_first() {
        let s = set(b">A\nKK\n>A\nLL\n>B\nMM\n");
        assert_eq!(s.len(), 2);
        assert_eq!(s.residues(s.lookup("A").unwrap()), Some(&b"KK"[..]));
        assert_eq!(s.duplicates.len(), 1);
    }

    #[test]
    fn version_insensitive_lookup() {
        let s = set(b">ENSP1.2\nAAAA\n>P_plain\nCC\n");
        assert_eq!(s.lookup("ENSP1.2"), Some(0), "exact");
        assert_eq!(s.lookup("ENSP1"), Some(0), "annotation carries no version");
        assert_eq!(s.lookup("NOPE"), None);
        assert!(s.has_exact("ENSP1.2"));
        assert!(!s.has_exact("ENSP1"));
    }

    #[test]
    fn unversioned_lookup_resolves_against_versioned_id() {
        let s = set(b">NP_001.3\nAAAA\n");
        assert_eq!(s.lookup("NP_001"), Some(0));
    }

    #[test]
    fn reads_gzip_named_paths_via_extension() {
        // Round-trip through flate2 to make sure the extension dispatch works.
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("ip_fasta_gz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("p.faa.gz");
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(b">G1\nMKV\n").unwrap();
        let gz = enc.finish().unwrap();
        std::fs::write(&path, gz).unwrap();
        let s = ProteinSet::read([path]).unwrap();
        assert_eq!(s.residues(0), Some(&b"MKV"[..]));
    }

    #[test]
    fn empty_input_is_an_error() {
        let dir = std::env::temp_dir().join(format!("ip_fasta_empty_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("none.faa");
        std::fs::write(&path, ">A\n\n").unwrap();
        let err = ProteinSet::read([path]).unwrap_err();
        assert!(matches!(err, Error::EmptyInput { .. }), "{err}");
    }

    fn write_fasta(dir: &std::path::Path, name: &str, body: &str) -> PathBuf {
        std::fs::create_dir_all(dir).unwrap();
        let p = dir.join(name);
        std::fs::write(&p, body).unwrap();
        p
    }

    #[test]
    fn lengths_only_and_parallel_reading_match_a_sequential_full_read() {
        let dir = std::env::temp_dir().join(format!("ip_prot_par_{}", std::process::id()));
        let files = vec![
            write_fasta(
                &dir,
                "a.faa",
                ">A1 first\nMKVL\nGA\n>A2\nMM\n>DUP\nAAA\n>E\n",
            ),
            write_fasta(&dir, "b.faa", ">B1\nWWW\n>DUP\nCCCC\n>B2.3\nHHH\n"),
            write_fasta(&dir, "c.faa", ">C1\nLLLLL\n>B2.4\nGG\n"),
        ];
        // Reference: the old sequential algorithm, via read_one into one set.
        let mut reference = ProteinSet::new();
        for f in &files {
            reference
                .read_one(util::open_buffered(f, 1 << 12).unwrap(), f)
                .unwrap();
        }
        let full = ProteinSet::read(&files).unwrap();
        let lens = ProteinSet::read_lengths(&files).unwrap();
        for set in [&full, &lens] {
            assert_eq!(set.len(), reference.len());
            assert_eq!(set.total_residues, reference.total_residues);
            assert_eq!(set.duplicates, reference.duplicates);
            assert_eq!(
                (set.empty_records, set.malformed),
                (reference.empty_records, reference.malformed)
            );
            for i in 0..reference.len() {
                assert_eq!(set.id(i), reference.id(i), "file order is kept");
                assert_eq!(set.get(i).unwrap().len, reference.get(i).unwrap().len);
            }
            for id in ["A1", "DUP", "B2.9", "C1", "nope"] {
                assert_eq!(set.lookup(id), reference.lookup(id), "{id}");
            }
        }
        assert!(full.keeps_residues() && !lens.keeps_residues());
        assert_eq!(full.residues(0).unwrap(), b"MKVLGA");
        assert!(
            lens.residues(0).is_none_or(|r| r.is_empty()),
            "lengths-only keeps no residues"
        );
        assert_eq!(lens.get(0).unwrap().len, 6);
    }
}
