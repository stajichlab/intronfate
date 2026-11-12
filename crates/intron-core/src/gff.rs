//! Line-oriented GFF3 / GTF2 reader.
//!
//! The format is simple enough that a hand-rolled parser is both faster and more forgiving than a
//! generic one, and forgiveness matters here: the goal is to project introns out of whatever an
//! annotation pipeline emitted (MAKER, Augustus/Braker, NCBI, Ensembl, protein2genome) without
//! rejecting a 3 GB file over one stray line.
//!
//! Records borrow from the reader's line buffer; attribute keys/values are kept in the reader as
//! byte spans, which keeps the hot loop allocation-free and the code free of `unsafe`. Callers must
//! copy anything they retain past the next [`GffReader::next_record`] call.
//!
//! ```
//! use intron_core::gff::GffReader;
//! let gff = b"I\ttest\tCDS\t100\t114\t.\t+\t0\tID=c1;Parent=t1;protein_id=P1\n\
//!            I\ttest\tCDS\t200\t221\t.\t+\t1\tID=c2;Parent=t1;protein_id=P1\n";
//! let mut reader = GffReader::new(&gff[..], "mem".into());
//! let mut n = 0;
//! while let Some(rec) = reader.next_record().unwrap() {
//!     assert_eq!(rec.attr("Parent"), Some("t1"));
//!     assert_eq!(rec.attr("protein_id"), Some("P1"));
//!     n += 1;
//! }
//! assert_eq!(n, 2);
//! ```

use std::io::BufRead;
use std::ops::Range;
use std::path::PathBuf;

use crate::error::{Error, Result};
use crate::model::Strand;
use crate::util;

/// Byte range of an attribute's key and value within the current line.
#[derive(Clone, Copy, Debug, Default)]
struct AttrSpan {
    key: (u32, u32),
    value: (u32, u32),
}

impl AttrSpan {
    fn slice<'a>(line: &'a str, r: &(u32, u32)) -> &'a str {
        &line[r.0 as usize..r.1 as usize]
    }
}

/// One parsed attribute, borrowed from the line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Attribute<'a> {
    pub key: &'a str,
    /// Raw value, still percent-encoded; call [`crate::util::unescape_gff`] when needed.
    pub value: &'a str,
}

/// One feature line.
#[derive(Clone, Copy, Debug)]
pub struct GffRecord<'a> {
    pub seqid: &'a str,
    pub source: &'a str,
    pub feature_type: &'a str,
    /// 1-based inclusive.
    pub start: u64,
    /// 1-based inclusive.
    pub end: u64,
    pub score: Option<f64>,
    pub strand: Strand,
    /// GFF3 `phase` / GTF `frame` for CDS features: the number of bases to remove from the start of
    /// this feature to reach the first base of the next complete codon, i.e.
    /// `(3 - coding_bases_before % 3) % 3`. Not the count of bases already emitted.
    pub frame: Option<u8>,
    /// Line the record was read from, so a diagnostic can quote it without re-borrowing the reader.
    pub line_no: usize,
    /// Path the record came from.
    pub source_path: &'a str,
    line: &'a str,
    attrs: &'a [AttrSpan],
}

impl<'a> GffRecord<'a> {
    /// Attribute key/value pairs in file order.
    pub fn attributes(&self) -> impl Iterator<Item = Attribute<'a>> + 'a {
        // Copy the two `&'a` fields out, so the closure never captures `&self`.
        let (line, attrs) = (self.line, self.attrs);
        attrs.iter().map(move |s| Attribute {
            key: AttrSpan::slice(line, &s.key),
            value: AttrSpan::slice(line, &s.value),
        })
    }

    /// Value for `key` (case-insensitive), ignoring empty and `.` placeholders.
    pub fn attr(&self, key: &str) -> Option<&'a str> {
        self.attributes().find_map(|a| {
            (a.key.eq_ignore_ascii_case(key) && !a.value.is_empty() && a.value != ".")
                .then_some(a.value)
        })
    }

    /// Value for the first key in `keys` that is present, respecting that order. Lets callers
    /// express "`transcript_id`, else `Parent`" fallbacks across differing pipeline conventions.
    pub fn attr_first_of(&self, keys: &[String]) -> Option<&'a str> {
        keys.iter().find_map(|k| self.attr(k))
    }

    /// All `Parent` entries (GFF3 permits comma-separated multi-parents).
    pub fn parents(&self) -> impl Iterator<Item = &'a str> {
        self.attr("Parent")
            .into_iter()
            .flat_map(|p| p.split(','))
            .map(str::trim)
            .filter(|p| !p.is_empty())
    }

    pub fn length(&self) -> u64 {
        self.end.saturating_sub(self.start).saturating_add(1)
    }

    /// Compact one-line description of the record, for diagnostics.
    pub fn summary(&self) -> String {
        format!(
            "{} {} {}-{} {}",
            self.feature_type, self.seqid, self.start, self.end, self.strand
        )
    }
}

/// Streaming reader over GFF3/GTF text.
#[derive(Debug)]
pub struct GffReader<R> {
    reader: R,
    path: PathBuf,
    /// Owned copy for `GffRecord::source_path`, which must borrow for the record's lifetime.
    path_string: String,
    line: String,
    /// Reused attribute storage: byte spans only, so `Self` carries no lifetime.
    spans: Vec<AttrSpan>,
    line_no: usize,
    /// Lines that could not be parsed. Counted, and fatal only when [`strict`](Self::strict).
    pub malformed: u64,
    pub records: u64,
    pub directives: u64,
    /// Lines that contained invalid UTF-8 and were read lossily.
    pub lossy: u64,
    scratch: Vec<u8>,
    pub strict: bool,
}

/// Parsed fixed columns, expressed as ranges plus scalars so they can cross the borrow boundary
/// from the mutable parse phase to the immutable record phase.
#[derive(Debug)]
struct Cols {
    seqid: Range<usize>,
    source: Range<usize>,
    feature_type: Range<usize>,
    start: u64,
    end: u64,
    score: Option<f64>,
    strand: Strand,
    frame: Option<u8>,
}

enum Step {
    Eof,
    Skip,
    Parsed(Cols),
    Bad(String),
}

impl<R: BufRead> GffReader<R> {
    pub fn new(reader: R, path: PathBuf) -> Self {
        let path_string = path.display().to_string();
        GffReader {
            reader,
            path,
            path_string,
            line: String::with_capacity(1024),
            spans: Vec::with_capacity(16),
            line_no: 0,
            malformed: 0,
            records: 0,
            directives: 0,
            lossy: 0,
            scratch: Vec::with_capacity(1024),
            strict: false,
        }
    }

    pub fn path(&self) -> &PathBuf {
        &self.path
    }

    pub fn line_no(&self) -> usize {
        self.line_no
    }

    /// Next feature record, or `None` at EOF. Directives (`##…`) and comments are skipped.
    pub fn next_record(&mut self) -> Result<Option<GffRecord<'_>>> {
        loop {
            // Phase 1: read + parse, holding mutable access to `self`. Nothing borrowed survives.
            let step = {
                self.spans.clear();
                self.line_no += 1;
                let (nread, lossy) =
                    util::read_line_lossy(&mut self.reader, &mut self.line, &mut self.scratch)?;
                self.lossy += u64::from(lossy);
                if nread == 0 {
                    Step::Eof
                } else {
                    let raw = self
                        .line
                        .trim_end_matches(['\n', '\r'])
                        .trim_start_matches(' ');
                    if raw.is_empty() {
                        Step::Skip
                    } else if raw.as_bytes()[0] == b'#' {
                        // Counted here rather than in phase 2, so the two skip reasons stay distinct.
                        self.directives += 1;
                        Step::Skip
                    } else {
                        match parse_columns(raw, &mut self.spans) {
                            Ok(cols) => Step::Parsed(cols),
                            Err(message) => Step::Bad(message),
                        }
                    }
                }
            };

            // Phase 2: re-borrow immutably; only this path hands references to the caller.
            match step {
                Step::Eof => return Ok(None),
                Step::Skip => {}
                Step::Parsed(c) => {
                    self.records += 1;
                    let line: &str = &self.line;
                    return Ok(Some(GffRecord {
                        seqid: &line[c.seqid],
                        source: &line[c.source],
                        feature_type: &line[c.feature_type],
                        start: c.start,
                        end: c.end,
                        score: c.score,
                        strand: c.strand,
                        frame: c.frame,
                        line_no: self.line_no,
                        source_path: &self.path_string,
                        line,
                        attrs: &self.spans,
                    }));
                }
                Step::Bad(message) => {
                    self.malformed += 1;
                    if self.strict {
                        let raw = self.line.trim_end();
                        return Err(Error::Parse {
                            message: format!("{message} in {:?}", util::ellipsis(raw, 200)),
                            location: Some(format!("{}:{}", self.path.display(), self.line_no)),
                        });
                    }
                }
            }
        }
    }
}

impl GffReader<Box<dyn BufRead>> {
    /// Open a plain or gzip-compressed file for streaming.
    pub fn open(path: impl Into<PathBuf>) -> Result<GffReader<Box<dyn BufRead>>> {
        let path = path.into();
        let br = util::open_buffered(&path, 1 << 18)?;
        Ok(GffReader::new(br, path))
    }
}

/// Byte offsets of the tab-separated columns within `line`.
fn column_bounds(line: &str) -> Vec<(usize, usize)> {
    let mut bounds = Vec::with_capacity(9);
    let mut prev = 0usize;
    for (idx, c) in line.char_indices() {
        if c == '\t' {
            bounds.push((prev, idx));
            prev = idx + 1;
        }
    }
    bounds.push((prev, line.len()));
    bounds
}

fn parse_columns(line: &str, attrs: &mut Vec<AttrSpan>) -> std::result::Result<Cols, String> {
    let bounds = column_bounds(line);
    if bounds.len() < 8 {
        return Err(format!(
            "expected 9 tab-separated columns, found {}",
            bounds.len()
        ));
    }
    let text = |i: usize| &line[bounds[i].0..bounds[i].1];

    if text(0).is_empty() {
        return Err("empty seqid".into());
    }
    if text(2).is_empty() {
        return Err("empty feature_type".into());
    }
    let start = parse_pos(text(3)).ok_or_else(|| format!("bad start {:?}", text(3)))?;
    let end = parse_pos(text(4)).ok_or_else(|| format!("bad end {:?}", text(4)))?;
    if end < start {
        return Err(format!("end {end} < start {start}"));
    }
    let strand = match text(6).as_bytes().first() {
        Some(b) => Strand::parse(*b).ok_or_else(|| format!("bad strand {:?}", text(6)))?,
        None => Strand::Unknown,
    };
    let score = match text(5) {
        "." | "" | "nan" | "None" | "NA" => None,
        s => s.parse::<f64>().ok(),
    };
    let frame = match text(7) {
        "." | "" => None,
        s => s
            .parse::<i64>()
            .ok()
            .and_then(|v| u8::try_from(v).ok())
            .filter(|v| *v <= 2),
    };

    if bounds.len() >= 9 {
        parse_attributes(text(8), bounds[8].0, attrs);
    }
    Ok(Cols {
        seqid: bounds[0].0..bounds[0].1,
        source: bounds[1].0..bounds[1].1,
        feature_type: bounds[2].0..bounds[2].1,
        start,
        end,
        score,
        strand,
        frame,
    })
}

fn parse_pos(s: &str) -> Option<u64> {
    let s = s.strip_prefix('+').unwrap_or(s);
    if s.is_empty() || s == "." {
        return None;
    }
    s.parse::<u64>().ok()
}

/// Parse GFF3 (`key=value;`) and GTF (`key "value";`) attribute syntaxes into spans over `line`.
///
/// `base` is the byte offset of `s` within the line, so emitted spans index the whole line.
fn parse_attributes(s: &str, base: usize, out: &mut Vec<AttrSpan>) {
    let mut abs = base;
    for chunk in s.split(';') {
        let lead = chunk.len() - chunk.trim_start().len();
        let trimmed = chunk.trim();
        let start = abs + lead;
        abs += chunk.len() + 1; // +1 for the consumed separator
        if trimmed.is_empty() {
            continue;
        }
        if let Some(span) = parse_attr_chunk(trimmed, start) {
            out.push(span);
        }
    }
}

/// Parse one whitespace-trimmed attribute chunk that begins at absolute byte offset `at`.
fn parse_attr_chunk(chunk: &str, at: usize) -> Option<AttrSpan> {
    let span_of = |off: usize, text: &str| (off as u32, (off + text.len()) as u32);

    if let Some(eq) = chunk.find('=') {
        let before = &chunk[..eq];
        let after = chunk[eq + 1..].trim_start();
        // GFF3 `key=value`, unless the value is quoted, which means GTF written as `key = "value"`.
        // GTF writes `key "value"`; a `=` inside the quoted value then follows whitespace or a quote
        // in what we took for the key. A GFF3 key is a single token, so anything else is GTF.
        let _ = after;
        let gtf_like = before.trim().contains(char::is_whitespace) || before.contains('"');
        if !gtf_like {
            let key = before.trim();
            let kl = before.len() - before.trim_start().len();
            let mut val = chunk[eq + 1..].trim();
            let mut vl = chunk[eq + 1..].len() - chunk[eq + 1..].trim_start().len();
            // GFF3 values are not meant to be quoted, but some writers do it: drop one pair.
            if let Some(inner) = val.strip_prefix('"').and_then(|v| v.strip_suffix('"')) {
                val = inner;
                vl += 1;
            }
            if key.is_empty() || val.is_empty() || val == "." {
                return None;
            }
            return Some(AttrSpan {
                key: span_of(at + kl, key),
                value: span_of(at + eq + 1 + vl, val),
            });
        }
    }

    // GTF: `key "value"` or bare `key value`.
    let ws = chunk.find(char::is_whitespace)?;
    let key = &chunk[..ws];
    let rest = chunk[ws..]
        .trim_start()
        .trim_start_matches('=')
        .trim_start();
    let offset = chunk.len() - rest.len();
    let val_raw = rest.trim_end();
    let quoted = val_raw.strip_prefix('"').and_then(|r| r.strip_suffix('"'));
    let (val, skip) = match quoted {
        Some(v) => (v, 1usize),
        None => (val_raw, 0usize),
    };
    if key.is_empty() || val.is_empty() || val == "." {
        return None;
    }
    Some(AttrSpan {
        key: span_of(at, key),
        value: span_of(at + offset + skip, val),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Snapshot a record's borrowed data into owned strings for assertions.
    #[derive(Clone, Debug, PartialEq)]
    struct Owned {
        seqid: String,
        feature_type: String,
        start: u64,
        end: u64,
        strand: Strand,
        frame: Option<u8>,
        score: Option<f64>,
        attrs: Vec<(String, String)>,
    }

    fn take(rec: &GffRecord<'_>) -> Owned {
        Owned {
            seqid: rec.seqid.to_string(),
            feature_type: rec.feature_type.to_string(),
            start: rec.start,
            end: rec.end,
            strand: rec.strand,
            frame: rec.frame,
            score: rec.score,
            attrs: rec
                .attributes()
                .map(|a| (a.key.to_string(), a.value.to_string()))
                .collect(),
        }
    }

    fn all(data: &[u8]) -> (Vec<Owned>, GffReader<Cursor<&[u8]>>) {
        let mut r = GffReader::new(Cursor::new(data), "mem".into());
        let mut out = Vec::new();
        while let Some(rec) = r.next_record().unwrap() {
            out.push(take(&rec));
        }
        (out, r)
    }

    #[test]
    fn parses_gff3() {
        let (recs, r) = all(
            b"chr1\tmaker\tCDS\t100\t114\t.\t+\t0\tID=c.1;Parent=t.1;protein_id=P.1;gene=GENE1\n",
        );
        assert_eq!(r.malformed, 0);
        assert_eq!(r.records, 1);
        let c = &recs[0];
        assert_eq!(c.seqid, "chr1");
        assert_eq!(c.feature_type, "CDS");
        assert_eq!((c.start, c.end, c.frame), (100, 114, Some(0)));
        assert_eq!(c.strand, Strand::Plus);
        assert_eq!(c.score, None);
        assert_eq!(
            c.attrs,
            vec![
                ("ID".to_string(), "c.1".to_string()),
                ("Parent".to_string(), "t.1".to_string()),
                ("protein_id".to_string(), "P.1".to_string()),
                ("gene".to_string(), "GENE1".to_string())
            ]
        );
    }

    #[test]
    fn accessors() {
        let data = b"s1\texonerate\tCDS\t5\t9\t3.5\t-\t2\tID=a;Parent=t1,t2;note=hi\n";
        let mut r = GffReader::new(&data[..], "m".into());
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(
            rec.attr("PARENT"),
            Some("t1,t2"),
            "keys are case-insensitive"
        );
        assert_eq!(rec.parents().collect::<Vec<_>>(), vec!["t1", "t2"]);
        assert_eq!(
            rec.attr_first_of(&["gene_id".into(), "Parent".into()]),
            Some("t1,t2")
        );
        assert_eq!(rec.attr("missing"), None);
        assert_eq!(
            (rec.score, rec.strand, rec.frame),
            (Some(3.5), Strand::Minus, Some(2))
        );
        assert_eq!(rec.length(), 5);
        assert_eq!(rec.summary(), "CDS s1 5-9 -");
    }

    #[test]
    fn quoted_gtf_values_with_equals() {
        let (recs, r) =
            all(b"chr1\thavana\tCDS\t9\t26\t.\t-\t0\tgene_id \"G\"; product \"DNA=binding\"; tag \"CCDS\";\n");
        assert_eq!(r.malformed, 0);
        let m = &recs[0].attrs;
        assert!(
            m.contains(&("gene_id".to_string(), "G".to_string())),
            "{m:?}"
        );
        assert!(
            m.contains(&("product".to_string(), "DNA=binding".to_string())),
            "quoted value keeps '=': {m:?}"
        );
        assert!(
            m.contains(&("tag".to_string(), "CCDS".to_string())),
            "{m:?}"
        );
    }

    #[test]
    fn unquoted_gtf_values_and_spaced_equals() {
        let (recs, r) =
            all(b"chr1\tx\texon\t9\t26\t.\t+\t.\tgene_id G1; product = \"kinase\";note x\n");
        assert_eq!(r.malformed, 0);
        let m = &recs[0].attrs;
        assert!(
            m.contains(&("gene_id".to_string(), "G1".to_string())),
            "{m:?}"
        );
        assert!(
            m.contains(&("product".to_string(), "kinase".to_string())),
            "{m:?}"
        );
        assert!(m.contains(&("note".to_string(), "x".to_string())), "{m:?}");
    }

    #[test]
    fn skips_directives_and_blank_lines() {
        let data = b"##gff-version 3\n\n#!sequence-region\tchr1\n#comment\nchr1\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n";
        let (recs, r) = all(data);
        assert_eq!(recs.len(), 1);
        assert_eq!(r.directives, 3);
        assert_eq!(r.malformed, 0);
    }

    #[test]
    fn malformed_lines_counted_or_fatal() {
        let data =
            b"chr1\tx\tCDS\tnotanumber\t9\t.\t+\t.\tID=a\nchr1\tx\tCDS\t1\t9\t.\t+\t.\tID=b\n";
        let (recs, r) = all(data);
        assert_eq!(recs.len(), 1, "bad line skipped, good line kept");
        assert_eq!(r.malformed, 1);

        let mut strict = GffReader::new(&data[..], "mem".into());
        strict.strict = true;
        let err = strict.next_record().unwrap_err();
        assert!(
            err.to_string().contains("mem:1"),
            "error names file and line: {err}"
        );
    }

    #[test]
    fn rejects_reversed_coords_and_bad_strand() {
        let (_, r) =
            all(b"chr1\tx\tCDS\t50\t9\t.\t+\t.\tID=a\nchr1\tx\tCDS\t9\t50\t.\tq\t.\tID=b\n");
        assert_eq!(r.malformed, 2);
    }

    #[test]
    fn tolerates_a_missing_attribute_column() {
        // Eight columns is legal input (the attribute column is simply absent); fewer is not.
        let (recs, r) = all(b"chr1\tx\tgene\t1\t9\t.\t+\t.\n");
        assert_eq!(r.malformed, 0);
        assert_eq!(recs[0].attrs, vec![]);
        assert_eq!(recs[0].frame, None);
        let (_, r) = all(b"chr1\tx\tgene\t1\t9\t.\t+\n");
        assert_eq!(r.malformed, 1, "seven columns cannot be parsed as GFF3");
    }

    #[test]
    fn crlf_and_utf8_values() {
        let data = "chr1\tx\tgene\t1\t9\t.\t+\t.\tName=ünïcode;note=a%3Bb\r\n".as_bytes();
        let (recs, r) = all(data);
        assert_eq!(r.malformed, 0);
        assert_eq!(
            recs[0].attrs[0],
            ("Name".to_string(), "ünïcode".to_string())
        );
        assert_eq!(util::unescape_gff(&recs[0].attrs[1].1), "a;b");
    }

    #[test]
    fn empty_attribute_values_are_dropped() {
        let (recs, _) = all(b"chr1\tx\tgene\t1\t9\t.\t+\t.\tID=g1;empty=;dot=.;Parent=t1\n");
        let keys = recs[0]
            .attrs
            .iter()
            .map(|(k, _)| k.as_str())
            .collect::<Vec<_>>();
        assert_eq!(
            keys,
            vec!["ID", "Parent"],
            "empty and '.' values carry no information"
        );
    }

    #[test]
    fn many_attributes_reuse_storage() {
        let mut line = String::from("chr1\tmk\tCDS\t1\t9\t.\t+\t0\tID=x");
        for i in 0..200 {
            line.push_str(&format!(";k{i}=v{i}"));
        }
        let (recs, _) = all(line.as_bytes());
        assert_eq!(recs[0].attrs.len(), 201);
        assert_eq!(recs[0].attrs[200], ("k199".to_string(), "v199".to_string()));
    }

    #[test]
    fn percent_escapes_and_spaces_in_gff3_values() {
        let (recs, _) =
            all(b"c1\tx\tCDS\t1\t9\t.\t+\t0\tproduct=hypothetical%20protein;note=a%3Bb=c\n");
        let m = &recs[0].attrs;
        assert_eq!(m[1], ("note".to_string(), "a%3Bb=c".to_string()));
        assert_eq!(util::unescape_gff(&m[0].1), "hypothetical protein");
    }
}
