//! Small IO/text helpers shared by the readers.

use std::borrow::Cow;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::error::{Error, Result};

fn hex_digit(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

/// Percent-decode GFF3 attribute escaping (`%3B` -> `;`, `%3D` -> `=`, `%20` -> space, …).
///
/// Borrows when the value needs no work, which is the overwhelmingly common case.
pub fn unescape_gff(value: &str) -> Cow<'_, str> {
    if value.is_ascii() && !value.contains('%') {
        return Cow::Borrowed(value);
    }
    if !value.contains('%') {
        return Cow::Borrowed(value);
    }
    let bytes = value.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            // Decode on bytes: slicing the `str` here would panic inside a multibyte character, and
            // `from_str_radix` would accept a sign (`%+f`).
            match (hex_digit(bytes[i + 1]), hex_digit(bytes[i + 2])) {
                (Some(hi), Some(lo)) => {
                    out.push(hi << 4 | lo);
                    i += 3;
                    continue;
                }
                // Not a valid escape: emit the '%' and let the byte copy handle any UTF-8 tail.
                _ => out.push(b'%'),
            }
        } else {
            out.push(bytes[i]);
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned().into()
}

/// `-` names standard input wherever an input file is expected.
pub fn is_stdin(path: &Path) -> bool {
    path.as_os_str() == "-"
}

static STDIN_TAKEN: AtomicBool = AtomicBool::new(false);

/// Claim standard input. It can be read once, so a second claim in the same run is an error that names the cause
/// (two inputs given as `-`, or one input that is read twice).
fn claim(flag: &AtomicBool) -> Result<()> {
    if flag.swap(true, Ordering::SeqCst) {
        return Err(Error::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "standard input (-) can supply only one input per run (it was already used, or two inputs are given as -)",
        )));
    }
    Ok(())
}

fn stdin_reader(
    flag: &AtomicBool,
    input: impl Read + 'static,
    buf_capacity: usize,
) -> Result<Box<dyn BufRead>> {
    claim(flag)?;
    // A stream cannot be reopened to look at its first bytes, so peek through the buffer: gzip is recognised by its magic
    // number, whatever the data is called.
    let mut br = BufReader::with_capacity(buf_capacity.max(1 << 16), input);
    let gz = br.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    Ok(if gz {
        Box::new(BufReader::with_capacity(
            buf_capacity.max(8 * 1024),
            flate2::read::MultiGzDecoder::new(br),
        ))
    } else {
        Box::new(br)
    })
}

/// The whole text of a file, or of standard input when `path` is `-`. Gzip is undone by name (`.gz`, `.bgz`) or by its
/// magic number, as in [`open_buffered`].
pub fn read_text(path: &Path) -> std::io::Result<String> {
    let mut reader = open_buffered(path, 1 << 16).map_err(|e| match e {
        Error::Io(e) => e,
        other => std::io::Error::other(other.to_string()),
    })?;
    let mut s = String::new();
    reader.read_to_string(&mut s)?;
    Ok(s)
}

/// Open a plain or gzip-compressed text file (`.gz`, `.bgz`) as a buffered byte reader. `-` is standard input (gzip is
/// recognised by its first two bytes); it can be used for one input per run.
///
/// Callers work in bytes and validate UTF-8 per record, which avoids a whole-file check on
/// multi-gigabyte annotations.
pub fn open_buffered(path: &Path, buf_capacity: usize) -> Result<Box<dyn BufRead>> {
    if is_stdin(path) {
        return stdin_reader(&STDIN_TAKEN, std::io::stdin(), buf_capacity);
    }
    let file = std::fs::File::open(path).map_err(|e| {
        Error::Io(std::io::Error::new(
            e.kind(),
            format!("could not open {}: {e}", path.display()),
        ))
    })?;
    let reader: Box<dyn Read> = if is_gzip(path) || looks_gzipped(path) {
        Box::new(flate2::read::MultiGzDecoder::new(BufReader::with_capacity(
            1 << 20,
            file,
        )))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::with_capacity(
        buf_capacity.max(8 * 1024),
        reader,
    )))
}

pub fn is_gzip(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|e| e.to_str()).map(|e| e.to_ascii_lowercase()),
        Some(ext) if ext == "gz" || ext == "bgz"
    )
}

/// `true` when the file starts with the gzip magic bytes, regardless of its name.
pub fn looks_gzipped(path: &Path) -> bool {
    let Ok(mut f) = std::fs::File::open(path) else {
        return false;
    };
    let mut magic = [0u8; 2];
    matches!(f.read_exact(&mut magic), Ok(()) if magic == [0x1f, 0x8b])
}

/// Strip a trailing `.N` sequence version, as used by Ensembl/RefSeq/GenPept.
///
/// Only `.` is treated as a version separator on purpose: `-` and `_` appear inside accession
/// bodies (`NP_001234`, `snoRNA-Zic4`), and this function is applied to both sides of a join, so it
/// must also be idempotent — `strip_version(strip_version(x)) == strip_version(x)`.
pub fn strip_version(id: &str) -> &str {
    let mut id = id;
    // Strip every trailing `.<digits>` segment, so the result is a fixed point (`X.1.2` -> `X`):
    // annotation ids are normalised before the join and sequence ids again inside it, and a
    // function that removed only one segment would give the two sides different keys.
    while let Some(cut) = id.rfind('.') {
        let tail = &id[cut + 1..];
        if cut > 0
            && !tail.is_empty()
            && tail.len() <= 12
            && tail.bytes().all(|b| b.is_ascii_digit())
        {
            id = &id[..cut];
        } else {
            break;
        }
    }
    id
}

/// Read one line into `buf` (cleared first), replacing invalid UTF-8 instead of failing.
///
/// Legacy annotation and alignment files routinely carry Latin-1 bytes in free-text fields; one such
/// byte must not abort a multi-gigabyte parse. Returns `(bytes_read, was_lossy)`; `0` bytes means EOF.
pub fn read_line_lossy<R: BufRead + ?Sized>(
    reader: &mut R,
    buf: &mut String,
    scratch: &mut Vec<u8>,
) -> std::io::Result<(usize, bool)> {
    buf.clear();
    scratch.clear();
    let n = reader.read_until(b'\n', scratch)?;
    match std::str::from_utf8(scratch) {
        Ok(s) => {
            buf.push_str(s);
            Ok((n, false))
        }
        Err(_) => {
            buf.push_str(&String::from_utf8_lossy(scratch));
            Ok((n, true))
        }
    }
}

/// Write `bytes` to `path` through a sibling temp file and a rename, so an interrupted run never
/// leaves a truncated table that a later `score` would read as if it were complete.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    prepare_parent(path)?;
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(format!(".tmp{}", std::process::id()));
    let tmp = std::path::PathBuf::from(tmp);
    std::fs::write(&tmp, bytes)?;
    if let Err(e) = std::fs::rename(&tmp, path) {
        let _ = std::fs::remove_file(&tmp);
        return Err(e.into());
    }
    Ok(())
}

/// Like [`write_atomic`], but the content is produced by `f` straight into a buffered file, so a table
/// with hundreds of millions of cells never exists in memory as one `String`.
pub fn write_atomic_with(
    path: &Path,
    f: impl FnOnce(&mut dyn std::io::Write) -> std::io::Result<()>,
) -> Result<()> {
    use std::io::Write;
    prepare_parent(path)?;
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(format!(".tmp{}", std::process::id()));
    let tmp = std::path::PathBuf::from(tmp);
    let outcome = (|| -> std::io::Result<()> {
        let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&tmp)?);
        f(&mut w)?;
        w.flush()
    })();
    match outcome.and_then(|()| std::fs::rename(&tmp, path)) {
        Ok(()) => Ok(()),
        Err(e) => {
            let _ = std::fs::remove_file(&tmp);
            Err(e.into())
        }
    }
}

/// Escape a taxon name for a whitespace- and comma-separated metadata list (`%`, `,` and whitespace).
pub fn encode_list_item(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '%' => out.push_str("%25"),
            ',' => out.push_str("%2C"),
            c if c.is_whitespace() => {
                let mut buf = [0u8; 4];
                for byte in c.encode_utf8(&mut buf).bytes() {
                    out.push_str(&format!("%{byte:02X}"));
                }
            }
            c => out.push(c),
        }
    }
    out
}

/// Inverse of [`encode_list_item`].
pub fn decode_list_item(s: &str) -> String {
    let b = s.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 2 < b.len() {
            if let (Some(h), Some(l)) = (hex_digit(b[i + 1]), hex_digit(b[i + 2])) {
                out.push(h << 4 | l);
                i += 3;
                continue;
            }
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// Ensure the parent directory of an output path exists.
pub fn prepare_parent(path: &Path) -> Result<()> {
    if let Some(dir) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    Ok(())
}

/// Truncate to at most `max` bytes on a char boundary, adding an ellipsis when cut.
pub fn ellipsis(s: &str, max: usize) -> String {
    if s.len() <= max {
        return s.to_string();
    }
    let mut end = max;
    while end > 0 && !s.is_char_boundary(end) {
        end -= 1;
    }
    let mut out = String::with_capacity(end + 3);
    out.push_str(&s[..end]);
    out.push('…');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gz(bytes: &[u8]) -> Vec<u8> {
        use std::io::Write;
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(bytes).unwrap();
        e.finish().unwrap()
    }

    fn all(mut r: Box<dyn BufRead>) -> String {
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        s
    }

    #[test]
    fn stdin_is_named_by_a_dash_only() {
        assert!(is_stdin(Path::new("-")));
        for p in ["", "--", "-x", "./-", "a-", "stdin"] {
            assert!(!is_stdin(Path::new(p)), "{p}");
        }
    }

    #[test]
    fn a_stream_is_read_whole_whether_plain_or_gzipped() {
        let text = "line one\nline two\n".repeat(50_000); // larger than one buffer
        let plain = stdin_reader(
            &AtomicBool::new(false),
            std::io::Cursor::new(text.clone().into_bytes()),
            1 << 12,
        )
        .unwrap();
        assert_eq!(all(plain), text);
        // gzip is recognised by its magic number alone; the name plays no part for a stream
        let zipped = stdin_reader(
            &AtomicBool::new(false),
            std::io::Cursor::new(gz(text.as_bytes())),
            1 << 12,
        )
        .unwrap();
        assert_eq!(all(zipped), text);
        // two members (a concatenated gzip file) are both read
        let mut two = gz(b"first\n");
        two.extend(gz(b"second\n"));
        let r = stdin_reader(&AtomicBool::new(false), std::io::Cursor::new(two), 1 << 12).unwrap();
        assert_eq!(all(r), "first\nsecond\n");
        // empty input is an empty reader, not an error
        let r = stdin_reader(
            &AtomicBool::new(false),
            std::io::Cursor::new(Vec::new()),
            1 << 12,
        )
        .unwrap();
        assert_eq!(all(r), "");
        // a plain file whose first byte is 0x1f but not a gzip header stays plain
        let r = stdin_reader(
            &AtomicBool::new(false),
            std::io::Cursor::new(vec![0x1f, b'a', b'\n']),
            1 << 12,
        )
        .unwrap();
        assert_eq!(all(r).as_bytes(), &[0x1f, b'a', b'\n']);
    }

    #[test]
    fn standard_input_can_be_claimed_once() {
        let flag = AtomicBool::new(false);
        assert!(claim(&flag).is_ok());
        let e = claim(&flag).unwrap_err().to_string();
        assert!(e.contains("only one input"), "{e}");
        // a reader built on a used flag is refused, and the refusal does not consume anything
        let e = stdin_reader(&flag, std::io::Cursor::new(b"x".to_vec()), 1 << 12)
            .err()
            .unwrap()
            .to_string();
        assert!(e.contains("only one input"), "{e}");
    }

    #[test]
    fn read_text_reads_a_file_and_reports_a_missing_one() {
        let p = std::env::temp_dir().join(format!("intronfate_util_text_{}", std::process::id()));
        std::fs::write(&p, "a\tb\n").unwrap();
        assert_eq!(read_text(&p).unwrap(), "a\tb\n");
        // gzip is undone whatever the file is called
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(b"a\tb\n").unwrap();
        std::fs::write(&p, enc.finish().unwrap()).unwrap();
        assert_eq!(read_text(&p).unwrap(), "a\tb\n");
        std::fs::remove_file(&p).unwrap();
        assert!(read_text(&p).is_err());
    }

    #[test]
    fn unescape_basics() {
        assert_eq!(unescape_gff("plain-id"), "plain-id");
        assert_eq!(unescape_gff("a%3Bb"), "a;b");
        assert_eq!(unescape_gff("x%20y"), "x y");
        assert_eq!(unescape_gff("bad%zz"), "bad%zz");
        assert_eq!(unescape_gff("ünïcode %3B ok"), "ünïcode ; ok");
    }

    #[test]
    fn list_items_round_trip_commas_spaces_and_percent() {
        for name in [
            "plain",
            "a,b",
            "s,1",
            "has space",
            "50%",
            "x%2Cy",
            "tab\there",
            "nb\u{a0}sp",
        ] {
            let enc = encode_list_item(name);
            assert!(
                !enc.contains(',') && !enc.contains(char::is_whitespace),
                "{enc}"
            );
            assert_eq!(decode_list_item(&enc), name);
        }
    }

    #[test]
    fn unescape_never_panics_and_decodes_only_real_hex() {
        assert_eq!(unescape_gff("a%3Bb"), "a;b");
        assert_eq!(
            unescape_gff("x%aé"),
            "x%aé",
            "a multibyte tail is not an escape"
        );
        assert_eq!(unescape_gff("%+f"), "%+f", "a sign is not a hex digit");
        assert_eq!(unescape_gff("100%"), "100%");
    }

    #[test]
    fn version_stripping() {
        assert_eq!(strip_version("ENSP00000123456.1"), "ENSP00000123456");
        assert_eq!(strip_version("NP_001234.2"), "NP_001234");
        assert_eq!(
            strip_version("NP_001234"),
            "NP_001234",
            "underscores are part of the accession"
        );
        assert_eq!(strip_version("snoRNA-Zic4"), "snoRNA-Zic4");
        assert_eq!(strip_version("no_version"), "no_version");
        assert_eq!(
            strip_version("P1."),
            "P1.",
            "an empty tail is not a version"
        );
        assert_eq!(strip_version(""), "");
        assert_eq!(
            strip_version("A.1.2"),
            "A",
            "all trailing numeric segments go"
        );
        assert_eq!(
            strip_version("evm.model.scaffold_1.23"),
            "evm.model.scaffold_1"
        );
        assert_eq!(strip_version(".5"), ".5", "never yields an empty id");
        let once = strip_version("ENSP1.4");
        assert_eq!(
            strip_version(once),
            once,
            "must be idempotent: applied to both join sides"
        );
    }

    #[test]
    fn ellipsis_on_boundary() {
        assert_eq!(
            ellipsis("héllo", 4),
            "hél…",
            "cut lands on the é boundary at 4 bytes"
        );
        assert_eq!(ellipsis("héllo", 2), "h…", "never split a multi-byte char");
        assert_eq!(ellipsis("abc", 5), "abc");
    }
}
