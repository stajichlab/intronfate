//! Small IO/text helpers shared by the readers.

use std::borrow::Cow;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

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

/// Open a plain or gzip-compressed text file (`.gz`, `.bgz`) as a buffered byte reader.
///
/// Callers work in bytes and validate UTF-8 per record, which avoids a whole-file check on
/// multi-gigabyte annotations.
pub fn open_buffered(path: &Path, buf_capacity: usize) -> Result<Box<dyn BufRead>> {
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
