//! Expanding input arguments: explicit files, list files, and folders filtered by extension.
//!
//! Thousands of alignments or genomes do not belong on a command line. Every multi-file input accepts the same three forms,
//! combined freely: `--x FILE` (repeatable), `--x-list LIST` (one path per line) and `--x-dir DIR` with `--x-ext EXT`.
//!
//! * A list file holds one path per line; blank lines and lines starting with `#` are skipped, and a relative path is resolved against
//!   the folder of the list file (so a list can travel with its data), not the working directory. A path that does not exist is an
//!   error that names the list file and the line.
//! * A folder contributes the files directly inside it (and, with `depth = 1`, inside its immediate sub-folders) whose name matches
//!   an extension, sorted by name. An extension matches a name that equals it or ends in `.<ext>`, ignoring case and a trailing
//!   `.gz`/`.bgz` (`aln` matches `OG1.aln` and `OG1.aln.gz`; `matrix.tsv` matches `matrix.tsv` and `OG1.matrix.tsv`).
//! * The result keeps the order: explicit files, then each list's lines, then each folder's sorted matches. A path that occurs more
//!   than once is kept once, at its first position. If nothing was requested the result is empty; if something was requested and
//!   nothing was found, that is an error.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};

/// What to expand. `default_exts` apply to folders when the user gave no extension.
#[derive(Debug)]
pub struct Inputs<'a> {
    /// Name used in messages, e.g. `alignment`.
    pub what: &'a str,
    pub files: &'a [PathBuf],
    pub lists: &'a [PathBuf],
    pub dirs: &'a [PathBuf],
    pub exts: &'a [String],
    pub default_exts: &'a [&'a str],
    /// 0: only the folder itself; 1: also its immediate sub-folders.
    pub depth: usize,
    /// Names of sub-folders never entered (e.g. `combined`, the output of `--combine`, which would double-count its own inputs).
    pub skip_dirs: &'a [&'a str],
}

fn bad(msg: String) -> Error {
    Error::Parse {
        message: msg,
        location: None,
    }
}

/// `true` when `name` is `ext` or ends in `.ext`, ignoring case and one trailing `.gz`/`.bgz`.
pub fn matches_ext(name: &str, ext: &str) -> bool {
    let lower = name.to_lowercase();
    let base = lower
        .strip_suffix(".gz")
        .or_else(|| lower.strip_suffix(".bgz"))
        .unwrap_or(&lower);
    let ext = ext.trim_start_matches('.').to_lowercase();
    !ext.is_empty() && (base == ext || base.ends_with(&format!(".{ext}")))
}

fn read_list(list: &Path, out: &mut Vec<PathBuf>, what: &str) -> Result<()> {
    let text = std::fs::read_to_string(list).map_err(|e| {
        bad(format!(
            "could not read {what} list {}: {e}",
            list.display()
        ))
    })?;
    let base = list.parent().unwrap_or_else(|| Path::new(""));
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let p = Path::new(line);
        let p = if p.is_absolute() {
            p.to_path_buf()
        } else {
            base.join(p)
        };
        if !p.is_file() {
            return Err(bad(format!(
                "{}:{}: {what} file {} does not exist",
                list.display(),
                i + 1,
                p.display()
            )));
        }
        out.push(p);
    }
    Ok(())
}

fn scan_dir(
    dir: &Path,
    exts: &[String],
    depth: usize,
    skip: &[&str],
    out: &mut Vec<PathBuf>,
) -> Result<()> {
    let entries = std::fs::read_dir(dir)
        .map_err(|e| bad(format!("could not read folder {}: {e}", dir.display())))?;
    let mut files = Vec::new();
    let mut subs = Vec::new();
    for e in entries {
        let p = e?.path();
        if p.is_dir() {
            let skipped = p
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| skip.contains(&n));
            if !skipped {
                subs.push(p);
            }
        } else if let Some(n) = p.file_name().and_then(|n| n.to_str()) {
            if exts.iter().any(|x| matches_ext(n, x)) {
                files.push(p);
            }
        }
    }
    files.sort();
    out.extend(files);
    if depth > 0 {
        subs.sort();
        for s in subs {
            scan_dir(&s, exts, depth - 1, skip, out)?;
        }
    }
    Ok(())
}

/// Expand `spec` into the ordered, de-duplicated list of files (see the module docs).
pub fn expand(spec: &Inputs) -> Result<Vec<PathBuf>> {
    let mut all: Vec<PathBuf> = spec.files.to_vec();
    for l in spec.lists {
        let before = all.len();
        read_list(l, &mut all, spec.what)?;
        if all.len() == before {
            return Err(bad(format!(
                "{} list {} names no files",
                spec.what,
                l.display()
            )));
        }
    }
    let exts: Vec<String> = if spec.exts.is_empty() {
        spec.default_exts.iter().map(|s| s.to_string()).collect()
    } else {
        spec.exts.to_vec()
    };
    for d in spec.dirs {
        let before = all.len();
        scan_dir(d, &exts, spec.depth, spec.skip_dirs, &mut all)?;
        if all.len() == before {
            return Err(bad(format!(
                "no {} files in {} (looked for extension(s): {}; use the matching --*-ext to change)",
                spec.what,
                d.display(),
                exts.join(", ")
            )));
        }
    }
    let mut seen = HashSet::new();
    all.retain(|p| seen.insert(std::fs::canonicalize(p).unwrap_or_else(|_| p.clone())));
    Ok(all)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("if_inputs_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }
    fn touch(p: &Path) {
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(p, "x").unwrap();
    }
    fn spec<'a>(
        files: &'a [PathBuf],
        lists: &'a [PathBuf],
        dirs: &'a [PathBuf],
        exts: &'a [String],
    ) -> Inputs<'a> {
        Inputs {
            what: "alignment",
            files,
            lists,
            dirs,
            exts,
            default_exts: &["aln", "fasta"],
            depth: 0,
            skip_dirs: &[],
        }
    }

    #[test]
    fn extension_rule() {
        assert!(matches_ext("OG1.aln", "aln"));
        assert!(matches_ext("OG1.ALN.gz", ".aln"));
        assert!(matches_ext("matrix.tsv", "matrix.tsv"));
        assert!(matches_ext("OG1.matrix.tsv", "matrix.tsv"));
        assert!(!matches_ext("OGaln", "aln"), "must be a whole extension");
        assert!(!matches_ext("OG1.aln.bak", "aln"));
        assert!(!matches_ext("x", ""));
    }

    #[test]
    fn nothing_requested_is_empty() {
        assert!(expand(&spec(&[], &[], &[], &[])).unwrap().is_empty());
    }

    #[test]
    fn list_skips_comments_and_resolves_relative_to_the_list() {
        let d = tmp("list");
        touch(&d.join("data/a.aln"));
        touch(&d.join("data/b.aln"));
        std::fs::write(
            d.join("fam.list"),
            "# families\n\ndata/b.aln\n  data/a.aln  \n",
        )
        .unwrap();
        let l = [d.join("fam.list")];
        let got = expand(&spec(&[], &l, &[], &[])).unwrap();
        assert_eq!(got, vec![d.join("data/b.aln"), d.join("data/a.aln")]);
    }

    #[test]
    fn list_errors_name_the_line() {
        let d = tmp("listerr");
        touch(&d.join("a.aln"));
        std::fs::write(d.join("l"), "a.aln\nmissing.aln\n").unwrap();
        let l = [d.join("l")];
        let e = expand(&spec(&[], &l, &[], &[])).unwrap_err().to_string();
        assert!(e.contains(":2:") && e.contains("missing.aln"), "{e}");
        std::fs::write(d.join("empty"), "# nothing\n").unwrap();
        let l2 = [d.join("empty")];
        assert!(expand(&spec(&[], &l2, &[], &[]))
            .unwrap_err()
            .to_string()
            .contains("names no files"));
    }

    #[test]
    fn folder_filters_sorts_and_dedups() {
        let d = tmp("dir");
        for n in ["b.aln", "a.aln", "c.fasta.gz", "notes.txt", "sub/z.aln"] {
            touch(&d.join(n));
        }
        let first = [d.join("a.aln")];
        let dirs = [d.clone()];
        let got = expand(&spec(&first, &[], &dirs, &[])).unwrap();
        assert_eq!(
            got,
            vec![d.join("a.aln"), d.join("b.aln"), d.join("c.fasta.gz")]
        );
        let fasta = ["fasta".to_string()];
        let only = expand(&spec(&[], &[], &dirs, &fasta)).unwrap();
        assert_eq!(only, vec![d.join("c.fasta.gz")]);
        let aln = ["aln".to_string()];
        let mut deep = spec(&[], &[], &dirs, &aln);
        deep.depth = 1;
        assert_eq!(expand(&deep).unwrap().len(), 3, "a, b and sub/z");
    }

    #[test]
    fn skipped_subfolders_are_not_entered() {
        let d = tmp("skip");
        for n in ["fam1/matrix.tsv", "fam2/matrix.tsv", "combined/matrix.tsv"] {
            touch(&d.join(n));
        }
        let dirs = [d.clone()];
        let ext = ["matrix.tsv".to_string()];
        let mut sp = spec(&[], &[], &dirs, &ext);
        sp.depth = 1;
        assert_eq!(expand(&sp).unwrap().len(), 3);
        sp.skip_dirs = &["combined"];
        assert_eq!(
            expand(&sp).unwrap(),
            vec![d.join("fam1/matrix.tsv"), d.join("fam2/matrix.tsv")]
        );
    }

    #[test]
    fn empty_folder_result_is_an_error_that_names_the_extensions() {
        let d = tmp("none");
        touch(&d.join("x.txt"));
        let dirs = [d];
        let e = expand(&spec(&[], &[], &dirs, &[])).unwrap_err().to_string();
        assert!(
            e.contains("aln, fasta") && e.contains("no alignment files"),
            "{e}"
        );
    }
}
