//! Where progress lines go, and how one table reaches standard output (`--stdout TABLE`, `convert --to-malin -o -`).
//!
//! Normal runs print their summary lines to stdout. When stdout carries data they must not, so a process-wide state
//! sends them to stderr (data written next to a directory) or drops them (stdout-only: the lines name a scratch
//! directory the user never asked for).

use std::fmt::Arguments;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU8, Ordering};

use anyhow::{bail, Context, Result};

const NORMAL: u8 = 0;
const STDERR: u8 = 1;
const SILENT: u8 = 2;

static STATE: AtomicU8 = AtomicU8::new(NORMAL);

/// Summary lines now go to stderr (stdout carries data).
pub fn progress_to_stderr() {
    STATE.store(STDERR, Ordering::Relaxed);
}

/// Summary lines are dropped (stdout carries data and no output directory was asked for).
pub fn progress_silent() {
    STATE.store(SILENT, Ordering::Relaxed);
}

/// Backend of `say!`: stdout normally, stderr or nothing once stdout carries data.
pub fn say(args: Arguments) {
    match STATE.load(Ordering::Relaxed) {
        NORMAL => println!("{args}"),
        STDERR => eprintln!("{args}"),
        _ => {}
    }
}

/// Run-level tables, written next to the per-alignment folders.
const RUN_LEVEL: &[&str] = &["introns", "unmatched", "diagnostics", "diagnostics_detail"];
const MATRIX_LEVEL: &[&str] = &[
    "sites",
    "taxa",
    "matrix",
    "presence",
    "pairwise",
    "sharing_matrix",
    "spectrum",
    "counts_sites",
];
const SCORE_LEVEL: &[&str] = &[
    "summary_score",
    "dollo_sites",
    "dollo_branches",
    "mp_sites",
    "dollo_weighted",
    "dollo_weighted_sites",
    "reconstructions",
    "likelihood_branches",
    "likelihood_tips",
    "likelihood_sites",
];

/// Table names `--stdout` accepts for a command (the file stems it can write).
pub fn table_names(command: &str) -> Vec<&'static str> {
    let mut v: Vec<&'static str> = Vec::new();
    match command {
        "project" => {
            v.extend(RUN_LEVEL);
            v.push("summary_project");
        }
        "matrix" => {
            v.extend(RUN_LEVEL);
            v.extend(MATRIX_LEVEL);
        }
        "score" => {
            v.extend(["pairwise", "sharing_matrix", "counts_sites"]);
            v.extend(SCORE_LEVEL);
        }
        "analyze" => {
            v.extend(RUN_LEVEL);
            v.extend(MATRIX_LEVEL);
            v.extend(SCORE_LEVEL);
        }
        _ => {}
    }
    v
}

/// Reject a table name the command never writes, before any work starts.
pub fn check_table_name(command: &str, name: &str) -> Result<()> {
    let names = table_names(command);
    if names.contains(&name) {
        return Ok(());
    }
    bail!(
        "--stdout {name}: `{command}` has no table of that name; choose one of: {}",
        names.join(", ")
    )
}

/// Scratch directory for stdout-only runs; removed on drop.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new() -> Result<Self> {
        let dir = std::env::temp_dir().join(format!("intronfate_stdout_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("could not create scratch directory {}", dir.display()))?;
        Ok(Scratch(dir))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The single `<name>.tsv` under `dir`: directly in it (run-level tables, `project`, `score`) or in exactly one
/// subfolder (the per-alignment folder of a one-alignment run).
pub fn find_table(dir: &Path, name: &str) -> Result<PathBuf> {
    let file = format!("{name}.tsv");
    let direct = dir.join(&file);
    if direct.is_file() {
        return Ok(direct);
    }
    let mut hits = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let p = entry?.path().join(&file);
        if p.is_file() {
            hits.push(p);
        }
    }
    match hits.len() {
        1 => Ok(hits.remove(0)),
        0 => {
            let mut have: Vec<String> = Vec::new();
            for e in std::fs::read_dir(dir)?.flatten() {
                let p = e.path();
                let dirs = if p.is_dir() {
                    std::fs::read_dir(&p)?.flatten().map(|x| x.path()).collect()
                } else {
                    vec![p]
                };
                for f in dirs {
                    if f.extension().is_some_and(|x| x == "tsv") {
                        if let Some(s) = f.file_stem().and_then(|s| s.to_str()) {
                            have.push(s.to_string());
                        }
                    }
                }
            }
            have.sort();
            have.dedup();
            bail!(
                "--stdout {name}: this run did not write that table (it depends on the options given); tables written: {}",
                have.join(", ")
            )
        }
        _ => bail!("--stdout {name}: more than one table of that name was written"),
    }
}

/// Write bytes to stdout. A reader that closes early (`| head`) is a normal stop, not an error.
pub fn emit_bytes(bytes: &[u8]) -> Result<()> {
    let mut out = std::io::stdout().lock();
    match out.write_all(bytes).and_then(|()| out.flush()) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e).context("could not write to standard output"),
    }
}

/// Copy a finished table file to stdout without holding it in memory.
pub fn emit_file(path: &Path) -> Result<()> {
    let mut f =
        std::fs::File::open(path).with_context(|| format!("could not open {}", path.display()))?;
    let mut out = std::io::stdout().lock();
    match std::io::copy(&mut f, &mut out).and_then(|_| out.flush()) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e).context("could not write to standard output"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unknown_names_list_the_valid_ones() {
        let e = check_table_name("project", "nope").unwrap_err().to_string();
        assert!(
            e.contains("introns") && e.contains("summary_project"),
            "{e}"
        );
        assert!(check_table_name("project", "introns").is_ok());
        assert!(check_table_name("score", "introns").is_err());
        assert!(check_table_name("analyze", "likelihood_sites").is_ok());
    }

    #[test]
    fn find_table_looks_in_the_dir_then_in_exactly_one_subfolder() {
        let d = std::env::temp_dir().join(format!("outmode_find_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(d.join("fam")).unwrap();
        std::fs::write(d.join("introns.tsv"), "a\n").unwrap();
        std::fs::write(d.join("fam").join("matrix.tsv"), "b\n").unwrap();
        assert_eq!(find_table(&d, "introns").unwrap(), d.join("introns.tsv"));
        assert_eq!(
            find_table(&d, "matrix").unwrap(),
            d.join("fam").join("matrix.tsv")
        );
        let e = find_table(&d, "mp_sites").unwrap_err().to_string();
        assert!(e.contains("did not write") && e.contains("matrix"), "{e}");
        std::fs::create_dir_all(d.join("fam2")).unwrap();
        std::fs::write(d.join("fam2").join("matrix.tsv"), "c\n").unwrap();
        assert!(find_table(&d, "matrix").is_err());
        let _ = std::fs::remove_dir_all(&d);
    }
}
