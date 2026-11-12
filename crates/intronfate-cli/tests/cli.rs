//! End-to-end tests that run the real `intronfate` binary.
//!
//! The `examples/` directories double as the golden-output suite: each has a `run.sh` and an
//! `expected/` folder, and the first test below re-runs every example and requires byte-identical
//! output. Hand-checked values in those examples are listed in `examples/*/README.md`.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_intronfate");

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("ip_cli_{}_{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().expect("binary runs")
}

fn ok(args: &[&str]) -> String {
    let o = run(args);
    assert!(
        o.status.success(),
        "`{}` failed: {}",
        args.join(" "),
        String::from_utf8_lossy(&o.stderr)
    );
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn copy_dir(from: &Path, to: &Path, skip: &[&str]) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let name = e.file_name();
        if skip.iter().any(|s| name.to_str() == Some(*s)) {
            continue;
        }
        let (src, dst) = (e.path(), to.join(&name));
        if src.is_dir() {
            copy_dir(&src, &dst, skip);
        } else {
            fs::copy(&src, &dst).unwrap();
        }
    }
}

fn files_under(root: &Path) -> Vec<PathBuf> {
    fn walk(dir: &Path, root: &Path, out: &mut Vec<PathBuf>) {
        for e in fs::read_dir(dir).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                walk(&p, root, out);
            } else {
                out.push(p.strip_prefix(root).unwrap().to_path_buf());
            }
        }
    }
    let mut v = Vec::new();
    walk(root, root, &mut v);
    v.sort();
    v
}

fn example(name: &str) -> PathBuf {
    repo().join("examples").join(name)
}

#[test]
fn every_example_reproduces_its_expected_output_exactly() {
    for name in ["01_quickstart", "02_missing_data", "03_malin"] {
        let work = scratch(name);
        copy_dir(&example(name), &work, &["out", "expected"]);
        let st = Command::new("sh")
            .arg("run.sh")
            .arg("out")
            .env("INTRONFATE", BIN)
            .current_dir(&work)
            .output()
            .unwrap();
        assert!(
            st.status.success(),
            "{name}/run.sh failed: {}",
            String::from_utf8_lossy(&st.stderr)
        );
        let expected = example(name).join("expected");
        let got = work.join("out");
        let (want_files, got_files) = (files_under(&expected), files_under(&got));
        assert_eq!(want_files, got_files, "{name}: set of output files differs");
        for f in want_files {
            let a = fs::read_to_string(expected.join(&f)).unwrap();
            let b = fs::read_to_string(got.join(&f)).unwrap();
            assert_eq!(a, b, "{name}: {} differs from expected", f.display());
        }
    }
}

#[test]
fn quickstart_numbers_match_the_hand_calculation() {
    let w = scratch("hand");
    let q = example("01_quickstart");
    let out = w.join("o");
    let sc = ok(&[
        "analyze",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        q.join("family.aln.fasta").to_str().unwrap(),
        "-t",
        q.join("family.nwk").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(sc.contains("dollo 2g/2l, mp 2 steps"), "{sc}");
    let introns = fs::read_to_string(out.join("introns.tsv")).unwrap();
    // A and B: offsets 12 and 21 (+ strand); C: 21 on the minus strand; D: 12.
    let offsets: Vec<(String, String)> = introns
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[10].to_string())
        })
        .collect();
    let want = [
        ("A", "12"),
        ("A", "21"),
        ("B", "12"),
        ("B", "21"),
        ("C", "21"),
        ("D", "12"),
    ];
    assert_eq!(
        offsets,
        want.iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect::<Vec<_>>()
    );
}

#[test]
fn a_gap_costs_no_loss_but_the_same_taxon_aligned_through_does() {
    let w = scratch("gap");
    let m = example("02_missing_data");
    let mut losses = Vec::new();
    for variant in ["gapped", "full"] {
        let gff = w.join(format!("{variant}.gff3"));
        let mut text = fs::read_to_string(m.join("abc.gff3")).unwrap();
        text.push_str(&fs::read_to_string(m.join(format!("d.{variant}.gff3"))).unwrap());
        fs::write(&gff, text).unwrap();
        let out = w.join(variant);
        let s = ok(&[
            "analyze",
            "--gff",
            gff.to_str().unwrap(),
            "-p",
            m.join(format!("proteins.{variant}.faa")).to_str().unwrap(),
            "--msa",
            m.join(format!("family.{variant}.aln.fasta"))
                .to_str()
                .unwrap(),
            "-t",
            m.join("family.nwk").to_str().unwrap(),
            "-m",
            "dollo",
            "-o",
            out.to_str().unwrap(),
        ]);
        losses.push(s);
    }
    assert!(losses[0].contains("dollo 2g/0l"), "{}", losses[0]);
    assert!(losses[1].contains("dollo 2g/1l"), "{}", losses[1]);
}

#[test]
fn matrix_then_score_equals_analyze() {
    let w = scratch("split");
    let q = example("01_quickstart");
    let (gff, faa, aln, nwk) = (
        q.join("family.gff3"),
        q.join("proteins.faa"),
        q.join("family.aln.fasta"),
        q.join("family.nwk"),
    );
    let m_out = w.join("m");
    ok(&[
        "matrix",
        "--gff",
        gff.to_str().unwrap(),
        "-p",
        faa.to_str().unwrap(),
        "--msa",
        aln.to_str().unwrap(),
        "-o",
        m_out.to_str().unwrap(),
    ]);
    let s_out = w.join("s");
    ok(&[
        "score",
        "--matrix",
        m_out.join("family.aln/matrix.tsv").to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "-m",
        "dollo,mp",
        "-o",
        s_out.to_str().unwrap(),
    ]);
    for f in ["dollo_sites.tsv", "mp_sites.tsv"] {
        let a = fs::read_to_string(example("01_quickstart").join("expected/family.aln").join(f))
            .unwrap();
        let b = fs::read_to_string(s_out.join(f)).unwrap();
        assert_eq!(
            a, b,
            "{f}: score from a matrix file differs from the one-pass analyze"
        );
    }
}

#[test]
fn min_observed_drops_sites_that_are_mostly_missing() {
    let w = scratch("sparse");
    let e = example("02_missing_data").join("expected/gapped/family.gapped.aln/matrix.tsv");
    let nwk = example("02_missing_data").join("family.nwk");
    let out = w.join("o");
    ok(&[
        "score",
        "--matrix",
        e.to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "-m",
        "dollo",
        "--min-observed",
        "4",
        "-o",
        out.to_str().unwrap(),
    ]);
    let summary = fs::read_to_string(out.join("summary_score.tsv")).unwrap();
    assert!(summary.contains("characters_too_sparse\t1"), "{summary}");
    assert!(summary.contains("characters_scored\t1"), "{summary}");
}

#[test]
fn gzip_inputs_are_read_transparently() {
    use flate2::{write::GzEncoder, Compression};
    let w = scratch("gz");
    let q = example("01_quickstart");
    let gz = w.join("family.gff3.gz");
    let mut enc = GzEncoder::new(fs::File::create(&gz).unwrap(), Compression::default());
    enc.write_all(&fs::read(q.join("family.gff3")).unwrap())
        .unwrap();
    enc.finish().unwrap();
    let out = w.join("o");
    let s = ok(&[
        "analyze",
        "--gff",
        gz.to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        q.join("family.aln.fasta").to_str().unwrap(),
        "-t",
        q.join("family.nwk").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(s.contains("dollo 2g/2l"), "{s}");
}

#[test]
fn demo_runs_and_checks_itself() {
    let w = scratch("demo");
    let s = ok(&["demo", "-o", w.join("d").to_str().unwrap()]);
    assert!(s.contains("expect 3 sites, 3 Dollo gains, 2 Dollo losses, 3 mp steps: got 3 sites, 3 gains, 2 losses, 3 steps"), "{s}");
}

#[test]
fn malin_export_requires_phase_and_shared_taxa() {
    let w = scratch("malin_err");
    // A 0/1 table imports without phase, so it cannot be exported again.
    let bin_table = repo().join("tests/data/malin_binary.table");
    ok(&[
        "convert",
        "--from-malin",
        bin_table.to_str().unwrap(),
        "-o",
        w.join("m").to_str().unwrap(),
    ]);
    let m = w.join("m/malin_binary/matrix.tsv");
    assert!(m.exists(), "{:?}", files_under(&w));
    let o = run(&[
        "convert",
        "--to-malin",
        m.to_str().unwrap(),
        "-o",
        w.join("x.table").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("no phase"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

#[test]
fn errors_are_clear_and_use_distinct_exit_codes() {
    // Runtime failure: exit 1, message on stderr, nothing on stdout.
    let o = run(&[
        "analyze",
        "--gff",
        "/no/such.gff3",
        "-p",
        "/no/such.faa",
        "--msa",
        "/no/such.aln",
        "-t",
        "/no/such.nwk",
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).starts_with("error:"));
    // Missing required input.
    let o = run(&["project", "-p", "x.faa"]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("--gff is required"));
    // Usage error: clap's exit code 2.
    assert_eq!(run(&["no-such-command"]).status.code(), Some(2));
    assert_eq!(run(&["score"]).status.code(), Some(2));
}

#[test]
fn a_tree_that_shares_no_taxon_with_the_matrix_fails_loudly() {
    let w = scratch("names");
    let q = example("01_quickstart");
    let nwk = w.join("wrong.nwk");
    fs::write(&nwk, "((X,Y),(Z,W));\n").unwrap();
    let out = w.join("o");
    let o = run(&[
        "analyze",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        q.join("family.aln.fasta").to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("no matrix taxon matches a tree tip"));
}

#[test]
fn two_alignments_with_the_same_name_are_refused_instead_of_overwritten() {
    let w = scratch("stem");
    let q = example("01_quickstart");
    fs::create_dir_all(w.join("a")).unwrap();
    fs::create_dir_all(w.join("b")).unwrap();
    fs::copy(q.join("family.aln.fasta"), w.join("a/fam.fasta")).unwrap();
    fs::copy(q.join("family.aln.fasta"), w.join("b/fam.fasta")).unwrap();
    let o = run(&[
        "matrix",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        w.join("a/fam.fasta").to_str().unwrap(),
        "--msa",
        w.join("b/fam.fasta").to_str().unwrap(),
        "-o",
        w.join("o").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("overwritten"));
}

#[test]
fn quickstart_variations_documented_in_the_example_readme_hold() {
    let w = scratch("variations");
    let q = example("01_quickstart");
    let base = |out: &Path, extra: &[&str]| {
        let mut a = vec![
            "analyze".to_string(),
            "--gff".into(),
            q.join("family.gff3").to_str().unwrap().into(),
            "-p".into(),
            q.join("proteins.faa").to_str().unwrap().into(),
            "--msa".into(),
            q.join("family.aln.fasta").to_str().unwrap().into(),
            "-t".into(),
            q.join("family.nwk").to_str().unwrap().into(),
            "-o".into(),
            out.to_str().unwrap().into(),
        ];
        a.extend(extra.iter().map(|s| s.to_string()));
        let refs: Vec<&str> = a.iter().map(String::as_str).collect();
        ok(&refs)
    };
    // Rooting on D changes Dollo (site 7:0 becomes a clade) but not Fitch.
    let s = base(&w.join("og"), &["--outgroup", "D"]);
    assert!(s.contains("dollo 2g/1l, mp 2 steps"), "{s}");
    let summary = fs::read_to_string(w.join("og/family.aln/summary_score.tsv")).unwrap();
    assert!(summary.contains("rooting\toutgroup:D"), "{summary}");
    // --min-taxa 4 leaves nothing to score.
    base(&w.join("mt"), &["--min-taxa", "4"]);
    let summary = fs::read_to_string(w.join("mt/family.aln/summary_score.tsv")).unwrap();
    assert!(
        summary.contains("characters_scored\t0") && summary.contains("characters_skipped\t2"),
        "{summary}"
    );
}

#[test]
fn scoring_with_dollo_but_no_tree_is_a_clear_error() {
    let w = scratch("notree");
    let q = example("01_quickstart");
    let o = run(&[
        "analyze",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        q.join("family.aln.fasta").to_str().unwrap(),
        "-o",
        w.join("o").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("--tree is required"));
    // ...but the tree-free tables can still be produced.
    ok(&[
        "analyze",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        q.join("family.aln.fasta").to_str().unwrap(),
        "-m",
        "counts",
        "-o",
        w.join("c").to_str().unwrap(),
    ]);
    assert!(w.join("c/family.aln/counts_sites.tsv").exists());
}

/// Rewrites example 1 so every species' protein has a different, family-specific id (as in a real
/// one-alignment-per-gene-family run), keeping the tree on species names.
fn renamed_quickstart(w: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    let q = example("01_quickstart");
    let rename = |s: String| {
        ["A", "B", "C", "D"].iter().fold(s, |acc, sp| {
            acc.replace(
                &format!("protein_id={sp}\n"),
                &format!("protein_id=fam1_{sp}\n"),
            )
        })
    };
    let gff = w.join("renamed.gff3");
    fs::write(
        &gff,
        rename(fs::read_to_string(q.join("family.gff3")).unwrap()),
    )
    .unwrap();
    let relabel = |text: String| {
        ["A", "B", "C", "D"].iter().fold(text, |acc, sp| {
            acc.replace(&format!(">{sp}\n"), &format!(">fam1_{sp}\n"))
        })
    };
    let faa = w.join("renamed.faa");
    fs::write(
        &faa,
        relabel(fs::read_to_string(q.join("proteins.faa")).unwrap()),
    )
    .unwrap();
    let aln = w.join("renamed.aln.fasta");
    fs::write(
        &aln,
        relabel(fs::read_to_string(q.join("family.aln.fasta")).unwrap()),
    )
    .unwrap();
    let map = w.join("taxa.tsv");
    fs::write(&map, "fam1_A\tA\nfam1_B\tB\nfam1_C\tC\nfam1_D\tD\n").unwrap();
    (gff, faa, aln, map)
}

#[test]
fn taxon_map_lets_family_specific_protein_ids_use_a_species_tree() {
    let w = scratch("taxonmap");
    let (gff, faa, aln, map) = renamed_quickstart(&w);
    let nwk = example("01_quickstart").join("family.nwk");
    // Without the map the tree tips (A..D) match nothing: a clear error that names the fix.
    let o = run(&[
        "analyze",
        "--gff",
        gff.to_str().unwrap(),
        "-p",
        faa.to_str().unwrap(),
        "--msa",
        aln.to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "-o",
        w.join("no").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("--taxon-map"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    // With it, the result equals the original example.
    let s = ok(&[
        "analyze",
        "--gff",
        gff.to_str().unwrap(),
        "-p",
        faa.to_str().unwrap(),
        "--msa",
        aln.to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "--taxon-map",
        map.to_str().unwrap(),
        "-o",
        w.join("yes").to_str().unwrap(),
    ]);
    assert!(s.contains("dollo 2g/2l, mp 2 steps"), "{s}");
    let m = fs::read_to_string(w.join("yes/renamed.aln/matrix.tsv")).unwrap();
    assert!(
        m.lines()
            .any(|l| l.starts_with("site\t") && l.ends_with("\tA\tB\tC\tD")),
        "{m}"
    );
}

#[test]
fn two_rows_mapped_to_one_taxon_are_refused() {
    let w = scratch("taxondup");
    let (gff, faa, aln, _) = renamed_quickstart(&w);
    let map = w.join("dup.tsv");
    fs::write(&map, "fam1_A\tX\nfam1_B\tX\n").unwrap();
    let o = run(&[
        "matrix",
        "--gff",
        gff.to_str().unwrap(),
        "-p",
        faa.to_str().unwrap(),
        "--msa",
        aln.to_str().unwrap(),
        "--taxon-map",
        map.to_str().unwrap(),
        "-o",
        w.join("o").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("both map to taxon"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

#[test]
fn malin_import_never_writes_outside_the_output_folder_or_overwrites() {
    let w = scratch("malin_stems");
    // `..` must not name the parent folder, and `a/b` / `A_B` (after sanitising, case-insensitively)
    // must not share one folder.
    let ok_table = w.join("dots.table");
    fs::write(&ok_table, "A\t31\nB\t00\n#MAP\t0\t..\t3\n#MAP\t1\tx\t7\n").unwrap();
    let out = w.join("o1");
    ok(&[
        "convert",
        "--from-malin",
        ok_table.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        out.join("alignment/matrix.tsv").exists(),
        "{:?}",
        files_under(&w)
    );
    assert!(!w.join("matrix.tsv").exists());

    let clash = w.join("clash.table");
    fs::write(&clash, "A\t31\nB\t00\n#MAP\t0\ta/b\t3\n#MAP\t1\tA_B\t7\n").unwrap();
    let o = run(&[
        "convert",
        "--from-malin",
        clash.to_str().unwrap(),
        "-o",
        w.join("o2").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("overwritten"));
}

#[test]
fn output_is_identical_whatever_the_thread_count() {
    // Rayon splits work differently at 1 and 8 threads; every output file must still match byte for byte.
    let q = example("02_missing_data");
    let w = scratch("threads");
    let mut text = fs::read_to_string(q.join("abc.gff3")).unwrap();
    text.push_str(&fs::read_to_string(q.join("d.full.gff3")).unwrap());
    let gff = w.join("f.gff3");
    fs::write(&gff, text).unwrap();
    let run_with = |threads: &str, out: &Path| {
        let o = Command::new(BIN)
            .env("RAYON_NUM_THREADS", threads)
            .args([
                "analyze",
                "--gff",
                gff.to_str().unwrap(),
                "-p",
                q.join("proteins.full.faa").to_str().unwrap(),
                "--msa",
                q.join("family.full.aln.fasta").to_str().unwrap(),
                "-t",
                q.join("family.nwk").to_str().unwrap(),
                "--reconstructions",
                "-o",
                out.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    };
    let (one, eight) = (w.join("t1"), w.join("t8"));
    run_with("1", &one);
    run_with("8", &eight);
    let (a, b) = (files_under(&one), files_under(&eight));
    assert_eq!(a, b);
    for f in a {
        let (x, y) = (
            fs::read(one.join(&f)).unwrap(),
            fs::read(eight.join(&f)).unwrap(),
        );
        // summary_score.tsv embeds the (different) output path; everything else must be identical.
        if f.file_name().and_then(|n| n.to_str()) == Some("summary_score.tsv") {
            continue;
        }
        assert_eq!(x, y, "{} differs between 1 and 8 threads", f.display());
    }
}

#[test]
fn verbose_prints_machine_readable_stage_timings() {
    let w = scratch("timers");
    let q = example("01_quickstart");
    let o = run(&[
        "-v",
        "analyze",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "--msa",
        q.join("family.aln.fasta").to_str().unwrap(),
        "-t",
        q.join("family.nwk").to_str().unwrap(),
        "-o",
        w.join("o").to_str().unwrap(),
    ]);
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    for stage in [
        "proteins_read",
        "gff_parse_and_join",
        "join_finish",
        "alignments_read",
        "alignments_project",
        "matrix_build",
    ] {
        let line = err
            .lines()
            .find(|l| l.starts_with(&format!("[time] {stage}\t")));
        let secs: f64 = line
            .unwrap_or_else(|| panic!("no [time] line for {stage}: {err}"))
            .split('\t')
            .nth(1)
            .unwrap()
            .parse()
            .unwrap();
        assert!(secs >= 0.0);
    }
    assert!(err
        .lines()
        .any(|l| l.starts_with("[size] gff_bytes_on_disk\t")));
}

#[cfg(target_os = "linux")]
#[test]
fn verbose_reports_resident_memory_per_stage_on_linux() {
    let w = scratch("memprobe");
    let q = example("01_quickstart");
    let o = run(&[
        "-v",
        "project",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "-o",
        w.join("o").to_str().unwrap(),
    ]);
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    let line = err
        .lines()
        .find(|l| l.starts_with("[mem] join_annotations\t"))
        .unwrap_or_else(|| panic!("{err}"));
    let f: Vec<u64> = line
        .split('\t')
        .skip(1)
        .map(|v| v.parse().unwrap())
        .collect();
    assert!(f.len() == 2 && f[0] > 0 && f[1] >= f[0], "{line}");
}

#[test]
fn annotation_diagnostics_reach_diagnostics_tsv() {
    // These were raised by the annotation reader but never reported (the loader dropped them).
    let w = scratch("anndiag");
    let gff = w.join("a.gff3");
    fs::write(
        &gff,
        "c\tt\tCDS\t100\t114\t.\t+\t0\tID=a1;Parent=tA;protein_id=PA\n\
c\tt\tCDS\t200\t214\t.\t+\t2\tID=a2;Parent=tA;protein_id=PA\n\
c\tt\texon\t300\t320\t.\t+\t.\tID=n1;Parent=tNoCds;protein_id=PN\n\
c\tt\texon\t400\t420\t.\t+\t.\tID=n2;Parent=tNoCds;protein_id=PN\n",
    )
    .unwrap();
    let faa = w.join("p.faa");
    fs::write(&faa, ">PA\nMKVLGAWEHTRSMKVLG\n>PN\nMKVLGA\n").unwrap();
    let out = w.join("o");
    ok(&[
        "project",
        "--gff",
        gff.to_str().unwrap(),
        "-p",
        faa.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let tally = fs::read_to_string(out.join("diagnostics.tsv")).unwrap();
    for kind in ["frame_mismatch", "no_cds"] {
        assert!(
            tally.lines().any(|l| l.starts_with(&format!("{kind}\t"))),
            "{kind} missing from:\n{tally}"
        );
    }
}

#[test]
fn a_clean_run_writes_no_diagnostics_file() {
    let w = scratch("cleandiag");
    let q = example("01_quickstart");
    let out = w.join("o");
    ok(&[
        "project",
        "--gff",
        q.join("family.gff3").to_str().unwrap(),
        "-p",
        q.join("proteins.faa").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        !out.join("diagnostics.tsv").exists(),
        "info-only runs write no report"
    );
}

/// The demo writes a GFF, proteins, an alignment and a tree; lay them out so every input form can be exercised.
fn demo_inputs(name: &str) -> PathBuf {
    let w = scratch(name);
    ok(&["demo", "-o", w.join("demo").to_str().unwrap()]);
    let d = w.join("in");
    fs::create_dir_all(d.join("aln")).unwrap();
    fs::copy(w.join("demo/demo.gff3"), d.join("a.gff3")).unwrap();
    fs::copy(w.join("demo/demo.faa"), d.join("a.faa")).unwrap();
    fs::copy(w.join("demo/demo.aln.fasta"), d.join("aln/one.aln.fasta")).unwrap();
    fs::copy(w.join("demo/demo.aln.fasta"), d.join("aln/two.aln.fasta")).unwrap();
    fs::write(d.join("aln/notes.txt"), "not an alignment").unwrap();
    fs::copy(w.join("demo/demo.nwk"), d.join("t.nwk")).unwrap();
    w
}

#[test]
fn msa_dir_and_msa_list_give_the_same_result_as_repeated_msa() {
    let w = demo_inputs("msa_forms");
    let d = w.join("in");
    let (gff, faa, nwk) = (
        d.join("a.gff3").to_str().unwrap().to_string(),
        d.join("a.faa").to_str().unwrap().to_string(),
        d.join("t.nwk").to_str().unwrap().to_string(),
    );
    let base = ["analyze", "--gff", &gff, "-p", &faa, "-t", &nwk];
    let go = |extra: &[&str], out: &str| {
        let mut a: Vec<&str> = base.to_vec();
        a.extend_from_slice(extra);
        let o = w.join(out);
        a.extend_from_slice(&["-o", o.to_str().unwrap()]);
        ok(&a);
        o
    };
    let (one, two) = (d.join("aln/one.aln.fasta"), d.join("aln/two.aln.fasta"));
    let explicit = go(
        &[
            "--msa",
            one.to_str().unwrap(),
            "--msa",
            two.to_str().unwrap(),
        ],
        "explicit",
    );
    let by_dir = go(&["--msa-dir", d.join("aln").to_str().unwrap()], "bydir");
    fs::write(
        d.join("fam.list"),
        "# two families\naln/one.aln.fasta\naln/two.aln.fasta\n",
    )
    .unwrap();
    let by_list = go(
        &["--msa-list", d.join("fam.list").to_str().unwrap()],
        "bylist",
    );
    for out in [&by_dir, &by_list] {
        for f in [
            "one.aln/matrix.tsv",
            "two.aln/matrix.tsv",
            "one.aln/dollo_branches.tsv",
        ] {
            assert_eq!(
                fs::read(explicit.join(f)).unwrap(),
                fs::read(out.join(f)).unwrap(),
                "{f} differs for {out:?}"
            );
        }
    }
    // an extension that matches nothing is an error naming what was tried
    let o = run(&[
        "matrix",
        "--gff",
        d.join("a.gff3").to_str().unwrap(),
        "-p",
        d.join("a.faa").to_str().unwrap(),
        "--msa-dir",
        d.join("aln").to_str().unwrap(),
        "--msa-ext",
        "stk",
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("no alignment files"));
}

#[test]
fn gff_and_protein_folders_and_lists_work() {
    let w = demo_inputs("ann_forms");
    let d = w.join("in");
    let o = w.join("out");
    ok(&[
        "project",
        "--gff-dir",
        d.to_str().unwrap(),
        "--proteins-dir",
        d.to_str().unwrap(),
        "-o",
        o.to_str().unwrap(),
    ]);
    assert!(o.join("introns.tsv").exists());
    fs::write(d.join("g.list"), "a.gff3\n").unwrap();
    fs::write(d.join("p.list"), "a.faa\n").unwrap();
    let o2 = w.join("out2");
    ok(&[
        "project",
        "--gff-list",
        d.join("g.list").to_str().unwrap(),
        "--proteins-list",
        d.join("p.list").to_str().unwrap(),
        "-o",
        o2.to_str().unwrap(),
    ]);
    assert_eq!(
        fs::read(o.join("introns.tsv")).unwrap(),
        fs::read(o2.join("introns.tsv")).unwrap()
    );
}

#[test]
fn score_matrix_dir_scores_every_family_into_its_own_folder() {
    let w = demo_inputs("score_dir");
    let d = w.join("in");
    let fam = w.join("fam");
    ok(&[
        "matrix",
        "--gff",
        d.join("a.gff3").to_str().unwrap(),
        "-p",
        d.join("a.faa").to_str().unwrap(),
        "--msa-dir",
        d.join("aln").to_str().unwrap(),
        "-o",
        fam.to_str().unwrap(),
    ]);
    let scored = w.join("scored");
    let out = ok(&[
        "score",
        "--matrix-dir",
        fam.to_str().unwrap(),
        "-t",
        d.join("t.nwk").to_str().unwrap(),
        "-o",
        scored.to_str().unwrap(),
    ]);
    assert!(out.contains("scored 2 matrices"), "{out}");
    for f in ["one.aln", "two.aln"] {
        assert!(
            scored.join(f).join("dollo_branches.tsv").exists(),
            "{:?}",
            files_under(&scored)
        );
    }
    // a single --matrix keeps the old layout: results directly in -o
    let single = w.join("single");
    ok(&[
        "score",
        "--matrix",
        fam.join("one.aln/matrix.tsv").to_str().unwrap(),
        "-t",
        d.join("t.nwk").to_str().unwrap(),
        "-o",
        single.to_str().unwrap(),
    ]);
    assert!(single.join("dollo_branches.tsv").exists());
    assert_eq!(
        fs::read(single.join("dollo_branches.tsv")).unwrap(),
        fs::read(scored.join("one.aln/dollo_branches.tsv")).unwrap()
    );
}

fn branch_totals(path: &Path) -> (u64, u64) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let (g, l) = (
        head.iter().position(|h| *h == "gains").unwrap(),
        head.iter().position(|h| *h == "losses").unwrap(),
    );
    lines.fold((0, 0), |(ag, al), line| {
        let c: Vec<&str> = line.split('\t').collect();
        (
            ag + c[g].parse::<u64>().unwrap(),
            al + c[l].parse::<u64>().unwrap(),
        )
    })
}

#[test]
fn combine_stacks_every_family_and_scores_the_total() {
    let w = demo_inputs("combine");
    let d = w.join("in");
    let (gff, faa, nwk, aln) = (
        d.join("a.gff3"),
        d.join("a.faa"),
        d.join("t.nwk"),
        d.join("aln"),
    );
    let out = w.join("out");
    let stdout = ok(&[
        "analyze",
        "--gff",
        gff.to_str().unwrap(),
        "-p",
        faa.to_str().unwrap(),
        "--msa-dir",
        aln.to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "--combine",
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        stdout.contains("combined: ") && stdout.contains("2 alignments"),
        "{stdout}"
    );
    let one = branch_totals(&out.join("one.aln/dollo_branches.tsv"));
    let two = branch_totals(&out.join("two.aln/dollo_branches.tsv"));
    let both = branch_totals(&out.join("combined/dollo_branches.tsv"));
    assert_eq!(
        both,
        (one.0 + two.0, one.1 + two.1),
        "combined totals are the sums"
    );
    assert!(both.0 > 0);
    let sites = fs::read_to_string(out.join("combined/dollo_sites.tsv")).unwrap();
    assert!(
        sites.contains("one.aln:") && sites.contains("two.aln:"),
        "{sites}"
    );
    let m = fs::read_to_string(out.join("combined/matrix.tsv")).unwrap();
    assert!(m
        .lines()
        .nth(1)
        .unwrap()
        .starts_with("site\tcolumn\tphase\tfamily\tfamily_column\t"));
    let n = |p: &str| {
        fs::read_to_string(out.join(p))
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .count()
            - 1
    };
    assert_eq!(
        n("combined/matrix.tsv"),
        n("one.aln/matrix.tsv") + n("two.aln/matrix.tsv")
    );
}

#[test]
fn score_combine_of_families_matches_analyze_combine() {
    let w = demo_inputs("combine_score");
    let d = w.join("in");
    let fam = w.join("fam");
    ok(&[
        "matrix",
        "--gff",
        d.join("a.gff3").to_str().unwrap(),
        "-p",
        d.join("a.faa").to_str().unwrap(),
        "--msa-dir",
        d.join("aln").to_str().unwrap(),
        "-o",
        fam.to_str().unwrap(),
    ]);
    let list = w.join("fam.list");
    fs::write(&list, "fam/one.aln/matrix.tsv\nfam/two.aln/matrix.tsv\n").unwrap();
    let nwk = d.join("t.nwk");
    let a = w.join("a");
    ok(&[
        "score",
        "--matrix-list",
        list.to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "--combine",
        "-o",
        a.to_str().unwrap(),
    ]);
    let b = w.join("b");
    ok(&[
        "analyze",
        "--gff",
        d.join("a.gff3").to_str().unwrap(),
        "-p",
        d.join("a.faa").to_str().unwrap(),
        "--msa-dir",
        d.join("aln").to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "--combine",
        "-o",
        b.to_str().unwrap(),
    ]);
    for f in ["dollo_branches.tsv", "dollo_sites.tsv", "matrix.tsv"] {
        let norm = |p: PathBuf| {
            fs::read_to_string(p)
                .unwrap()
                .lines()
                .filter(|l| !l.starts_with('#'))
                .collect::<Vec<_>>()
                .join("\n")
        };
        assert_eq!(norm(a.join(f)), norm(b.join("combined").join(f)), "{f}");
    }
    // a tolerance that could reach across alignments is refused
    let o = run(&[
        "score",
        "--matrix-list",
        list.to_str().unwrap(),
        "-t",
        nwk.to_str().unwrap(),
        "--combine",
        "--tolerance",
        "10000",
        "-o",
        w.join("c").to_str().unwrap(),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("not below"));
}

#[test]
fn combined_only_writes_just_the_combined_folder_and_the_same_totals() {
    let w = demo_inputs("combined_only");
    let d = w.join("in");
    let args = |out: &Path, extra: &[&str]| {
        let mut a: Vec<String> = [
            "analyze",
            "--gff",
            d.join("a.gff3").to_str().unwrap(),
            "-p",
            d.join("a.faa").to_str().unwrap(),
            "--msa-dir",
            d.join("aln").to_str().unwrap(),
            "-t",
            d.join("t.nwk").to_str().unwrap(),
            "--combine",
            "-o",
            out.to_str().unwrap(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        a.extend(extra.iter().map(|s| s.to_string()));
        a
    };
    let (full, lean) = (w.join("full"), w.join("lean"));
    for (out, extra) in [(&full, vec![]), (&lean, vec!["--combined-only"])] {
        let a = args(out, &extra);
        ok(&a.iter().map(|s| s.as_str()).collect::<Vec<_>>());
    }
    assert!(full.join("one.aln/matrix.tsv").exists());
    assert!(!lean.join("one.aln").exists() && !lean.join("two.aln").exists());
    assert!(lean.join("introns.tsv").exists() && lean.join("unmatched.tsv").exists());
    for f in [
        "matrix.tsv",
        "dollo_branches.tsv",
        "dollo_sites.tsv",
        "pairwise.tsv",
    ] {
        assert_eq!(
            fs::read(full.join("combined").join(f)).unwrap(),
            fs::read(lean.join("combined").join(f)).unwrap(),
            "{f}"
        );
    }
    // --combined-only means nothing without --combine: a usage error
    let o = run(&[
        "matrix",
        "--gff",
        "a.gff3",
        "-p",
        "a.faa",
        "--msa",
        "x.aln",
        "--combined-only",
    ]);
    assert_eq!(o.status.code(), Some(2));
}
