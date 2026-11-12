//! End-to-end checks that alignment gaps are scored as missing data, not as intron loss.
//!
//! These drive the public pipeline (join -> project -> matrix -> file -> score) on a tiny hand-built
//! family so the expected numbers can be verified by reading the fixture.

use std::fs;
use std::io::Cursor;
use std::path::Path;

use intron_core::matrix::{Cell, SiteMatrix};
use intron_core::model::{AnnotatedIntron, Strand, TranscriptAnnotation};
use intron_core::msa::Alignment;
use intron_core::project::{join_annotations, project_into_alignment, JoinOptions, ProjectOptions};
use intron_core::score::{prepare, score_dollo, score_mp};
use intron_core::seqio::ProteinSet;
use intron_core::tables::matrix_table;
use intron_core::tree::Tree;

fn tx(protein: &str, offsets: &[u32]) -> TranscriptAnnotation {
    TranscriptAnnotation {
        transcript_id: format!("{protein}.t"),
        gene_id: None,
        protein_id: Some(protein.into()),
        seqid: "chr".into(),
        strand: Strand::Plus,
        cds_len: 3 * 12,
        cds_segments: offsets.len() as u32 + 1,
        introns: offsets
            .iter()
            .enumerate()
            .map(|(i, &o)| AnnotatedIntron {
                start: 100 + 50 * i as u64,
                end: 120 + 50 * i as u64,
                rank: i as u32,
                cds_offset: Some(o),
            })
            .collect(),
        frame_mismatches: 0,
        source: "test".into(),
    }
}

/// Run the whole pipeline for four taxa. `d_row`/`d_protein` let a test choose whether D is gapped
/// where the other three share an intron.
fn matrix(d_protein: &str, d_row: &str, d_offsets: &[u32], flank: u32) -> SiteMatrix {
    let fasta = format!(">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLGAWEHTRS\n>D\n{d_protein}\n");
    let mut ps = ProteinSet::new();
    ps.read_one(Cursor::new(fasta.into_bytes()), Path::new("p"))
        .unwrap();
    let aln = Alignment::from_fasta_str(
        &format!(">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLGAWEHTRS\n>D\n{d_row}\n"),
        "fam",
    )
    .unwrap();
    // A, B and C share a phase-1 intron inside residue index 5 (offset 16).
    let txs = vec![
        tx("A", &[16]),
        tx("B", &[16]),
        tx("C", &[16]),
        tx("D", d_offsets),
    ];
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&txs, &ps, &join).unwrap();
    let opts = ProjectOptions {
        flank,
        ..Default::default()
    };
    let ap = project_into_alignment(&aln, &proj, &ps, &opts);
    SiteMatrix::build(&ap)
}

fn site_index(m: &SiteMatrix, col: u32, phase: u8) -> usize {
    m.sites
        .iter()
        .position(|s| s.col == col && s.phase == phase)
        .unwrap_or_else(|| panic!("no site {col}:{phase}"))
}

fn dollo_losses(m: &SiteMatrix) -> u64 {
    let t = Tree::from_newick("((A,B),(C,D));").unwrap();
    let (data, _) = prepare(m, &t, &[], 1).unwrap();
    score_dollo(&data, false).total_losses
}

const FULL: &str = "MKVLGAWEHTRS";
/// D lacks the residues that A/B/C have in columns 4-5 (G, A), so the shared intron's residue
/// (index 5, column 5) is a gap in D.
const GAPPED_ROW: &str = "MKVL--WEHTRS";
const GAPPED_PROTEIN: &str = "MKVLWEHTRS";

#[test]
fn gap_at_the_shared_intron_is_missing_and_costs_no_loss() {
    // D has its own intron at residue 7 of its protein (col 9, phase 0) so it is annotated.
    let m = matrix(GAPPED_PROTEIN, GAPPED_ROW, &[21], 1);
    let s = site_index(&m, 5, 1);
    assert_eq!(m.sites[s].cell(0), Cell::Present);
    assert_eq!(
        m.sites[s].cell(3),
        Cell::Missing,
        "D has no residue at column 5"
    );
    assert_eq!(
        dollo_losses(&m),
        0,
        "a gapped taxon must not be scored as a loss"
    );
}

#[test]
fn the_same_taxon_aligned_through_the_site_is_a_real_loss() {
    // Identical, except D's row now has a residue at column 5, so its absence is an observation.
    let m = matrix(FULL, FULL, &[27], 1);
    let s = site_index(&m, 5, 1);
    assert_eq!(m.sites[s].cell(3), Cell::Absent);
    assert_eq!(
        dollo_losses(&m),
        1,
        "observed absence on D's branch is one loss"
    );
}

#[test]
fn missing_cells_survive_the_matrix_file_and_score_identically() {
    let m = matrix(GAPPED_PROTEIN, GAPPED_ROW, &[21], 1);
    let dir = std::env::temp_dir().join(format!("ip_percell_{}", std::process::id()));
    let path = dir.join("m.tsv");
    matrix_table(&m).write(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("format=2"), "{text}");
    assert!(
        text.lines().any(|l| l.contains("\t111?")),
        "D must be written as ?: {text}"
    );
    let back = SiteMatrix::read_tsv(&path).unwrap();
    let s = site_index(&back, 5, 1);
    assert_eq!(back.sites[s].cell(3), Cell::Missing);
    assert_eq!(back.pattern(s), "111?");
    assert_eq!(dollo_losses(&back), dollo_losses(&m));
}

#[test]
fn flank_zero_still_requires_the_anchoring_residue_only() {
    // D's gap sits on column 4 only (G): the intron's own residue (column 5) is present, so with
    // flank 0 D is observable and absent (a loss); flank 1 needs column 4 too, so D becomes missing.
    let row = "MKVL-AWEHTRS";
    let protein = "MKVLAWEHTRS";
    let strict = matrix(protein, row, &[24], 1);
    let loose = matrix(protein, row, &[24], 0);
    let s = site_index(&strict, 5, 1);
    assert_eq!(strict.sites[s].cell(3), Cell::Missing);
    let l = site_index(&loose, 5, 1);
    assert_eq!(loose.sites[l].cell(3), Cell::Absent);
    assert_eq!(dollo_losses(&strict), 0);
    assert_eq!(dollo_losses(&loose), 1);
}

#[test]
fn phase_zero_beside_an_indel_is_missing_not_split_across_columns() {
    // A, B, C share a boundary intron after residue 5 (offset 18 -> column 6, phase 0). D lacks
    // residue 5 (column 5), so its residues around the boundary are in columns 4 and 6: not adjacent.
    // Anchoring "to the following residue" would still put D's intron at column 6 if it had one, but
    // its *absence* at the boundary cannot be observed.
    let row = "MKVLG-WEHTRS";
    let protein = "MKVLGWEHTRS";
    let fasta = format!(">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLGAWEHTRS\n>D\n{protein}\n");
    let mut ps = ProteinSet::new();
    ps.read_one(Cursor::new(fasta.into_bytes()), Path::new("p"))
        .unwrap();
    let aln = Alignment::from_fasta_str(
        &format!(">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLGAWEHTRS\n>D\n{row}\n"),
        "fam",
    )
    .unwrap();
    let txs = vec![
        tx("A", &[18]),
        tx("B", &[18]),
        tx("C", &[18]),
        tx("D", &[30]),
    ];
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&txs, &ps, &join).unwrap();
    let ap = project_into_alignment(&aln, &proj, &ps, &ProjectOptions::default());
    let m = SiteMatrix::build(&ap);
    let s = site_index(&m, 6, 0);
    assert_eq!(m.sites[s].cell(0), Cell::Present);
    assert_eq!(m.sites[s].cell(3), Cell::Missing);
}

#[test]
fn a_carrier_not_aligned_around_its_site_is_masked_and_counted() {
    // C carries the intron but its row is gapped next to the anchoring residue. With flank 1 that
    // present call is not an observation: it becomes `?` and is reported, not kept as presence.
    let fasta = ">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLAWEHTRS\n>D\nMKVLGAWEHTRS\n";
    let mut ps = ProteinSet::new();
    ps.read_one(Cursor::new(fasta.as_bytes().to_vec()), Path::new("p"))
        .unwrap();
    let aln = Alignment::from_fasta_str(
        ">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVL-AWEHTRS\n>D\nMKVLGAWEHTRS\n",
        "fam",
    )
    .unwrap();
    let txs = vec![
        tx("A", &[16]),
        tx("B", &[16]),
        tx("C", &[13]),
        tx("D", &[16]),
    ];
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&txs, &ps, &join).unwrap();
    let ap = project_into_alignment(&aln, &proj, &ps, &ProjectOptions::default());
    let m = SiteMatrix::build(&ap);
    assert_eq!(m.masked_present, 1);
    let s = site_index(&m, 5, 1);
    assert_eq!(m.sites[s].cell(2), Cell::Missing);
}

#[test]
fn an_all_missing_taxon_changes_no_score() {
    let base = matrix(FULL, FULL, &[27], 1);
    let t3 = Tree::from_newick("((A,B),C);").unwrap();
    let t4 = Tree::from_newick("((A,B),(C,D));").unwrap();
    // D with no annotation at all is all-missing; compare to a tree that omits D.
    let m = {
        let fasta = ">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLGAWEHTRS\n";
        let mut ps = ProteinSet::new();
        ps.read_one(Cursor::new(fasta.as_bytes().to_vec()), Path::new("p"))
            .unwrap();
        let aln = Alignment::from_fasta_str(
            ">A\nMKVLGAWEHTRS\n>B\nMKVLGAWEHTRS\n>C\nMKVLGAWEHTRS\n>D\nMKVLGAWEHTRS\n",
            "fam",
        )
        .unwrap();
        let txs = vec![tx("A", &[16]), tx("B", &[16]), tx("C", &[16])];
        let join = JoinOptions {
            length_slack: u64::MAX,
            ..Default::default()
        };
        let proj = join_annotations(&txs, &ps, &join).unwrap();
        SiteMatrix::build(&project_into_alignment(
            &aln,
            &proj,
            &ps,
            &ProjectOptions::default(),
        ))
    };
    let (d3, _) = prepare(&m, &t3, &[], 1).unwrap();
    let (d4, _) = prepare(&m, &t4, &[], 1).unwrap();
    assert_eq!(
        score_dollo(&d3, false).total_losses,
        score_dollo(&d4, false).total_losses
    );
    assert_eq!(
        score_mp(&d3, false).total_steps,
        score_mp(&d4, false).total_steps
    );
    let _ = base;
}

#[test]
fn legacy_matrix_files_without_question_marks_still_read() {
    let dir = std::env::temp_dir().join(format!("ip_legacy_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("old.tsv");
    std::fs::write(
        &p,
        "# alignment=x alen=10 phase_mode=strict taxa=3 sites=2 unannotated=C residues=-\n\
site\tcolumn\tphase\tn_taxa\tn_annotated\tpattern\tA\tB\tC\n\
3:0\t3\t0\t1\t2\t100\t1\t0\t0\n\
9:1\t9\t1\t2\t2\t110\t1\t1\t0\n",
    )
    .unwrap();
    let m = SiteMatrix::read_tsv(&p).unwrap();
    assert_eq!(
        m.sites[0].cell(2),
        Cell::Missing,
        "unannotated C is missing, not absent"
    );
    assert_eq!(m.sites[0].cell(1), Cell::Absent);
}

#[test]
fn presence_in_an_unannotated_taxon_is_rejected_and_unsorted_rows_are_fixed() {
    let dir = std::env::temp_dir().join(format!("ip_bad_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bad = dir.join("bad.tsv");
    std::fs::write(
        &bad,
        "# unannotated=B\nsite\tcolumn\tphase\tn_taxa\tn_observable\tpattern\tA\tB\n1:0\t1\t0\t2\t1\t11\t1\t1\n",
    )
    .unwrap();
    assert!(SiteMatrix::read_tsv(&bad).is_err());

    let unsorted = dir.join("unsorted.tsv");
    std::fs::write(
        &unsorted,
        "site\tcolumn\tphase\tn_taxa\tn_observable\tpattern\tA\tB\n\
20:0\t20\t0\t2\t2\t11\t1\t1\n3:0\t3\t0\t2\t2\t11\t1\t1\n5:0\t5\t0\t2\t2\t11\t1\t1\n",
    )
    .unwrap();
    let m = SiteMatrix::read_tsv(&unsorted).unwrap();
    assert_eq!(
        m.sites.iter().map(|s| s.col).collect::<Vec<_>>(),
        vec![3, 5, 20]
    );
    assert_eq!(m.pair_stat(0, 1, 0, true).shared, 3);
}

#[test]
fn a_row_that_is_not_the_joined_protein_is_missing_data_not_a_shifted_projection() {
    // The alignment row lost the residue "K" (e.g. column-trimmed), so protein offsets no longer
    // line up with its columns. Projecting through it would put the site on the wrong column.
    let fasta = ">A\nMKVLGAWE\n>B\nMKVLGAWE\n";
    let mut ps = ProteinSet::new();
    ps.read_one(Cursor::new(fasta.as_bytes().to_vec()), Path::new("p"))
        .unwrap();
    let aln = Alignment::from_fasta_str(">A\nMKVLGAWE\n>B\nM-VLGAWE\n", "fam").unwrap();
    let txs = vec![tx("A", &[15]), tx("B", &[15])];
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&txs, &ps, &join).unwrap();
    let ap = project_into_alignment(&aln, &proj, &ps, &ProjectOptions::default());
    assert_eq!(ap.coverage.row_mismatch, vec!["B".to_string()]);
    assert!(!ap.taxa[1].has_annotation);
    assert_eq!(ap.taxa[1].unannotated_reason, Some("row_sequence_mismatch"));
    assert!(ap.taxa[1].sites.is_empty());
}

#[test]
fn an_annotated_intronless_protein_is_not_reported_as_unannotated() {
    let mut ps = ProteinSet::new();
    ps.read_one(
        Cursor::new(b">A\nMKVLGAWE\n>B\nMKVLGAWE\n".to_vec()),
        Path::new("p"),
    )
    .unwrap();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&[tx("A", &[15]), tx("B", &[])], &ps, &join).unwrap();
    assert!(
        proj.proteins_without_annotation.is_empty(),
        "{:?}",
        proj.proteins_without_annotation
    );
    assert_eq!(proj.joined_proteins, 2);
}

#[test]
fn clamp_does_not_pile_non_tail_sites_onto_the_last_column() {
    use intron_core::project::TailPolicy;
    // The row covers only residues 4.. (Stockholm-style window); a site before the window is not a
    // tail overhang and must be dropped, not snapped to the C-terminal column.
    let mut ps = ProteinSet::new();
    ps.read_one(Cursor::new(b">A\nMKVLGAWE\n".to_vec()), Path::new("p"))
        .unwrap();
    let aln = Alignment::from_stockholm_str(
        "# STOCKHOLM 1.0\n#=GS A SF 5\n#=GS A SE 8\nA GAWE\nB GAWE\n//\n",
        "fam",
    )
    .unwrap();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&[tx("A", &[3])], &ps, &join).unwrap();
    let opts = ProjectOptions {
        tail: TailPolicy::Clamp,
        ..Default::default()
    };
    let ap = project_into_alignment(&aln, &proj, &ps, &opts);
    assert!(
        ap.taxa[0].sites.is_empty(),
        "a site before the aligned window must be dropped under Clamp: {:?}",
        ap.taxa[0].sites
    );
}

#[test]
fn an_internal_node_labelled_root_is_not_mistaken_for_the_real_root() {
    use intron_core::matrix::{Bitset, Site};
    use intron_core::project::PhaseMode;
    use intron_core::tables::tree_with_gain_counts;
    // Rooting `(A,B,(C,D)cd)root;` on A makes the old root an ordinary internal node that is still
    // labelled "root". B and C carry the site, so the gain belongs on that node, not on the new root.
    let taxa: Vec<String> = ["A", "B", "C", "D"].iter().map(|s| s.to_string()).collect();
    let (mut presence, mut observable) = (Bitset::new(4), Bitset::new(4));
    for t in 0..4 {
        observable.set(t);
    }
    presence.set(1);
    presence.set(2);
    let m = SiteMatrix::from_parts(
        "x",
        taxa,
        vec![Site {
            col: 1,
            phase: 0,
            presence,
            observable,
        }],
        PhaseMode::Strict,
    );
    let tree = Tree::from_newick("(A,B,(C,D)cd)root;").unwrap();
    let (data, _) = prepare(&m, &tree, &["A".to_string()], 1).unwrap();
    let d = score_dollo(&data, false);
    let gained = d.sites[0].gain_node.clone().unwrap();
    assert_ne!(
        gained, "root",
        "the gain is on the old root node, which must not be called `root`"
    );
    assert!(gained.starts_with("root#"), "{gained}");
    let real_root = d.branches.iter().find(|b| b.node == 0).unwrap();
    assert_eq!(real_root.name, "root");
    assert_eq!(real_root.gains, 0);
    let nwk = tree_with_gain_counts(&data, &d);
    assert!(
        nwk.matches("[gains=1]").count() == 1 && !nwk.trim_end().ends_with(")[gains=1];"),
        "the gain must not be drawn on the true root: {nwk}"
    );
    let names: Vec<&str> = d.branches.iter().map(|b| b.name.as_str()).collect();
    let mut dedup = names.clone();
    dedup.sort_unstable();
    dedup.dedup();
    assert_eq!(
        dedup.len(),
        names.len(),
        "branch names must be unique: {names:?}"
    );
}

/// Family where A carries a phase-1 intron at column 5 and B is gapped only at column 6.
fn gap_beside_the_site(opts: ProjectOptions) -> SiteMatrix {
    let mut ps = ProteinSet::new();
    ps.read_one(
        Cursor::new(b">A\nMKVLGAWE\n>B\nMKVLGAE\n".to_vec()),
        Path::new("p"),
    )
    .unwrap();
    let aln = Alignment::from_fasta_str(">A\nMKVLGAWE\n>B\nMKVLGA-E\n", "fam").unwrap();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&[tx("A", &[16]), tx("B", &[3])], &ps, &join).unwrap();
    SiteMatrix::build(&project_into_alignment(&aln, &proj, &ps, &opts))
}

#[test]
fn ignore_phase_does_not_turn_a_gap_into_a_loss() {
    use intron_core::project::PhaseMode;
    // Strict: the phase-1 window covers columns 4-6, B is gapped at 6 -> `?`.
    let strict = gap_beside_the_site(ProjectOptions::default());
    let s = site_index(&strict, 5, 1);
    assert_eq!(strict.sites[s].cell(1), Cell::Missing);
    // Ignore: the merged site must not be judged by the narrower phase-0 window, which would
    // make B an observed absence (a loss) that the strict analysis rightly calls missing.
    let ignore = gap_beside_the_site(ProjectOptions {
        phase_mode: PhaseMode::Ignore,
        ..Default::default()
    });
    let s = ignore.sites.iter().position(|x| x.col == 5).unwrap();
    assert_eq!(ignore.sites[s].cell(0), Cell::Present);
    assert_eq!(
        ignore.sites[s].cell(1),
        Cell::Missing,
        "`?` must stay `?` under --ignore-phase"
    );
}

#[test]
fn ignore_phase_keeps_a_carrier_that_strict_keeps() {
    use intron_core::project::PhaseMode;
    // flank 0: a phase-1 carrier only needs its own residue. Under --ignore-phase the carrier must
    // not be dropped because the merged site is (wrongly) treated as a phase-0 boundary needing column c-1.
    let mut ps = ProteinSet::new();
    ps.read_one(
        Cursor::new(b">A\nMKVLAWE\n>B\nMKVLGAWE\n".to_vec()),
        Path::new("p"),
    )
    .unwrap();
    let aln = Alignment::from_fasta_str(">A\nMKVL-AWE\n>B\nMKVLGAWE\n", "fam").unwrap();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&[tx("A", &[13]), tx("B", &[16])], &ps, &join).unwrap();
    for mode in [PhaseMode::Strict, PhaseMode::Ignore] {
        let opts = ProjectOptions {
            phase_mode: mode,
            flank: 0,
            ..Default::default()
        };
        let m = SiteMatrix::build(&project_into_alignment(&aln, &proj, &ps, &opts));
        assert_eq!(
            m.masked_present, 0,
            "{mode:?}: no carrier may be masked at flank 0"
        );
        let s = m.sites.iter().position(|x| x.col == 5).unwrap();
        assert_eq!(m.sites[s].cell(0), Cell::Present, "{mode:?}");
    }
}

#[test]
fn clamped_sites_survive_the_default_flank() {
    use intron_core::project::TailPolicy;
    // A's row ends at column 5 (`MKVLGA--`): an intron anchored to residue 6 lies beyond it and is
    // clamped onto column 5. Its right-hand flank is past the row by construction, so the default
    // flank must not mask it.
    let mut ps = ProteinSet::new();
    ps.read_one(
        Cursor::new(b">A\nMKVLGAWE\n>B\nMKVLGAWE\n".to_vec()),
        Path::new("p"),
    )
    .unwrap();
    let aln = Alignment::from_fasta_str(">A\nMKVLGA--\n>B\nMKVLGAWE\n", "fam").unwrap();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&[tx("A", &[22]), tx("B", &[16])], &ps, &join);
    let proj = proj.unwrap();
    let opts = ProjectOptions {
        tail: TailPolicy::Clamp,
        ..Default::default()
    };
    let ap = project_into_alignment(&aln, &proj, &ps, &opts);
    let m = SiteMatrix::build(&ap);
    assert_eq!(ap.taxa[0].report.clamped, 1);
    assert_eq!(
        m.masked_present, 0,
        "the clamped carrier must not be reported as masked"
    );
    assert!(
        m.sites.iter().any(|s| s.presence.test(0)),
        "A's clamped intron is kept"
    );
}

#[test]
fn isoforms_with_nested_intron_sets_are_flagged() {
    let mut ps = ProteinSet::new();
    ps.read_one(Cursor::new(b">A\nMKVLGAWEHTRS\n".to_vec()), Path::new("p"))
        .unwrap();
    let mut t2 = tx("A", &[3]);
    t2.transcript_id = "A.t2".into();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    // {3,6} and {3}: the smaller set is a subset of the larger, which a "union is bigger" test misses.
    let proj = join_annotations(&[tx("A", &[3, 6]), t2], &ps, &join).unwrap();
    assert_eq!(proj.diagnostics.count_of("isoform_intron_conflict"), 1);
}

#[test]
fn taxon_names_with_commas_or_spaces_survive_the_matrix_file() {
    use intron_core::matrix::{Bitset, Site};
    use intron_core::project::PhaseMode;
    let taxa: Vec<String> = ["s,1", "has space", "ok"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let (mut p, mut o) = (Bitset::new(3), Bitset::new(3));
    p.set(2);
    o.set(2);
    let m = SiteMatrix::from_parts(
        "x",
        taxa,
        vec![Site {
            col: 4,
            phase: 0,
            presence: p,
            observable: o,
        }],
        PhaseMode::Strict,
    );
    assert_eq!(m.annotated, vec![false, false, true]);
    let dir = std::env::temp_dir().join(format!("ip_names_{}", std::process::id()));
    let path = dir.join("m.tsv");
    matrix_table(&m).write(&path).unwrap();
    let back = SiteMatrix::read_tsv(&path).unwrap();
    assert_eq!(
        back.annotated,
        vec![false, false, true],
        "unannotated taxa with odd names must stay unannotated"
    );
    assert_eq!(back.n_annotated, 1);
}

#[test]
fn duplicate_taxon_names_and_huge_columns_are_rejected_not_mangled() {
    use intron_core::malin::{from_malin, to_malin};
    use intron_core::matrix::{Bitset, Site};
    use intron_core::project::PhaseMode;
    let dir = std::env::temp_dir().join(format!("ip_dup_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let p = dir.join("dup.tsv");
    fs::write(
        &p,
        "site\tcolumn\tphase\tn_taxa\tn_observable\tpattern\tA\tA\n1:0\t1\t0\t1\t2\t10\t1\t0\n",
    )
    .unwrap();
    assert!(SiteMatrix::read_tsv(&p).is_err());
    assert!(from_malin("A\t01\nA\t00\n", "x").is_err());
    // A column near u32::MAX cannot be written as a MALIN position (3 * column overflows).
    let (mut pr, mut ob) = (Bitset::new(1), Bitset::new(1));
    pr.set(0);
    ob.set(0);
    let m = SiteMatrix::from_parts(
        "big",
        vec!["A".into()],
        vec![Site {
            col: u32::MAX / 2,
            phase: 0,
            presence: pr,
            observable: ob,
        }],
        PhaseMode::Strict,
    );
    assert!(to_malin(&[&m]).is_err());
    // Two matrices with one name would merge silently when the table is re-imported.
    assert!(to_malin(&[&m, &m]).is_err());
}

#[test]
fn duplicate_rows_in_matrix_tsv_and_malin_merge_without_losing_introns() {
    use intron_core::malin::from_malin;
    let dir = std::env::temp_dir().join(format!("ip_dupmerge_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    // Two rows for site 5:0: A present in the first, B present in the second -> both carry it.
    let p = dir.join("m.tsv");
    fs::write(
        &p,
        "site\tcolumn\tphase\tn_taxa\tn_observable\tpattern\tA\tB\tC\n\
5:0\t5\t0\t1\t2\t10?\t1\t0\t?\n5:0\t5\t0\t1\t2\t?10\t?\t1\t0\n",
    )
    .unwrap();
    let m = SiteMatrix::read_tsv(&p).unwrap();
    assert_eq!(m.n_sites(), 1);
    assert_eq!(
        m.pattern(0),
        "110",
        "presence and observability are both unioned"
    );
    assert_eq!(m.taxon_sites[0], vec![0]);
    assert_eq!(m.taxon_sites[1], vec![0]);
    // The same for two MALIN sites that map to one position.
    let t = "A\t30\nB\t03\nC\t00\n#MAP\t0\tK\t3\n#MAP\t1\tK\t3\n";
    let ms = from_malin(t, "x").unwrap();
    assert_eq!(ms[0].n_sites(), 1);
    assert_eq!(ms[0].pattern(0), "110");
}

#[test]
fn clamp_applies_to_the_first_residue_past_the_end_and_not_before() {
    use intron_core::project::TailPolicy;
    // Row covers residues 0-5. Offset 18 is codon 6 (the first one past the end): clamped to the last
    // covered column. Offset 15 is codon 5 (the last covered residue): placed normally, not clamped.
    let mut ps = ProteinSet::new();
    ps.read_one(
        Cursor::new(b">A\nMKVLGAWE\n>B\nMKVLGAWE\n".to_vec()),
        Path::new("p"),
    )
    .unwrap();
    let aln = Alignment::from_fasta_str(">A\nMKVLGA--\n>B\nMKVLGAWE\n", "fam").unwrap();
    let join = JoinOptions {
        length_slack: u64::MAX,
        ..Default::default()
    };
    let proj = join_annotations(&[tx("A", &[15, 18]), tx("B", &[3])], &ps, &join).unwrap();
    let opts = ProjectOptions {
        tail: TailPolicy::Clamp,
        ..Default::default()
    };
    let ap = project_into_alignment(&aln, &proj, &ps, &opts);
    assert_eq!(ap.taxa[0].report.placed, 1);
    assert_eq!(ap.taxa[0].report.clamped, 1);
    assert_eq!(ap.taxa[0].sites.iter().filter(|s| s.clamped).count(), 1);
}

#[test]
fn counts_report_missing_separately_for_unannotated_taxa() {
    use intron_core::matrix::{Bitset, Site};
    use intron_core::project::PhaseMode;
    use intron_core::score::score_counts;
    let (mut p, mut o) = (Bitset::new(4), Bitset::new(4));
    p.set(0);
    o.set(0);
    o.set(1);
    let m = SiteMatrix::from_parts(
        "x",
        ["A", "B", "C", "D"].iter().map(|s| s.to_string()).collect(),
        vec![Site {
            col: 2,
            phase: 0,
            presence: p,
            observable: o,
        }],
        PhaseMode::Strict,
    );
    let c = &score_counts(&m)[0];
    assert_eq!((c.n_present, c.n_absent, c.n_missing), (1, 1, 2));
    assert!(
        (c.frequency - 0.5).abs() < 1e-12,
        "missing cells are excluded from the frequency"
    );
}

#[test]
fn streaming_the_annotation_into_the_join_equals_collecting_it_first() {
    use intron_core::annotation::{AnnotationSource, GffSource};
    use intron_core::project::JoiningSink;
    let dir = std::env::temp_dir().join(format!("ip_stream_join_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let mut paths = Vec::new();
    let mut fasta = String::new();
    for f in 0..5u32 {
        let mut gff = String::new();
        for g in 0..30u32 {
            let id = format!("f{f}g{g}");
            let strand = if (f + g) % 2 == 0 { '+' } else { '-' };
            let len2 = 3 * (g % 5 + 1);
            for (a, b) in [(100u32, 114u32), (200, 199 + len2)] {
                gff += &format!("c{f}\tt\tCDS\t{a}\t{b}\t.\t{strand}\t0\tID={id}.c;Parent={id};protein_id=P{id}\n");
            }
            fasta += &format!(">P{id}\n{}\n", "M".repeat(5 + len2 as usize / 3));
        }
        let p = dir.join(format!("f{f}.gff3"));
        fs::write(&p, gff).unwrap();
        paths.push(p);
    }
    let faa = dir.join("p.faa");
    fs::write(&faa, fasta).unwrap();
    let proteins = ProteinSet::read([&faa]).unwrap();
    let opts = JoinOptions::default();

    let set = GffSource::new(paths.clone(), Default::default())
        .load()
        .unwrap();
    let collected = join_annotations(&set.transcripts, &proteins, &opts).unwrap();

    let mut sink = JoiningSink::new(&proteins, &opts);
    GffSource::new(paths, Default::default())
        .stream_into(&mut sink)
        .unwrap();
    let (streamed, ann_diag, n_tx) = sink.finish().unwrap();

    assert_eq!(n_tx, set.transcripts.len());
    assert_eq!(ann_diag.count_of("source_summary"), 1);
    assert_eq!(streamed.rows.len(), collected.rows.len());
    for (a, b) in streamed.rows.iter().zip(&collected.rows) {
        assert_eq!(streamed.protein_id(a), collected.protein_id(b));
        assert_eq!(streamed.transcript_id(a), collected.transcript_id(b));
        assert_eq!(
            (a.site, a.genomic_start, a.genomic_end, a.coding_rank),
            (b.site, b.genomic_start, b.genomic_end, b.coding_rank)
        );
    }
    assert_eq!(streamed.joined_proteins, collected.joined_proteins);
    assert_eq!(streamed.diagnostics.tally(), collected.diagnostics.tally());
}
