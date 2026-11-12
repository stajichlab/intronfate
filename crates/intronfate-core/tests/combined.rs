//! A combined matrix must score exactly like the sum of its families, survive a file round trip, and keep sites of different
//! alignments apart (same column in two families is two sites, not one).

use intronfate_core::matrix::{Bitset, FamilyBlock, Site, SiteMatrix, FAMILY_PAD};
use intronfate_core::project::PhaseMode;
use intronfate_core::score::{prepare, score_dollo};
use intronfate_core::tables::write_matrix_tsv;
use intronfate_core::tree::Tree;

/// A matrix from `(column, phase, carriers, observed)` over `taxa`; `observed` lists the taxa that were annotated and aligned.
fn matrix(name: &str, taxa: &[&str], sites: &[(u32, u8, &[&str], &[&str])]) -> SiteMatrix {
    let idx = |n: &str| taxa.iter().position(|t| *t == n).unwrap();
    let s: Vec<Site> = sites
        .iter()
        .map(|(col, phase, carriers, observed)| {
            let mut presence = Bitset::new(taxa.len());
            let mut observable = Bitset::new(taxa.len());
            for c in *carriers {
                presence.set(idx(c));
                observable.set(idx(c));
            }
            for o in *observed {
                observable.set(idx(o));
            }
            Site {
                col: *col,
                phase: *phase,
                presence,
                observable,
            }
        })
        .collect();
    SiteMatrix::from_parts(
        name,
        taxa.iter().map(|t| t.to_string()).collect(),
        s,
        PhaseMode::Strict,
    )
}

fn families() -> (SiteMatrix, SiteMatrix) {
    // family f1 has taxa A B C; f2 has only A B D (C absent from it, D absent from f1)
    let f1 = matrix(
        "f1",
        &["A", "B", "C"],
        &[
            (5, 0, &["A", "B"], &["C"]),
            (9, 1, &["A", "C"], &["B"]),
            (20, 2, &["C"], &["A", "B"]),
        ],
    );
    let f2 = matrix(
        "f2",
        &["A", "B", "D"],
        &[(5, 0, &["A", "D"], &["B"]), (7, 0, &["B", "D"], &["A"])],
    );
    (f1, f2)
}

fn tree() -> Tree {
    Tree::from_newick("(((A,B),C),D);").unwrap()
}

#[test]
fn combined_scores_equal_the_sum_of_the_families() {
    let (f1, f2) = families();
    let c = SiteMatrix::combine(&[f1.clone(), f2.clone()], "combined").unwrap();
    assert_eq!(c.taxa, ["A", "B", "C", "D"]);
    assert_eq!(c.n_sites(), 5);
    let t = tree();
    let score = |m: &SiteMatrix| {
        let (d, _) = prepare(m, &t, &[], 1).unwrap();
        let r = score_dollo(&d, false);
        (r.total_gains, r.total_losses)
    };
    let (a, b, both) = (score(&f1), score(&f2), score(&c));
    assert_eq!(both, (a.0 + b.0, a.1 + b.1), "{a:?} {b:?} {both:?}");
    // per branch too
    let branches = |m: &SiteMatrix| {
        let (d, _) = prepare(m, &t, &[], 1).unwrap();
        let mut v: Vec<(String, u64, u64)> = score_dollo(&d, false)
            .branches
            .iter()
            .map(|b| (b.tips.join(","), b.gains, b.losses))
            .collect();
        v.sort();
        v
    };
    let mut summed: std::collections::BTreeMap<String, (u64, u64)> = Default::default();
    for m in [&f1, &f2] {
        for (k, g, l) in branches(m) {
            let e = summed.entry(k).or_default();
            e.0 += g;
            e.1 += l;
        }
    }
    for (k, g, l) in branches(&c) {
        let want = summed.get(&k).copied().unwrap_or_default();
        assert_eq!((g, l), want, "branch {k}");
    }
}

#[test]
fn missing_taxa_stay_missing_and_sites_keep_their_family() {
    let (f1, f2) = families();
    let c = SiteMatrix::combine(&[f1, f2], "combined").unwrap();
    let d = c.taxa.iter().position(|t| t == "D").unwrap();
    // D is absent from f1: its cells at f1's three sites are `?`, not `0`
    for i in 0..3 {
        assert_eq!(c.sites[i].cell(d).symbol(), '?');
    }
    // C is absent from f2
    let cc = c.taxa.iter().position(|t| t == "C").unwrap();
    for i in 3..5 {
        assert_eq!(c.sites[i].cell(cc).symbol(), '?');
    }
    assert_eq!(c.site_label(0), "f1:5:0");
    assert_eq!(
        c.site_label(3),
        "f2:5:0",
        "same column in f2 is a different site"
    );
    assert_eq!(c.site_column(4), 7);
    assert_eq!(
        c.families,
        vec![
            FamilyBlock {
                name: "f1".into(),
                offset: 0
            },
            FamilyBlock {
                name: "f2".into(),
                offset: 21 + FAMILY_PAD
            },
        ]
    );
    // columns are unique and ascending
    assert!(c
        .sites
        .windows(2)
        .all(|w| (w[0].col, w[0].phase) < (w[1].col, w[1].phase)));
    // pairwise sharing is per family: A and B share f1's (5,0) only; f2's (5,0) is carried by A and D
    let pairs = c.pairwise(0, true);
    let ab = pairs
        .iter()
        .find(|p| p.taxon_a == "A" && p.taxon_b == "B")
        .unwrap();
    assert_eq!(ab.shared, 1);
    let ad = pairs
        .iter()
        .find(|p| p.taxon_a == "A" && p.taxon_b == "D")
        .unwrap();
    assert_eq!(ad.shared, 1, "only f2's site, not matched across families");
}

#[test]
fn file_round_trip_keeps_families_and_scores() {
    let (f1, f2) = families();
    let c = SiteMatrix::combine(&[f1, f2], "combined").unwrap();
    let dir = std::env::temp_dir().join(format!("if_combined_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("matrix.tsv");
    write_matrix_tsv(&c, &path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("families=2 family_pad=10000"), "{text}");
    assert!(text
        .lines()
        .nth(1)
        .unwrap()
        .starts_with("site\tcolumn\tphase\tfamily\tfamily_column\t"));
    let r = SiteMatrix::read_tsv(&path).unwrap();
    assert_eq!(r.families, c.families);
    assert_eq!(r.site_family, c.site_family);
    assert_eq!(r.n_sites(), c.n_sites());
    for i in 0..c.n_sites() {
        assert_eq!(r.site_label(i), c.site_label(i));
        for t in 0..c.taxa.len() {
            assert_eq!(r.sites[i].cell(t), c.sites[i].cell(t));
        }
    }
    let t = tree();
    let tot = |m: &SiteMatrix| {
        let (d, _) = prepare(m, &t, &[], 1).unwrap();
        let s = score_dollo(&d, false);
        (s.total_gains, s.total_losses)
    };
    assert_eq!(tot(&r), tot(&c));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn combining_refuses_ambiguous_inputs() {
    let (f1, f2) = families();
    let mut twin = f2.clone();
    twin.alignment = "f1".into();
    assert!(SiteMatrix::combine(&[f1.clone(), twin], "x")
        .unwrap_err()
        .to_string()
        .contains("both named"));
    let mut other = f2.clone();
    other.flank = 3;
    assert!(SiteMatrix::combine(&[f1.clone(), other], "x")
        .unwrap_err()
        .to_string()
        .contains("phase mode or flank"));
    let c = SiteMatrix::combine(&[f1.clone(), f2], "c").unwrap();
    assert!(SiteMatrix::combine(&[c, f1], "x")
        .unwrap_err()
        .to_string()
        .contains("already a combined"));
    assert!(SiteMatrix::combine(&[], "x").is_err());
}
