//! Convergence check of the likelihood fit on a real matrix (ignored by default; run on a cluster):
//!   INTRONFATE_ML_MATRIX=matrix.tsv INTRONFATE_ML_TREE=tree.nwk cargo test --release -p intronfate-core --test ml_realdata -- --ignored --nocapture
//! Fits with the default options from several starting points, and once with a far stricter tolerance, under no error and under a
//! uniform error (sensitivity 0.95, false rate 0.001); prints lnL, E-steps and the total expected gains and losses. A sound fit gives
//! the same optimum from every start.
use intronfate_core::likelihood::{fit, ErrorRates, FitOptions, Params};
use intronfate_core::matrix::SiteMatrix;
use intronfate_core::score::{prepare_with, PrepareOptions};
use intronfate_core::tree::Tree;

#[test]
#[ignore]
fn likelihood_converges_to_the_same_optimum_on_real_data() {
    let (Ok(m), Ok(t)) = (
        std::env::var("INTRONFATE_ML_MATRIX"),
        std::env::var("INTRONFATE_ML_TREE"),
    ) else {
        eprintln!("set INTRONFATE_ML_MATRIX and INTRONFATE_ML_TREE");
        return;
    };
    let m = SiteMatrix::read_tsv(&m).unwrap();
    let t = Tree::read(&t).unwrap();
    let (d, _) = prepare_with(
        &m,
        &t,
        &[],
        &PrepareOptions {
            min_taxa: 1,
            min_observed: 0,
        },
    )
    .unwrap();
    let n = d.tips.len();
    for (name, err) in [
        ("none", ErrorRates::none(n)),
        (
            "uniform",
            ErrorRates {
                sensitivity: vec![0.95; n],
                false_rate: vec![0.001; n],
            },
        ),
    ] {
        let mut lnls = Vec::new();
        for (label, init, tol, max) in [
            ("default", None, 1e-6, 50_000),
            (
                "start2",
                Some(Params::uniform(&d.tree, 0.2, 0.05, 0.8)),
                1e-6,
                50_000,
            ),
            (
                "start3",
                Some(Params::uniform(&d.tree, 0.01, 0.4, 0.3)),
                1e-6,
                50_000,
            ),
            ("strict", None, 1e-10, 500_000),
        ] {
            let r = fit(
                &d,
                &err,
                &FitOptions {
                    max_iter: max,
                    tol,
                    init,
                },
            )
            .unwrap();
            let g: f64 = r
                .nodes
                .iter()
                .skip(1)
                .filter(|x| x.has_data)
                .map(|x| x.expected_gains)
                .sum();
            let l: f64 = r
                .nodes
                .iter()
                .filter(|x| x.has_data)
                .map(|x| x.expected_losses)
                .sum();
            eprintln!(
                "{name:8} {label:8} lnL {:.4} E-steps {:6} converged {} gains {g:.1} losses {l:.1} unseen {:.1} root {:.4}",
                r.log_lik, r.iterations, r.converged, r.unobserved_positions, r.params.root_present
            );
            lnls.push(r.log_lik);
        }
        let spread = lnls.iter().cloned().fold(f64::MIN, f64::max)
            - lnls.iter().cloned().fold(f64::MAX, f64::min);
        eprintln!("{name}: lnL spread over starts {spread:.4}");
    }
}
