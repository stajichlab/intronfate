//! Quantitative check of the likelihood model against MALIN's own estimates (ignored by default; the data are MALIN's licensed test set,
//! so run it where you have them, for example on the HPCC):
//!   INTRONFATE_MALIN_TABLE=<dir>/rogozin_introns_table INTRONFATE_MALIN_RATES=<dir>/Rog-ecdysozoa.rates \
//!   cargo test --release -p intronfate-core --test malin_rates -- --ignored --nocapture
//! `Rog-ecdysozoa.rates` is the output of MALIN's `calculateConstantRates` (one rate class, no rate variation) on that table and tree.
//! The test (1) fits this tool's model (`-m likelihood`, no annotation error, one class) to the same table and tree, (2) evaluates this
//! tool's log-likelihood at MALIN's parameters, and (3) prints the two parameter sets side by side. If both tools maximise the same
//! likelihood, (1) is at least (2) and the difference is small; the printed branch probabilities agree where the branch is identifiable
//! (the two root branches and the root probability are not separately identifiable, SPEC §5.5).
use intronfate_core::likelihood::{fit, log_likelihood, ErrorRates, FitOptions, Params};
use intronfate_core::malin::{from_malin, parse_rates};
use intronfate_core::score::{prepare_with, PrepareOptions};
use intronfate_core::tree::Tree;

#[test]
#[ignore]
fn likelihood_fit_matches_malins_constant_rates() {
    let (Ok(table), Ok(rates)) = (
        std::env::var("INTRONFATE_MALIN_TABLE"),
        std::env::var("INTRONFATE_MALIN_RATES"),
    ) else {
        eprintln!("set INTRONFATE_MALIN_TABLE and INTRONFATE_MALIN_RATES");
        return;
    };
    let ms = from_malin(&std::fs::read_to_string(&table).unwrap(), "malin").unwrap();
    assert_eq!(
        ms.len(),
        1,
        "expected one matrix from a table without #MAP lines"
    );
    let mr = parse_rates(&std::fs::read_to_string(&rates).unwrap()).unwrap();
    assert_eq!(
        mr.classes.len(),
        1,
        "this check is for a one-class (constant rates) file"
    );
    let tree = Tree::from_newick(&mr.tree).unwrap();
    let (d, _) = prepare_with(
        &ms[0],
        &tree,
        &[],
        &PrepareOptions {
            min_taxa: 1,
            min_observed: 0,
        },
    )
    .unwrap();
    let err = ErrorRates::none(d.tips.len());
    let ours = fit(
        &d,
        &err,
        &FitOptions {
            tol: 1e-10,
            param_tol: 1e-9,
            max_iter: 500_000,
            ..Default::default()
        },
    )
    .unwrap();
    // MALIN's parameters in this tool's Params (by node name)
    let mut pm = Params::uniform(&d.tree, 0.5, 0.5, mr.classes[0].root);
    for v in 0..d.tree.len() {
        if d.tree.parent(v).is_none() {
            continue;
        }
        let name = d
            .tree
            .label(v)
            .expect("every node is named in a MALIN .rates tree");
        let n = mr
            .nodes
            .iter()
            .find(|n| n.name == name)
            .unwrap_or_else(|| panic!("node {name} not in the .rates file"));
        pm.gain[v] = n.p_gain;
        pm.loss[v] = n.p_loss;
    }
    let l_malin = log_likelihood(&d, &err, &pm);
    println!(
        "sites {}   lnL at this tool's MLE {:.4} ({} E-steps, converged {})   lnL at MALIN's parameters {:.4}   difference {:+.4}",
        d.characters.len(),
        ours.log_lik,
        ours.iterations,
        ours.converged,
        l_malin,
        ours.log_lik - l_malin
    );
    println!(
        "root presence: ours {:.4}  MALIN {:.4}",
        ours.params.root_present, mr.classes[0].root
    );
    println!(
        "{:<12} {:>9} {:>9} {:>9} {:>9}   fixed",
        "branch", "gain", "gain(M)", "loss", "loss(M)"
    );
    for v in 0..d.tree.len() {
        if d.tree.parent(v).is_none() {
            continue;
        }
        let name = d.tree.label(v).unwrap();
        println!(
            "{:<12} {:>9.5} {:>9.5} {:>9.5} {:>9.5}   {:?}",
            name,
            ours.params.gain[v],
            pm.gain[v],
            ours.params.loss[v],
            pm.loss[v],
            ours.nodes[v].fixed
        );
    }
    assert!(ours.converged);
    // the maximum cannot be below the value at MALIN's parameters if the two tools maximise the same function
    assert!(
        ours.log_lik >= l_malin - 1e-6,
        "this tool's maximum is below its value at MALIN's parameters"
    );
}
