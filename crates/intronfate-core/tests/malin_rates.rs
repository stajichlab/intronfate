//! Quantitative checks of the likelihood model against MALIN's own estimates (ignored by default; the data are MALIN's licensed test
//! set, so run it where you have them, for example on the HPCC):
//!
//!   INTRONFATE_MALIN_TABLE=<dir>/rogozin_introns_table INTRONFATE_MALIN_RATES=<dir>/Rog-ecdysozoa.rates \
//!   cargo test --release -p intronfate-core --test malin_rates -- --ignored --nocapture
//!
//! `Rog-ecdysozoa.rates` is the output of MALIN's `calculateConstantRates` (one rate class); `Rog-ecdysozoa1000.rates` and `Chromo23.rates` are
//! `calculateRASRates` outputs (classes with their own probabilities and loss/gain multipliers). For each, the test (1) fits this tool's model
//! with the same classes (free levels, SPEC §5.5d) to the same table and tree, (2) evaluates this tool's log-likelihood at MALIN's printed
//! parameters, and (3) prints the two parameter sets side by side. If both tools maximise the same likelihood, (1) is at least (2) and the
//! difference is small. Optional: `INTRONFATE_MALIN_MIN_OBSERVED=N` keeps positions observed in at least N taxa (MALIN's table filter),
//! `INTRONFATE_MALIN_LL=x` is MALIN's own reported log-likelihood, printed for comparison, and
//! `INTRONFATE_MALIN_COUNT_ONLY=1` only prints how many positions each minimum-observed filter keeps, and
//! `INTRONFATE_MALIN_NEVER_SEEN=global` models the never-seen positions as MALIN does (one probability p0 that a position has no intron
//! at any genome; SPEC §5.5) instead of conditioning each position on its own pattern of missing cells.
use intronfate_core::likelihood::{
    fit, log_likelihood_with, no_intron_probability, ErrorRates, FitOptions, NeverSeen, Params,
    Rates,
};
use intronfate_core::malin::{from_malin, parse_rates};
use intronfate_core::matrix::SiteMatrix;
use intronfate_core::score::{prepare_with, PrepareOptions};
use intronfate_core::tree::Tree;

fn run(table: &str, rates: &str) {
    let ms = from_malin(&std::fs::read_to_string(table).unwrap(), "malin").unwrap();
    let m = if ms.len() == 1 {
        ms.into_iter().next().unwrap()
    } else {
        SiteMatrix::combine(&ms, "malin").unwrap()
    };
    let mr = parse_rates(&std::fs::read_to_string(rates).unwrap()).unwrap();
    let (gain_lv, loss_lv) = if mr.classes.len() > 1 {
        mr.levels().unwrap()
    } else {
        (Vec::new(), Vec::new())
    };
    let tree = Tree::from_newick(&mr.tree).unwrap();
    let min_observed: usize = std::env::var("INTRONFATE_MALIN_MIN_OBSERVED")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    if std::env::var("INTRONFATE_MALIN_COUNT_ONLY").is_ok() {
        for n in 0..=tree.n_tips() {
            let (d, _) = prepare_with(
                &m,
                &tree,
                &[],
                &PrepareOptions {
                    min_taxa: 1,
                    min_observed: n,
                },
            )
            .unwrap();
            println!("min observed taxa {n}: {} positions", d.characters.len());
        }
        return;
    }
    let (d, report) = prepare_with(
        &m,
        &tree,
        &[],
        &PrepareOptions {
            min_taxa: 1,
            min_observed,
        },
    )
    .unwrap();
    println!(
        "{} sites in the table, {} scored (min observed taxa {min_observed}): {report:?}",
        m.sites.len(),
        d.characters.len()
    );
    let err = ErrorRates::none(d.tips.len());
    let never_seen = match std::env::var("INTRONFATE_MALIN_NEVER_SEEN").as_deref() {
        Ok("global") => NeverSeen::Global,
        _ => NeverSeen::PerPattern,
    };
    println!("never-seen positions: {never_seen:?}");
    let opt = FitOptions {
        never_seen,
        tol: 1e-9,
        param_tol: 1e-8,
        max_iter: 500_000,
        free_gain_levels: gain_lv.len().max(1),
        free_loss_levels: loss_lv.len().max(1),
        ..Default::default()
    };
    let ours = fit(&d, &err, &opt).unwrap();
    // with the larger root branch left free, as in MALIN (it is only weakly identifiable under rate classes, SPEC §5.5)
    let open = FitOptions {
        free_root_branch: true,
        ..opt.clone()
    };
    // MALIN's parameters in this tool's Params (by node name); the node probabilities refer to the reference levels
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
    pm.rates = Rates::free(gain_lv.clone(), loss_lv.clone());
    let l_malin = log_likelihood_with(&d, &err, &pm, never_seen);
    let (p0, n) = (
        no_intron_probability(&d, &err, &pm),
        d.characters.len() as f64,
    );
    println!(
        "at MALIN's parameters: p0 {p0:.6}, predicted never-seen positions {:.2}, uncorrected lnL {:.4}",
        n * p0 / (1.0 - p0),
        l_malin + n * (1.0 - p0).ln()
    );
    println!(
        "lnL at this tool's maximum {:.4} ({} E-steps, converged {})   lnL at MALIN's parameters {:.4}   difference {:+.4}",
        ours.log_lik, ours.iterations, ours.converged, l_malin, ours.log_lik - l_malin
    );
    if let Ok(ll) = std::env::var("INTRONFATE_MALIN_LL") {
        println!("MALIN's own reported log-likelihood: {ll}");
    }
    println!(
        "root presence: ours {:.4}  MALIN {:.4}",
        ours.params.root_present, mr.classes[0].root
    );
    println!("MALIN loss levels (prob, mult): {loss_lv:?}   gain levels: {gain_lv:?}");
    println!(
        "ours  loss levels: {:?}   gain levels: {:?}",
        ours.params.rates.free_loss, ours.params.rates.free_gain
    );
    println!(
        "{:<14} {:>9} {:>9} {:>9} {:>9}   fixed",
        "branch", "gain", "gain(M)", "loss", "loss(M)"
    );
    for v in 0..d.tree.len() {
        if d.tree.parent(v).is_none() {
            continue;
        }
        let name = d.tree.label(v).unwrap();
        println!(
            "{:<14} {:>9.5} {:>9.5} {:>9.5} {:>9.5}   {:?}",
            name,
            ours.params.gain[v],
            pm.gain[v],
            ours.params.loss[v],
            pm.loss[v],
            ours.nodes[v].fixed
        );
    }
    // a fit with the root branch free, started at MALIN's own parameters, cannot end below the value at those parameters
    let relaxed = fit(
        &d,
        &err,
        &FitOptions {
            init: Some(pm.clone()),
            ..open
        },
    )
    .unwrap();
    println!(
        "root branch free, started at MALIN's parameters: lnL {:.4} ({} E-steps, converged {}); difference from lnL at MALIN's parameters {:+.4}",
        relaxed.log_lik,
        relaxed.iterations,
        relaxed.converged,
        relaxed.log_lik - l_malin
    );
    println!(
        "  loss levels {:?}  gain levels {:?}  root {:.4}",
        relaxed.params.rates.free_loss, relaxed.params.rates.free_gain, relaxed.params.root_present
    );
    assert!(
        relaxed.log_lik >= l_malin - 1e-6,
        "the fit started at MALIN's parameters ended below them"
    );
    // with the root rule (the default) the maximum can only be lower, by the cost of the restriction
    assert!(ours.log_lik <= relaxed.log_lik + 1e-3);
}

#[test]
#[ignore]
fn likelihood_fit_matches_malins_estimates() {
    let (Ok(table), Ok(rates)) = (
        std::env::var("INTRONFATE_MALIN_TABLE"),
        std::env::var("INTRONFATE_MALIN_RATES"),
    ) else {
        eprintln!("set INTRONFATE_MALIN_TABLE and INTRONFATE_MALIN_RATES");
        return;
    };
    run(&table, &rates);
}
