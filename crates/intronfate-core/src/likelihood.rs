//! Likelihood model of intron presence/absence that includes annotation error.
//!
//! # Model
//!
//! Every position (character) evolves independently along a rooted tree with two states, absent (0) and present (1).
//!
//! - **Root**: present with probability `root_present`.
//! - **Branch** leading to node `v`: from absent, the position is gained with probability `gain[v]`; from present, it is lost with
//!   probability `loss[v]`. Branch-specific transition probabilities, as in MALIN's model without rate variation (Csűrös 2008); gains
//!   and losses may both happen anywhere, so this is not Dollo.
//! - **Annotation error** at tip `t`: a true intron is annotated with probability `sensitivity[t]` (else it is a *miss*); at a position
//!   where the genome truly has no intron, an intron is annotated with probability `false_rate[t]` (a *false intron*). A `?` cell is
//!   missing and contributes nothing. The error rates are **inputs** (e.g. from a comparison with a reference annotation): for one
//!   genome a miss and a loss on its terminal branch produce the same data, so they cannot be estimated from the matrix itself.
//! - **Ascertainment**: a position exists in the matrix only because at least one genome shows an intron there. Each position's
//!   likelihood is therefore divided by `1 − Z`, where `Z` is the probability that, with the same missing cells, no genome shows one.
//!
//! # The root
//!
//! With two branches at the root, the root's presence probability and the two root branches are not separately identifiable
//! (five parameters for a joint distribution of the two children's states with three degrees of freedom), and fitting drifts
//! along that ridge. The branch from the root to its **larger** child (the ingroup) is therefore fixed to "no change": the model's
//! root state is the ingroup ancestor's state, `root_present` is the probability that the ingroup ancestor has the intron, and all
//! change between the outgroup and the ingroup ancestor is placed on the outgroup's branch ([`fixed_root_branch`]).
//!
//! # Fitting
//!
//! Expectation–maximisation. The positions that were never seen (all observed cells 0) are treated as missing data: for a position
//! with no-show probability `Z`, on average `Z / (1 − Z)` unseen positions share its pattern of missing cells, and their expected
//! transitions are added to the counts with that weight (the usual EM for a zero-truncated sample). The M-step is then closed form:
//! `gain = E[0→1] / E[parent 0]`, `loss = E[1→0] / E[parent 1]`, `root_present = E[root 1] / E[all]`. Identical columns are scored
//! once. Each iteration does not decrease the conditional log-likelihood.
//!
//! # Outputs
//!
//! Posterior expectations over the observed positions: gains and losses on each branch, the number of positions present at each
//! node (an ancestral-state estimate with uncertainty folded in), and per tip the expected misses and false introns; also the
//! conditional log-likelihood per position, for comparing trees.
//!
//! Not modelled (v1): rate variation across positions (MALIN's gamma categories), which makes ancestral counts and losses depend on
//! how heterogeneous the positions are; dependence between neighbouring positions.

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::score::{ScoringData, ABSENT, PRESENT};
use crate::tree::{NodeId, Tree};

/// Probabilities are kept inside `[CLAMP, 1 − CLAMP]` so that every transition stays possible.
const CLAMP: f64 = 1e-9;

/// Per-tip annotation error, in [`ScoringData::tips`] order.
#[derive(Clone, Debug, PartialEq)]
pub struct ErrorRates {
    /// P(annotated | the genome has the intron).
    pub sensitivity: Vec<f64>,
    /// P(annotated | the genome has no intron at this position).
    pub false_rate: Vec<f64>,
}

impl ErrorRates {
    /// No error: every intron annotated, nothing invented.
    pub fn none(n_tips: usize) -> Self {
        ErrorRates {
            sensitivity: vec![1.0; n_tips],
            false_rate: vec![0.0; n_tips],
        }
    }

    pub fn validate(&self, n_tips: usize) -> Result<(), String> {
        if self.sensitivity.len() != n_tips || self.false_rate.len() != n_tips {
            return Err(format!(
                "error rates for {} / {} tips, the tree has {n_tips}",
                self.sensitivity.len(),
                self.false_rate.len()
            ));
        }
        for (i, (&s, &f)) in self.sensitivity.iter().zip(&self.false_rate).enumerate() {
            if !(s > 0.0 && s <= 1.0) {
                return Err(format!(
                    "sensitivity of tip {i} is {s}; it must be in (0, 1]"
                ));
            }
            if !(0.0..1.0).contains(&f) || f >= s {
                return Err(format!("false-intron rate of tip {i} is {f}; it must be in [0, 1) and below the sensitivity ({s})"));
            }
        }
        Ok(())
    }
}

/// Per-position false-intron rate from annotation precision.
///
/// With `carried` positions annotated in a genome, of which a fraction `precision` are real, `carried·(1 − precision)` are false; the
/// genome truly lacks an intron at about `observed − carried·precision/sensitivity` of its `observed` (non-`?`) positions. The rate is
/// the first over the second.
pub fn false_rate_from_precision(
    carried: f64,
    observed: f64,
    sensitivity: f64,
    precision: f64,
) -> f64 {
    let false_ones = carried * (1.0 - precision);
    let truly_absent = (observed - carried * precision / sensitivity).max(false_ones.max(1.0));
    (false_ones / truly_absent).clamp(0.0, 0.5)
}

/// Model parameters. `gain` and `loss` are indexed by node id; the root's entries are unused.
#[derive(Clone, Debug, PartialEq)]
pub struct Params {
    pub gain: Vec<f64>,
    pub loss: Vec<f64>,
    pub root_present: f64,
}

impl Params {
    pub fn uniform(tree: &Tree, gain: f64, loss: f64, root_present: f64) -> Self {
        Params {
            gain: vec![gain; tree.len()],
            loss: vec![loss; tree.len()],
            root_present,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FitOptions {
    pub max_iter: usize,
    /// Stop when an iteration improves the log-likelihood by less than this.
    pub tol: f64,
    /// Starting values (default: gain 0.05, loss 0.2, root 0.5 everywhere).
    pub init: Option<Params>,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            max_iter: 10_000,
            tol: 1e-4,
            init: None,
        }
    }
}

/// Posterior expectations for one node and the branch above it.
#[derive(Clone, Debug, Default)]
pub struct NodeStats {
    pub gain_prob: f64,
    pub loss_prob: f64,
    /// Expected gains on the branch above the node, over the observed positions. For the root: positions present at the root.
    pub expected_gains: f64,
    /// Expected losses on the branch above the node.
    pub expected_losses: f64,
    /// Expected number of observed positions present at the node.
    pub expected_present: f64,
}

/// Posterior expectations for one tip.
#[derive(Clone, Debug, Default)]
pub struct TipStats {
    /// Positions observed absent (0) that the genome is expected to really have.
    pub expected_missed: f64,
    /// Positions observed present (1) that are expected to be false introns.
    pub expected_false: f64,
}

#[derive(Clone, Debug)]
pub struct FitResult {
    pub params: Params,
    /// Conditional log-likelihood (given that each position is seen at least once).
    pub log_lik: f64,
    pub iterations: usize,
    pub converged: bool,
    /// Indexed by node id.
    pub nodes: Vec<NodeStats>,
    /// In [`ScoringData::tips`] order.
    pub tips: Vec<TipStats>,
    /// Conditional log-likelihood of each character, in [`ScoringData::characters`] order.
    pub site_log_lik: Vec<f64>,
    /// Expected number of positions that exist but were never seen (all observed cells absent).
    pub unobserved_positions: f64,
}

/// Distinct columns with their multiplicity, and each character's column.
struct Patterns {
    columns: Vec<Vec<u8>>,
    counts: Vec<f64>,
    of_char: Vec<usize>,
}

fn compress(data: &ScoringData) -> Patterns {
    let mut index: FxHashMap<&[u8], usize> = FxHashMap::default();
    let mut p = Patterns {
        columns: Vec::new(),
        counts: Vec::new(),
        of_char: Vec::with_capacity(data.characters.len()),
    };
    for ch in &data.characters {
        let k = *index.entry(ch.states.as_slice()).or_insert_with(|| {
            p.columns.push(ch.states.clone());
            p.counts.push(0.0);
            p.columns.len() - 1
        });
        p.counts[k] += 1.0;
        p.of_char.push(k);
    }
    p
}

/// P(observation | true state) at a tip, as `[given absent, given present]`.
fn emission(obs: u8, s: f64, f: f64) -> [f64; 2] {
    match obs {
        PRESENT => [f, s],
        ABSENT => [1.0 - f, 1.0 - s],
        _ => [1.0, 1.0],
    }
}

fn transition(p: &Params, v: NodeId) -> [[f64; 2]; 2] {
    let (g, l) = (p.gain[v], p.loss[v]);
    [[1.0 - g, g], [l, 1.0 - l]]
}

/// Per-thread scratch space for one column.
struct Work {
    up: Vec<[f64; 2]>,
    msg: Vec<[f64; 2]>,
    down: Vec<[f64; 2]>,
}

impl Work {
    fn new(n: usize) -> Self {
        Work {
            up: vec![[0.0; 2]; n],
            msg: vec![[0.0; 2]; n],
            down: vec![[0.0; 2]; n],
        }
    }
}

/// Everything about the tree that the passes need.
struct Ctx<'a> {
    tree: &'a Tree,
    /// Tip index (into the column) of each node, if it is a scored tip.
    tip_of: Vec<Option<usize>>,
    err: &'a ErrorRates,
}

impl<'a> Ctx<'a> {
    fn new(data: &'a ScoringData, err: &'a ErrorRates) -> Self {
        let mut tip_of = vec![None; data.tree.len()];
        for (i, &t) in data.tips.iter().enumerate() {
            tip_of[t] = Some(i);
        }
        Ctx {
            tree: &data.tree,
            tip_of,
            err,
        }
    }

    fn emit(&self, v: NodeId, col: &[u8]) -> [f64; 2] {
        match self.tip_of[v] {
            Some(i) => emission(col[i], self.err.sensitivity[i], self.err.false_rate[i]),
            None => [1.0, 1.0],
        }
    }

    /// Felsenstein pruning with per-node rescaling. Returns ln P(column) (unconditional).
    fn upward(&self, p: &Params, col: &[u8], w: &mut Work) -> f64 {
        let mut log_scale = 0.0;
        for &v in self.tree.postorder() {
            let mut l = self.emit(v, col);
            for &c in self.tree.children(v) {
                l[0] *= w.msg[c][0];
                l[1] *= w.msg[c][1];
            }
            let mx = l[0].max(l[1]);
            if mx <= 0.0 {
                return f64::NEG_INFINITY;
            }
            l = [l[0] / mx, l[1] / mx];
            log_scale += mx.ln();
            w.up[v] = l;
            if self.tree.parent(v).is_some() {
                let t = transition(p, v);
                w.msg[v] = [
                    t[0][0] * l[0] + t[0][1] * l[1],
                    t[1][0] * l[0] + t[1][1] * l[1],
                ];
            }
        }
        let r = self.tree.root();
        let tot = (1.0 - p.root_present) * w.up[r][0] + p.root_present * w.up[r][1];
        log_scale + tot.ln()
    }

    /// Outside pass after [`Ctx::upward`] on the same column; adds `weight` × posterior expectations to `acc`.
    fn accumulate(
        &self,
        p: &Params,
        col: &[u8],
        w: &mut Work,
        weight: f64,
        real: bool,
        acc: &mut Acc,
    ) {
        let tree = self.tree;
        let r = tree.root();
        w.down[r] = [1.0 - p.root_present, p.root_present];
        for &v in tree.postorder().iter().rev() {
            let dv = w.down[v];
            let post = norm([dv[0] * w.up[v][0], dv[1] * w.up[v][1]]);
            if real {
                acc.present[v] += weight * post[1];
                if let Some(i) = self.tip_of[v] {
                    match col[i] {
                        ABSENT => acc.missed[i] += weight * post[1],
                        PRESENT => acc.false_[i] += weight * post[0],
                        _ => {}
                    }
                }
            }
            if v == r {
                acc.root[0] += weight * post[0];
                acc.root[1] += weight * post[1];
            }
            let kids = tree.children(v);
            for &c in kids {
                // everything at v except the subtree of c
                let mut ex = dv;
                for &s in kids {
                    if s != c {
                        ex[0] *= w.msg[s][0];
                        ex[1] *= w.msg[s][1];
                    }
                }
                let t = transition(p, c);
                w.down[c] = norm([
                    ex[0] * t[0][0] + ex[1] * t[1][0],
                    ex[0] * t[0][1] + ex[1] * t[1][1],
                ]);
                let uc = w.up[c];
                let mut pair = [
                    [ex[0] * t[0][0] * uc[0], ex[0] * t[0][1] * uc[1]],
                    [ex[1] * t[1][0] * uc[0], ex[1] * t[1][1] * uc[1]],
                ];
                let z: f64 = pair.iter().flatten().sum();
                for row in pair.iter_mut() {
                    for x in row.iter_mut() {
                        *x /= z;
                    }
                }
                for (a, row) in pair.iter().enumerate() {
                    for (b, &x) in row.iter().enumerate() {
                        acc.n[c][a][b] += weight * x;
                    }
                }
                if real {
                    acc.gains[c] += weight * pair[0][1];
                    acc.losses[c] += weight * pair[1][0];
                }
            }
        }
    }
}

fn norm(x: [f64; 2]) -> [f64; 2] {
    let s = x[0] + x[1];
    [x[0] / s, x[1] / s]
}

/// The column with every observed cell turned to "absent": the event "nobody shows an intron" for this pattern of missing cells.
fn phantom(col: &[u8]) -> Vec<u8> {
    col.iter()
        .map(|&c| if c == PRESENT { ABSENT } else { c })
        .collect()
}

#[derive(Clone)]
struct Acc {
    /// Expected transitions per branch, observed + unseen positions (for the M-step).
    n: Vec<[[f64; 2]; 2]>,
    root: [f64; 2],
    /// Observed positions only.
    gains: Vec<f64>,
    losses: Vec<f64>,
    present: Vec<f64>,
    missed: Vec<f64>,
    false_: Vec<f64>,
    log_lik: f64,
    unobserved: f64,
}

impl Acc {
    fn new(n: usize, tips: usize) -> Self {
        Acc {
            n: vec![[[0.0; 2]; 2]; n],
            root: [0.0; 2],
            gains: vec![0.0; n],
            losses: vec![0.0; n],
            present: vec![0.0; n],
            missed: vec![0.0; tips],
            false_: vec![0.0; tips],
            log_lik: 0.0,
            unobserved: 0.0,
        }
    }

    fn add(mut self, o: Acc) -> Acc {
        for (a, b) in self.n.iter_mut().zip(&o.n) {
            for i in 0..2 {
                for j in 0..2 {
                    a[i][j] += b[i][j];
                }
            }
        }
        for i in 0..2 {
            self.root[i] += o.root[i];
        }
        for (a, b) in [
            (&mut self.gains, &o.gains),
            (&mut self.losses, &o.losses),
            (&mut self.present, &o.present),
        ] {
            a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
        }
        self.missed
            .iter_mut()
            .zip(&o.missed)
            .for_each(|(x, y)| *x += y);
        self.false_
            .iter_mut()
            .zip(&o.false_)
            .for_each(|(x, y)| *x += y);
        self.log_lik += o.log_lik;
        self.unobserved += o.unobserved;
        self
    }
}

/// Conditional ln-likelihood of one column: ln P(column) − ln(1 − Z).
fn column_terms(ctx: &Ctx, p: &Params, col: &[u8], w: &mut Work) -> (f64, f64) {
    let ln_l = ctx.upward(p, col, w);
    let ph = phantom(col);
    let mut w2 = Work::new(ctx.tree.len());
    let z = ctx.upward(p, &ph, &mut w2).exp().min(1.0 - 1e-15);
    (ln_l - (1.0 - z).ln(), z)
}

/// One E-step over all distinct columns.
fn e_step(ctx: &Ctx, pats: &Patterns, p: &Params) -> Acc {
    let n = ctx.tree.len();
    let tips = ctx.err.sensitivity.len();
    (0..pats.columns.len())
        .into_par_iter()
        .fold(
            || (Acc::new(n, tips), Work::new(n), Work::new(n)),
            |(mut acc, mut w, mut wp), k| {
                let col = &pats.columns[k];
                let cnt = pats.counts[k];
                let ln_l = ctx.upward(p, col, &mut w);
                ctx.accumulate(p, col, &mut w, cnt, true, &mut acc);
                let ph = phantom(col);
                let z = ctx.upward(p, &ph, &mut wp).exp().min(1.0 - 1e-15);
                let wt = cnt * z / (1.0 - z);
                ctx.accumulate(p, &ph, &mut wp, wt, false, &mut acc);
                acc.log_lik += cnt * (ln_l - (1.0 - z).ln());
                acc.unobserved += wt;
                (acc, w, wp)
            },
        )
        .map(|(acc, _, _)| acc)
        .reduce(|| Acc::new(n, tips), Acc::add)
}

fn m_step(tree: &Tree, acc: &Acc, p: &mut Params, fixed: Option<NodeId>) {
    for v in 0..tree.len() {
        if tree.parent(v).is_none() || Some(v) == fixed {
            continue;
        }
        let c = acc.n[v];
        let from0 = c[0][0] + c[0][1];
        let from1 = c[1][0] + c[1][1];
        if from0 > 0.0 {
            p.gain[v] = (c[0][1] / from0).clamp(CLAMP, 1.0 - CLAMP);
        }
        if from1 > 0.0 {
            p.loss[v] = (c[1][0] / from1).clamp(CLAMP, 1.0 - CLAMP);
        }
    }
    let t = acc.root[0] + acc.root[1];
    if t > 0.0 {
        p.root_present = (acc.root[1] / t).clamp(CLAMP, 1.0 - CLAMP);
    }
}

/// The root branch that is fixed to "no change": the root's larger child (by tips; ties to the first), when the root has exactly
/// two children. `None` for a root with three or more children, where the parameters are identifiable as they are.
pub fn fixed_root_branch(tree: &Tree) -> Option<NodeId> {
    let kids = tree.children(tree.root());
    if kids.len() != 2 {
        return None;
    }
    let n = |v: NodeId| tree.leaves_under(v).len();
    Some(if n(kids[1]) > n(kids[0]) {
        kids[1]
    } else {
        kids[0]
    })
}

/// Conditional log-likelihood of the data under fixed parameters (no fitting).
pub fn log_likelihood(data: &ScoringData, err: &ErrorRates, p: &Params) -> f64 {
    let ctx = Ctx::new(data, err);
    let mut w = Work::new(data.tree.len());
    data.characters
        .iter()
        .map(|ch| column_terms(&ctx, p, &ch.states, &mut w).0)
        .sum()
}

/// Fit the model by EM and return the parameters with the posterior expectations at the optimum.
pub fn fit(data: &ScoringData, err: &ErrorRates, opt: &FitOptions) -> Result<FitResult, String> {
    err.validate(data.tips.len())?;
    if data.characters.is_empty() {
        return Err("no characters to fit".into());
    }
    let tree = &data.tree;
    let ctx = Ctx::new(data, err);
    let pats = compress(data);
    let mut p = opt
        .init
        .clone()
        .unwrap_or_else(|| Params::uniform(tree, 0.05, 0.2, 0.5));
    let fixed = fixed_root_branch(tree);
    if let Some(f) = fixed {
        p.gain[f] = CLAMP;
        p.loss[f] = CLAMP;
    }
    let mut prev = f64::NEG_INFINITY;
    let mut iterations = 0;
    let mut converged = false;
    let mut acc = e_step(&ctx, &pats, &p);
    while iterations < opt.max_iter {
        if acc.log_lik - prev < opt.tol && iterations > 0 {
            converged = true;
            break;
        }
        prev = acc.log_lik;
        m_step(tree, &acc, &mut p, fixed);
        acc = e_step(&ctx, &pats, &p);
        iterations += 1;
    }
    let nodes = (0..tree.len())
        .map(|v| {
            let root = tree.parent(v).is_none();
            NodeStats {
                gain_prob: if root { p.root_present } else { p.gain[v] },
                loss_prob: if root { 0.0 } else { p.loss[v] },
                expected_gains: if root { acc.present[v] } else { acc.gains[v] },
                expected_losses: if root { 0.0 } else { acc.losses[v] },
                expected_present: acc.present[v],
            }
        })
        .collect();
    let tips = (0..data.tips.len())
        .map(|i| TipStats {
            expected_missed: acc.missed[i],
            expected_false: acc.false_[i],
        })
        .collect();
    let mut w = Work::new(tree.len());
    let per_pattern: Vec<f64> = pats
        .columns
        .iter()
        .map(|c| column_terms(&ctx, &p, c, &mut w).0)
        .collect();
    let site_log_lik = pats.of_char.iter().map(|&k| per_pattern[k]).collect();
    Ok(FitResult {
        params: p,
        log_lik: acc.log_lik,
        iterations,
        converged,
        nodes,
        tips,
        site_log_lik,
        unobserved_positions: acc.unobserved,
    })
}

/// Percentile intervals over bootstrap replicates, per node.
#[derive(Clone, Debug, Default)]
pub struct Bootstrap {
    pub replicates: usize,
    /// Resampling unit: "family" (all positions of an alignment together) or "position".
    pub unit: &'static str,
    /// (2.5 %, 97.5 %) per node id.
    pub gains: Vec<(f64, f64)>,
    pub losses: Vec<(f64, f64)>,
    pub present: Vec<(f64, f64)>,
}

/// The family of a character label `<family>:<column>:<phase>` (combined matrices), or `None` for `<column>:<phase>`.
fn family_of(label: &str) -> Option<&str> {
    let mut it = label.rsplitn(3, ':');
    let (_, _, fam) = (it.next(), it.next(), it.next());
    fam
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let i = (q * (sorted.len() - 1) as f64).round() as usize;
    sorted[i.min(sorted.len() - 1)]
}

// Node ids index several parallel arrays here.
#[allow(clippy::needless_range_loop)]
/// Nonparametric bootstrap of the expected gains, losses and present counts.
///
/// Resamples whole families (all positions of one alignment) when the labels carry one and there are at least two, since positions
/// in one gene are not independent; otherwise single positions. Each replicate is refitted by EM starting from `start` (the fit to
/// the full data), with the same error rates.
pub fn bootstrap(
    data: &ScoringData,
    err: &ErrorRates,
    opt: &FitOptions,
    start: &Params,
    replicates: usize,
    seed: u64,
) -> Result<Bootstrap, String> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut index: FxHashMap<&str, usize> = FxHashMap::default();
    for (i, ch) in data.characters.iter().enumerate() {
        if let Some(f) = family_of(&ch.label) {
            let k = *index.entry(f).or_insert_with(|| {
                groups.push(Vec::new());
                groups.len() - 1
            });
            groups[k].push(i);
        }
    }
    let unit = if groups.len() >= 2
        && groups.iter().map(Vec::len).sum::<usize>() == data.characters.len()
    {
        "family"
    } else {
        groups = (0..data.characters.len()).map(|i| vec![i]).collect();
        "position"
    };
    let n = data.tree.len();
    let mut draws: Vec<[Vec<f64>; 3]> = vec![[Vec::new(), Vec::new(), Vec::new()]; n];
    let mut state = seed ^ 0x5DEECE66D;
    let mut next = || {
        state = state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    };
    for _ in 0..replicates {
        let mut chars = Vec::with_capacity(data.characters.len());
        for _ in 0..groups.len() {
            let g = &groups[(next() % groups.len() as u64) as usize];
            chars.extend(g.iter().map(|&i| data.characters[i].clone()));
        }
        let rep = ScoringData {
            characters: chars,
            ..data.clone()
        };
        let r = fit(
            &rep,
            err,
            &FitOptions {
                init: Some(start.clone()),
                ..opt.clone()
            },
        )?;
        for v in 0..n {
            draws[v][0].push(r.nodes[v].expected_gains);
            draws[v][1].push(r.nodes[v].expected_losses);
            draws[v][2].push(r.nodes[v].expected_present);
        }
    }
    let ci = |x: &mut Vec<f64>| {
        x.sort_by(|a, b| a.total_cmp(b));
        (percentile(x, 0.025), percentile(x, 0.975))
    };
    let mut b = Bootstrap {
        replicates,
        unit,
        ..Default::default()
    };
    for d in draws.iter_mut() {
        b.gains.push(ci(&mut d[0]));
        b.losses.push(ci(&mut d[1]));
        b.present.push(ci(&mut d[2]));
    }
    Ok(b)
}

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;
    use crate::score::{Character, MISSING};

    fn tree(nwk: &str) -> Tree {
        Tree::from_newick(nwk).unwrap()
    }

    /// ScoringData straight from columns given in tip-label order.
    fn data(nwk: &str, cols: &[&[u8]]) -> ScoringData {
        let t = tree(nwk);
        let tips = t.tips();
        let tip_names: Vec<String> = tips
            .iter()
            .map(|&v| t.label(v).unwrap().to_string())
            .collect();
        let characters = cols
            .iter()
            .enumerate()
            .map(|(i, c)| Character {
                label: format!("{i}:0"),
                col: i as u32,
                phase: 0,
                states: c.to_vec(),
            })
            .collect();
        ScoringData {
            annotated_tip: vec![true; tips.len()],
            tree: t,
            tips,
            tip_names,
            characters,
        }
    }

    /// Sum over every assignment of true states to all nodes: (P(column), P(no observed intron), posteriors).
    fn brute(
        d: &ScoringData,
        err: &ErrorRates,
        p: &Params,
        col: &[u8],
    ) -> (f64, f64, Vec<f64>, Vec<[[f64; 2]; 2]>) {
        let t = &d.tree;
        let n = t.len();
        let mut tip_of = vec![None; n];
        for (i, &v) in d.tips.iter().enumerate() {
            tip_of[v] = Some(i);
        }
        let mut lik = 0.0;
        let mut z = 0.0;
        let mut present = vec![0.0; n];
        let mut pairs = vec![[[0.0; 2]; 2]; n];
        for mask in 0u32..(1 << n) {
            let st = |v: usize| ((mask >> v) & 1) as usize;
            let mut prob = if st(t.root()) == 1 {
                p.root_present
            } else {
                1.0 - p.root_present
            };
            for v in 0..n {
                if let Some(u) = t.parent(v) {
                    prob *= transition(p, v)[st(u)][st(v)];
                }
            }
            let mut e = 1.0;
            let mut e0 = 1.0;
            for v in 0..n {
                if let Some(i) = tip_of[v] {
                    e *= emission(col[i], err.sensitivity[i], err.false_rate[i])[st(v)];
                    let ph = if col[i] == MISSING { MISSING } else { ABSENT };
                    e0 *= emission(ph, err.sensitivity[i], err.false_rate[i])[st(v)];
                }
            }
            lik += prob * e;
            z += prob * e0;
            for v in 0..n {
                present[v] += prob * e * st(v) as f64;
                if let Some(u) = t.parent(v) {
                    pairs[v][st(u)][st(v)] += prob * e;
                }
            }
        }
        for x in present.iter_mut() {
            *x /= lik;
        }
        for pr in pairs.iter_mut() {
            for row in pr.iter_mut() {
                for x in row.iter_mut() {
                    *x /= lik;
                }
            }
        }
        (lik, z, present, pairs)
    }

    fn some_params(t: &Tree) -> Params {
        let mut p = Params::uniform(t, 0.1, 0.3, 0.6);
        for v in 0..t.len() {
            p.gain[v] = 0.03 + 0.04 * v as f64;
            p.loss[v] = 0.5 - 0.035 * v as f64;
        }
        p
    }

    const NWK: &str = "((A,B),(C,(D,E)));";

    fn some_errors() -> ErrorRates {
        ErrorRates {
            sensitivity: vec![0.9, 1.0, 0.7, 0.95, 0.8],
            false_rate: vec![0.01, 0.0, 0.05, 0.02, 0.0],
        }
    }

    #[test]
    fn pruning_and_posteriors_match_brute_force() {
        let cols: &[&[u8]] = &[
            &[1, 0, 1, MISSING, 0],
            &[0, 0, 0, 1, 1],
            &[1, 1, 1, 1, 1],
            &[MISSING, MISSING, 1, 0, MISSING],
        ];
        let d = data(NWK, cols);
        let err = some_errors();
        let p = some_params(&d.tree);
        let ctx = Ctx::new(&d, &err);
        for col in cols {
            let (lik, z, present, pairs) = brute(&d, &err, &p, col);
            let mut w = Work::new(d.tree.len());
            let ln = ctx.upward(&p, col, &mut w);
            assert!((ln - lik.ln()).abs() < 1e-10, "{ln} vs {}", lik.ln());
            let (cond, z2) = column_terms(&ctx, &p, col, &mut w);
            assert!((z - z2).abs() < 1e-12 && (cond - (lik.ln() - (1.0 - z).ln())).abs() < 1e-10);
            let mut acc = Acc::new(d.tree.len(), 5);
            ctx.upward(&p, col, &mut w);
            ctx.accumulate(&p, col, &mut w, 1.0, true, &mut acc);
            for v in 0..d.tree.len() {
                assert!((acc.present[v] - present[v]).abs() < 1e-10, "node {v}");
                if d.tree.parent(v).is_some() {
                    for a in 0..2 {
                        for b in 0..2 {
                            assert!(
                                (acc.n[v][a][b] - pairs[v][a][b]).abs() < 1e-10,
                                "edge {v} {a}{b}"
                            );
                        }
                    }
                    assert!(
                        (acc.gains[v] - pairs[v][0][1]).abs() < 1e-10
                            && (acc.losses[v] - pairs[v][1][0]).abs() < 1e-10
                    );
                }
            }
        }
    }

    #[test]
    fn conditional_probabilities_of_all_seen_patterns_sum_to_one() {
        let d = data(NWK, &[&[1, 0, 0, 0, 0]]);
        let err = some_errors();
        let p = some_params(&d.tree);
        let ctx = Ctx::new(&d, &err);
        let mut w = Work::new(d.tree.len());
        let mut total = 0.0;
        for mask in 1u32..32 {
            let col: Vec<u8> = (0..5).map(|i| ((mask >> i) & 1) as u8).collect();
            total += column_terms(&ctx, &p, &col, &mut w).0.exp();
        }
        assert!((total - 1.0).abs() < 1e-12, "{total}");
        // with a missing cell, the patterns over the other four cells sum to one as well
        let mut total = 0.0;
        for mask in 1u32..16 {
            let mut col: Vec<u8> = (0..4).map(|i| ((mask >> i) & 1) as u8).collect();
            col.insert(2, MISSING);
            total += column_terms(&ctx, &p, &col, &mut w).0.exp();
        }
        assert!((total - 1.0).abs() < 1e-12, "{total}");
    }

    /// splitmix64
    struct Rng(u64);
    impl Rng {
        fn f(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Simulate `n` columns that show at least one intron (unseen ones are dropped, as in real data), with 5 % of cells missing.
    fn simulate(t: &Tree, p: &Params, err: &ErrorRates, n: usize, seed: u64) -> Vec<Vec<u8>> {
        simulate_counted(t, p, err, n, seed).0
    }

    /// As [`simulate`], also returning the true number of gains and losses on each branch among the kept columns.
    fn simulate_counted(
        t: &Tree,
        p: &Params,
        err: &ErrorRates,
        n: usize,
        seed: u64,
    ) -> (Vec<Vec<u8>>, Vec<f64>, Vec<f64>) {
        let tips = t.tips();
        let mut rng = Rng(seed);
        let mut out = Vec::new();
        let (mut gains, mut losses) = (vec![0.0; t.len()], vec![0.0; t.len()]);
        let order: Vec<NodeId> = t.postorder().iter().rev().copied().collect();
        while out.len() < n {
            let mut st = vec![0u8; t.len()];
            for &v in &order {
                st[v] = match t.parent(v) {
                    None => (rng.f() < p.root_present) as u8,
                    Some(u) if st[u] == 0 => (rng.f() < p.gain[v]) as u8,
                    Some(_) => (rng.f() >= p.loss[v]) as u8,
                };
            }
            let col: Vec<u8> = tips
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    if rng.f() < 0.05 {
                        MISSING
                    } else if st[v] == 1 {
                        (rng.f() < err.sensitivity[i]) as u8
                    } else {
                        (rng.f() < err.false_rate[i]) as u8
                    }
                })
                .collect();
            if col.contains(&PRESENT) {
                out.push(col);
                for v in 0..t.len() {
                    if let Some(u) = t.parent(v) {
                        gains[v] += (st[u] == 0 && st[v] == 1) as u8 as f64;
                        losses[v] += (st[u] == 1 && st[v] == 0) as u8 as f64;
                    }
                }
            }
        }
        (out, gains, losses)
    }

    fn sim_data(nwk: &str, cols: Vec<Vec<u8>>) -> ScoringData {
        let refs: Vec<&[u8]> = cols.iter().map(|c| c.as_slice()).collect();
        data(nwk, &refs)
    }

    #[test]
    fn em_recovers_simulated_parameters_and_never_decreases_the_likelihood() {
        let t = tree(NWK);
        let mut truth = Params::uniform(&t, 0.05, 0.15, 0.7);
        let c = t.tip_index("C").unwrap();
        truth.loss[c] = 0.35;
        truth.gain[c] = 0.12;
        let f = fixed_root_branch(&t).unwrap(); // (C,(D,E)): the larger root child
        truth.gain[f] = CLAMP;
        truth.loss[f] = CLAMP;
        let err = some_errors();
        let (cols, true_gains, true_losses) = simulate_counted(&t, &truth, &err, 60_000, 7);
        let d = sim_data(NWK, cols);
        // the log-likelihood never decreases from one iteration to the next
        let mut last = f64::NEG_INFINITY;
        for it in [1, 2, 5, 10, 20, 40] {
            let r = fit(
                &d,
                &err,
                &FitOptions {
                    max_iter: it,
                    tol: 0.0,
                    init: None,
                },
            )
            .unwrap();
            assert!(
                r.log_lik >= last - 1e-6,
                "iteration {it}: {} < {last}",
                r.log_lik
            );
            last = r.log_lik;
        }
        let r = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!(r.converged);
        // What is reported: expected gains and losses per branch, against the events that really happened in the simulation.
        // Losses are well determined. With annotation error in the model, the number of gains is not: it trades off against the
        // number of never-seen positions and the root probability. Over three seeds of 60,000 positions the fitted/true total
        // gains were 1.24, 1.00 and 0.67 (losses 0.97-1.04), narrowing to 0.97-1.04 at 1.5 million positions. So gains get a
        // wide tolerance here and need a bootstrap interval in practice.
        for v in 0..t.len() {
            if t.parent(v).is_none() {
                continue;
            }
            let (eg, el) = (r.nodes[v].expected_gains, r.nodes[v].expected_losses);
            assert!(
                (eg - true_gains[v]).abs() < 0.6 * true_gains[v] + 100.0,
                "gains {v}: {eg:.0} vs {}",
                true_gains[v]
            );
            assert!(
                (el - true_losses[v]).abs() < 0.10 * true_losses[v] + 300.0,
                "losses {v}: {el:.0} vs {}",
                true_losses[v]
            );
            assert!(
                (r.params.loss[v] - truth.loss[v]).abs() < 0.05,
                "loss {v}: {} vs {}",
                r.params.loss[v],
                truth.loss[v]
            );
            assert!(
                (r.params.gain[v] - truth.gain[v]).abs() < 0.04,
                "gain {v}: {} vs {}",
                r.params.gain[v],
                truth.gain[v]
            );
        }
        assert!((r.params.root_present - truth.root_present).abs() < 0.05);
        let (sg, tg): (f64, f64) = (
            r.nodes.iter().skip(1).map(|n| n.expected_gains).sum(),
            true_gains.iter().sum(),
        );
        let (sl, tl): (f64, f64) = (
            r.nodes.iter().skip(1).map(|n| n.expected_losses).sum(),
            true_losses.iter().sum(),
        );
        assert!((sl / tl - 1.0).abs() < 0.05, "total losses {sl:.0} vs {tl}");
        assert!((sg / tg - 1.0).abs() < 0.35, "total gains {sg:.0} vs {tg}");
        // the fitted model is at least as likely as the truth on these data
        assert!(r.log_lik >= log_likelihood(&d, &err, &truth) - 1e-6);
        // per-site log-likelihoods add up to the total
        assert!((r.site_log_lik.iter().sum::<f64>() - r.log_lik).abs() < 1e-6 * r.log_lik.abs());
    }

    #[test]
    fn ignoring_misses_inflates_terminal_losses() {
        let t = tree(NWK);
        let mut truth = Params::uniform(&t, 0.05, 0.1, 0.7);
        let f = fixed_root_branch(&t).unwrap();
        truth.gain[f] = CLAMP;
        truth.loss[f] = CLAMP;
        let c = t.tip_index("C").unwrap();
        let ci = t.tips().iter().position(|&v| v == c).unwrap();
        let err = some_errors(); // C: sensitivity 0.7
        let d = sim_data(NWK, simulate(&t, &truth, &err, 60_000, 11));
        let with = fit(&d, &err, &FitOptions::default()).unwrap();
        let without = fit(&d, &ErrorRates::none(5), &FitOptions::default()).unwrap();
        assert!(
            (with.params.loss[c] - 0.1).abs() < 0.04,
            "with error model: {}",
            with.params.loss[c]
        );
        // a miss rate of 0.3 looks like ~30 % extra terminal loss when errors are ignored
        assert!(
            without.params.loss[c] > 0.3,
            "without error model: {}",
            without.params.loss[c]
        );
        assert!(
            with.tips[ci].expected_missed > 0.0 && without.tips[ci].expected_missed.abs() < 1e-9
        );
    }

    #[test]
    fn the_larger_root_child_is_fixed() {
        let t = tree(NWK);
        let f = fixed_root_branch(&t).unwrap();
        assert_eq!(t.leaves_under(f).len(), 3);
        assert_eq!(fixed_root_branch(&tree("(A,B,(C,D));")), None);
        let d = sim_data(
            NWK,
            simulate(
                &t,
                &Params::uniform(&t, 0.05, 0.1, 0.7),
                &some_errors(),
                2000,
                3,
            ),
        );
        let r = fit(&d, &some_errors(), &FitOptions::default()).unwrap();
        assert_eq!((r.params.gain[f], r.params.loss[f]), (CLAMP, CLAMP));
        assert!(r.nodes[f].expected_gains < 1e-3 && r.nodes[f].expected_losses < 1e-3);
    }

    #[test]
    fn bootstrap_resamples_families_and_brackets_the_estimate() {
        assert_eq!(family_of("OG1:12:0"), Some("OG1"));
        assert_eq!(family_of("a:b:12:0"), Some("a:b"));
        assert_eq!(family_of("12:0"), None);
        let t = tree(NWK);
        let mut truth = Params::uniform(&t, 0.05, 0.15, 0.7);
        let f = fixed_root_branch(&t).unwrap();
        truth.gain[f] = CLAMP;
        truth.loss[f] = CLAMP;
        let err = some_errors();
        let mut d = sim_data(NWK, simulate(&t, &truth, &err, 3000, 5));
        for (i, ch) in d.characters.iter_mut().enumerate() {
            ch.label = format!("fam{}:{}:0", i / 10, i); // 300 families of 10 positions
        }
        let r = fit(&d, &err, &FitOptions::default()).unwrap();
        let b = bootstrap(&d, &err, &FitOptions::default(), &r.params, 20, 1).unwrap();
        assert_eq!((b.unit, b.replicates), ("family", 20));
        let c = t.tip_index("C").unwrap();
        let (lo, hi) = b.losses[c];
        assert!(
            lo < hi
                && lo <= r.nodes[c].expected_losses * 1.02
                && hi >= r.nodes[c].expected_losses * 0.98,
            "{lo} {hi} {}",
            r.nodes[c].expected_losses
        );
        // without family labels, positions are resampled
        let d2 = sim_data(NWK, simulate(&t, &truth, &err, 500, 5));
        assert_eq!(
            bootstrap(&d2, &err, &FitOptions::default(), &r.params, 2, 1)
                .unwrap()
                .unit,
            "position"
        );
    }

    #[test]
    fn identical_columns_are_scored_once_with_their_multiplicity() {
        let a: &[u8] = &[1, 0, 1, 0, 0];
        let b: &[u8] = &[0, 1, 1, 1, 0];
        let d1 = data(NWK, &[a, a, a, b]);
        let err = some_errors();
        let p = some_params(&d1.tree);
        let pats = compress(&d1);
        assert_eq!(
            (
                pats.columns.len(),
                pats.counts.clone(),
                pats.of_char.clone()
            ),
            (2, vec![3.0, 1.0], vec![0, 0, 0, 1])
        );
        let ctx = Ctx::new(&d1, &err);
        let acc = e_step(&ctx, &pats, &p);
        assert!((acc.log_lik - log_likelihood(&d1, &err, &p)).abs() < 1e-10);
    }

    #[test]
    fn error_rates_are_checked_and_precision_converts() {
        assert!(ErrorRates::none(3).validate(3).is_ok());
        assert!(ErrorRates {
            sensitivity: vec![0.0],
            false_rate: vec![0.0]
        }
        .validate(1)
        .is_err());
        assert!(ErrorRates {
            sensitivity: vec![0.5],
            false_rate: vec![0.6]
        }
        .validate(1)
        .is_err());
        assert!(ErrorRates::none(2).validate(3).is_err());
        // 1,000 annotated, 90 % real, sensitivity 0.9 -> 100 false among ~ 10,000 - 1,000 truly absent positions
        let f = false_rate_from_precision(1000.0, 10_000.0, 0.9, 0.9);
        assert!((f - 100.0 / 9000.0).abs() < 1e-12, "{f}");
        assert_eq!(false_rate_from_precision(1000.0, 10_000.0, 0.9, 1.0), 0.0);
    }

    #[test]
    fn polytomies_and_unscored_tips_are_handled() {
        // a root with three children and a tip (E) that is in the tree but has no data
        let t = tree("(A,B,(C,(D,E)));");
        let tips = t.tips();
        let keep: Vec<NodeId> = tips
            .iter()
            .copied()
            .filter(|&v| t.label(v) != Some("E"))
            .collect();
        let names = keep
            .iter()
            .map(|&v| t.label(v).unwrap().to_string())
            .collect();
        let d = ScoringData {
            annotated_tip: vec![true; 4],
            tips: keep,
            tip_names: names,
            characters: vec![Character {
                label: "1:0".into(),
                col: 1,
                phase: 0,
                states: vec![1, 0, 1, 1],
            }],
            tree: t.clone(),
        };
        let err = ErrorRates {
            sensitivity: vec![0.9, 0.8, 1.0, 0.95],
            false_rate: vec![0.0, 0.01, 0.02, 0.0],
        };
        let p = some_params(&t);
        // E counts as missing: the same as a 5-tip column with E = '?'
        let full = data("(A,B,(C,(D,E)));", &[&[1, 0, 1, 1, MISSING]]);
        let err5 = ErrorRates {
            sensitivity: vec![0.9, 0.8, 1.0, 0.95, 1.0],
            false_rate: vec![0.0, 0.01, 0.02, 0.0, 0.0],
        };
        assert!((log_likelihood(&d, &err, &p) - log_likelihood(&full, &err5, &p)).abs() < 1e-12);
        let (lik, z, _, _) = brute(&full, &err5, &p, &[1, 0, 1, 1, MISSING]);
        assert!((log_likelihood(&full, &err5, &p) - (lik.ln() - (1.0 - z).ln())).abs() < 1e-10);
    }
}
