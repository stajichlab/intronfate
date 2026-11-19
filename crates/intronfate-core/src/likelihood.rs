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
//! # Rate variation across positions
//!
//! Positions differ in how readily they are gained and lost. With [`FitOptions::categories`] `K > 1`, each position falls into one
//! of `K` equally likely rate classes whose multipliers `r_k` are the category means of a gamma distribution with mean 1 and shape
//! `α` (Yang 1994's discrete gamma; MALIN's rate variation, Csűrös & Miklós 2009, uses one such factor on both rates). For that
//! each branch is read as a two-state continuous-time process: its rate-1 probabilities `(g, l)` define the equilibrium frequency
//! `π = g / (g + l)` and the decay `d = 1 − g − l`, and in class `k` the branch has `g_k = π (1 − d^{r_k})`,
//! `l_k = (1 − π)(1 − d^{r_k})`. This needs `g + l < 1` on every branch, which the fit keeps. A position's likelihood is the mean
//! over the classes, and so is the probability `Z` that nobody shows it. `α` is estimated (or fixed) along with the branches; small
//! `α` means strong heterogeneity (many nearly frozen positions and a few fast ones). With `K = 1` the model is exactly the one
//! above, and `g + l` is not restricted.
//!
//! EM with classes: the E-step gives expected transition counts per branch and class, and the M-step maximises the expected
//! complete-data log-likelihood by coordinate ascent (the shape by a one-dimensional search, then each branch's `(π, d)`), accepting
//! only moves that raise it, which keeps the log-likelihood from decreasing (generalised EM).
//!
//! # Separate rate factors for gain and loss
//!
//! With [`FitOptions::gain_categories`] and [`FitOptions::loss_categories`] the gain rate and the loss rate of a position are scaled by
//! different factors, each a mean-1 discrete gamma with its own shape (MALIN's rate variation has separate loss and gain levels). The
//! position classes are the grid of gain classes × loss classes, equally likely, with the loss class as the outer index. A branch's
//! rate-1 probabilities `(g, l)` give `π = g / (g + l)` and `d = 1 − g − l` as above; with factors `(r_g, r_l)` the gain rate is
//! `π·(−ln d)·r_g` and the loss rate `(1 − π)·(−ln d)·r_l`, so the class has `w = π·r_g + (1 − π)·r_l`, `π' = π·r_g / w` and
//! `g' = π'(1 − d^w)`, `l' = (1 − π')(1 − d^w)` (the transition matrix of that two-state process over the branch). Equal factors give
//! the shared form above. Shared classes and separate classes are not combined, and there are at most 16 classes in all.
//!
//! # Free levels (MALIN-style classes)
//!
//! With [`FitOptions::free_gain_levels`] and [`FitOptions::free_loss_levels`] each level of the gain (loss) dimension has its own
//! probability and rate multiplier instead of being a gamma slice. The first level is the reference (multiplier 1; the branch
//! probabilities refer to it), a position class is a pair of levels and its probability the product of the two, as in MALIN's `|class`
//! lines. The EM updates the level probabilities in closed form from the expected class counts, each other multiplier by a
//! one-dimensional search, and the branches as before. Because the likelihood is multimodal in the levels, [`fit`] starts from several
//! spreads of the multipliers (or from the levels carried by `FitOptions::init`) and keeps the best. With rate classes the larger root
//! branch is only weakly identifiable; it stays fixed unless [`FitOptions::free_root_branch`] is set.
//!
//! # Never-seen positions: per pattern or MALIN's global correction
//!
//! By default each position is conditioned on being seen given its own pattern of missing cells ([`NeverSeen::PerPattern`]). MALIN
//! instead corrects the whole table with one probability `p0` that a position has no intron at any genome, all genomes being
//! observable ([`NeverSeen::Global`]): corrected log-likelihood = uncorrected − `N ln(1 − p0)`, and `N p0 / (1 − p0)` never-seen positions
//! (MALIN's "predicted missing"). The two agree when no cell is missing; with `?` cells they differ, and only the global form reproduces
//! MALIN's printed log-likelihood (SPEC §6b).
//!
//! Not modelled: dependence between neighbouring positions.

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::score::{presence_probability, ScoringData, ABSENT, MISSING, PRESENT};
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

/// Rate variation across positions. Either one factor shared by gain and loss (`categories` equally likely classes from a mean-1 gamma
/// with this `shape`), or separate factors: `gain_categories` classes for the gain rate and `loss_categories` for the loss rate, each its
/// own mean-1 gamma, combined as a grid of equally likely classes (as MALIN's rate variation has separate loss and gain levels). The
/// shared and the separate form are not combined: at most one of them has more than one class.
#[derive(Clone, Debug, PartialEq)]
pub struct Rates {
    /// Shared classes; 1 = none.
    pub categories: usize,
    /// Gamma shape `α` of the shared factor (unused with one shared class).
    pub shape: f64,
    /// Classes of the gain-rate factor; 1 = none.
    pub gain_categories: usize,
    /// Gamma shape of the gain-rate factor.
    pub gain_shape: f64,
    /// Classes of the loss-rate factor; 1 = none.
    pub loss_categories: usize,
    /// Gamma shape of the loss-rate factor.
    pub loss_shape: f64,
    /// Free gain levels `(probability, multiplier)` (MALIN's gain levels): empty or one level = none. The first level is the reference
    /// (multiplier 1, the one the branch probabilities refer to); the probabilities sum to 1. Not combined with the gamma classes.
    pub free_gain: Vec<(f64, f64)>,
    /// Free loss levels, as `free_gain`.
    pub free_loss: Vec<(f64, f64)>,
}

impl Default for Rates {
    fn default() -> Self {
        Rates {
            categories: 1,
            shape: 1.0,
            gain_categories: 1,
            gain_shape: 1.0,
            loss_categories: 1,
            loss_shape: 1.0,
            free_gain: Vec::new(),
            free_loss: Vec::new(),
        }
    }
}

impl Rates {
    /// Shared factor with `categories` classes of shape `shape`.
    pub fn shared(categories: usize, shape: f64) -> Self {
        Rates {
            categories,
            shape,
            ..Default::default()
        }
    }

    /// Separate gain and loss factors.
    pub fn separate(
        gain_categories: usize,
        gain_shape: f64,
        loss_categories: usize,
        loss_shape: f64,
    ) -> Self {
        Rates {
            gain_categories,
            gain_shape,
            loss_categories,
            loss_shape,
            ..Default::default()
        }
    }

    /// Free levels: `gain` and `loss` are `(probability, multiplier)` lists, the first entry of each the reference level.
    pub fn free(gain: Vec<(f64, f64)>, loss: Vec<(f64, f64)>) -> Self {
        Rates {
            free_gain: gain,
            free_loss: loss,
            ..Default::default()
        }
    }

    /// `true` when gain and loss have separate gamma factors.
    pub fn is_separate(&self) -> bool {
        self.gain_categories > 1 || self.loss_categories > 1
    }

    /// `true` when the classes are free levels (own probabilities and multipliers).
    pub fn is_free(&self) -> bool {
        self.free_gain.len() > 1 || self.free_loss.len() > 1
    }

    /// Number of position classes.
    pub fn n_classes(&self) -> usize {
        if self.is_free() {
            self.free_gain.len().max(1) * self.free_loss.len().max(1)
        } else if self.is_separate() {
            self.gain_categories.max(1) * self.loss_categories.max(1)
        } else {
            self.categories.max(1)
        }
    }
}

/// Bounds of the gamma shape.
pub const SHAPE_MIN: f64 = 0.05;
pub const SHAPE_MAX: f64 = 200.0;
/// Most rate categories accepted.
pub const MAX_CATEGORIES: usize = 16;

/// Model parameters. `gain` and `loss` are indexed by node id; the root's entries are unused. With rate classes they are the
/// probabilities at rate 1 (the mean rate).
#[derive(Clone, Debug, PartialEq)]
pub struct Params {
    pub gain: Vec<f64>,
    pub loss: Vec<f64>,
    pub root_present: f64,
    pub rates: Rates,
}

impl Params {
    /// The same probabilities on every branch, no rate variation.
    pub fn uniform(tree: &Tree, gain: f64, loss: f64, root_present: f64) -> Self {
        Params {
            gain: vec![gain; tree.len()],
            loss: vec![loss; tree.len()],
            root_present,
            rates: Rates::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct FitOptions {
    /// Exact limit on E-steps (at least 1). Near the limit the fit switches from accelerated to plain EM steps so the count never exceeds it.
    pub max_iter: usize,
    /// Stop when an acceleration cycle improves the log-likelihood by less than this ...
    pub tol: f64,
    /// ... and moves no parameter by more than this.
    pub param_tol: f64,
    /// Starting values (default: gain 0.05, loss 0.2, root 0.5 everywhere, shape 1). Their `rates.categories` is ignored: the fit
    /// uses [`FitOptions::categories`].
    pub init: Option<Params>,
    /// Rate classes across positions (1 = none).
    pub categories: usize,
    /// Fix the gamma shape at this value instead of estimating it.
    pub shape: Option<f64>,
    /// Separate rate classes for the gain rate (1 = none); not combined with `categories > 1`.
    pub gain_categories: usize,
    /// Fix the gamma shape of the gain-rate factor instead of estimating it.
    pub gain_shape: Option<f64>,
    /// Separate rate classes for the loss rate (1 = none).
    pub loss_categories: usize,
    /// Fix the gamma shape of the loss-rate factor instead of estimating it.
    pub loss_shape: Option<f64>,
    /// Free gain levels (MALIN-style: each level has its own probability and multiplier, the first is the reference with multiplier 1);
    /// 1 = none. Not combined with the gamma classes.
    pub free_gain_levels: usize,
    /// Free loss levels (see `free_gain_levels`).
    pub free_loss_levels: usize,
    /// Do not fix the larger root branch (see [`fixed_branches_with`]). Only meaningful with rate classes; slows the fit.
    pub free_root_branch: bool,
    /// How never-seen positions are modelled (see [`NeverSeen`]).
    pub never_seen: NeverSeen,
    /// Print one line per acceleration cycle to stderr (log-likelihood, shape, time in E- and M-steps).
    pub trace: bool,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            max_iter: 50_000,
            tol: 1e-6,
            param_tol: 1e-7,
            init: None,
            categories: 1,
            shape: None,
            gain_categories: 1,
            gain_shape: None,
            loss_categories: 1,
            loss_shape: None,
            free_gain_levels: 1,
            free_loss_levels: 1,
            free_root_branch: false,
            never_seen: NeverSeen::PerPattern,
            trace: false,
        }
    }
}

/// ln Γ(x) for x > 0 (Lanczos, g = 7, n = 9; relative error about 1e-15).
fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // reflection
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, &c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Regularised lower incomplete gamma P(a, x) (series below a + 1, continued fraction above; Numerical Recipes 6.2).
fn gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let lead = a * x.ln() - x - ln_gamma(a);
    if x < a + 1.0 {
        let (mut ap, mut sum, mut del) = (a, 1.0 / a, 1.0 / a);
        for _ in 0..10_000 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-16 {
                break;
            }
        }
        (sum.ln() + lead).exp().min(1.0)
    } else {
        let tiny = 1e-300;
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..10_000 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + an / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-16 {
                break;
            }
        }
        (1.0 - (h.ln() + lead).exp()).max(0.0)
    }
}

/// The `q` quantile of Gamma(shape a, rate 1), by bisection on ln x.
fn gamma_quantile(a: f64, q: f64) -> f64 {
    let (mut lo, mut hi) = (-700.0f64, 12.0f64); // e^12 ≈ 1.6e5, far above any quantile used here (a ≤ 200 + 1)
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if gamma_p(a, mid.exp()) < q {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-13 {
            break;
        }
    }
    (0.5 * (lo + hi)).exp()
}

/// The rate multipliers of the classes: the means of `categories` equal-probability slices of a gamma with mean 1 and shape `shape`
/// (Yang 1994), rescaled to average exactly 1. `[1.0]` for one class.
pub fn category_rates(rates: &Rates) -> Vec<f64> {
    let k = rates.categories.max(1);
    if k == 1 {
        return vec![1.0];
    }
    let a = rates.shape;
    // with Gamma(a, rate a) for the rates, P(rate < c/a) = P(a, c) and the partial mean up to c/a is P(a + 1, c)
    let cuts: Vec<f64> = (1..k)
        .map(|i| gamma_quantile(a, i as f64 / k as f64))
        .collect();
    let mut prev = 0.0;
    let mut r: Vec<f64> = (0..k)
        .map(|i| {
            let next = if i + 1 < k {
                gamma_p(a + 1.0, cuts[i])
            } else {
                1.0
            };
            let m = (k as f64 * (next - prev)).max(1e-10);
            prev = next;
            m
        })
        .collect();
    let mean = r.iter().sum::<f64>() / k as f64;
    r.iter_mut().for_each(|x| *x /= mean);
    r
}

/// Transition matrix of a branch with rate-1 probabilities `(g, l)` at rate multiplier `r` (see the module documentation).
/// `r == 1` returns `(g, l)` themselves.
#[cfg(test)]
fn transition_at(g: f64, l: f64, r: f64) -> [[f64; 2]; 2] {
    transition_gl(g, l, r, r)
}

/// Transition matrix of a branch with rate-1 probabilities `(g, l)` when the gain rate is multiplied by `rg` and the loss rate by `rl`.
///
/// The branch is a two-state continuous-time process with gain rate `a` and loss rate `b` over a unit time, so `π = a/(a+b)` and
/// `d = e^{−(a+b)} = 1 − g − l`. With factors, `a' = a·rg` and `b' = b·rl`: `a + b = −ln d` is split by `π`, hence
/// `a' + b' = −ln d · (π·rg + (1 − π)·rl)`, `π' = π·rg / (π·rg + (1 − π)·rl)`, and the branch has `g' = π'(1 − d^{π·rg + (1−π)·rl})`,
/// `l' = (1 − π')(1 − d^{π·rg + (1−π)·rl})`. With `rg = rl = r` this is the shared-factor form `g' = π(1 − d^r)`.
fn transition_gl(g: f64, l: f64, rg: f64, rl: f64) -> [[f64; 2]; 2] {
    if rg == 1.0 && rl == 1.0 {
        return [[1.0 - g, g], [l, 1.0 - l]];
    }
    let s = (g + l).min(1.0 - CLAMP);
    let pi = g / (g + l);
    let w = pi * rg + (1.0 - pi) * rl; // the factor on the total rate
    let e = -(w * (-s).ln_1p()).exp_m1(); // 1 − (1 − s)^w
    let pi2 = pi * rg / w;
    let (gr, lr) = (pi2 * e, (1.0 - pi2) * e);
    [[1.0 - gr, gr], [lr, 1.0 - lr]]
}

/// A position class: its probability and the factors on the gain and the loss rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Class {
    pub w: f64,
    pub rg: f64,
    pub rl: f64,
}

/// The position classes with their probabilities. Gamma classes (shared or separate) are equally likely; free levels have their own
/// probabilities, and a class's probability is the product of its gain level's and its loss level's (loss outer, gain inner).
pub fn class_list(rates: &Rates) -> Vec<Class> {
    if rates.is_free() {
        let one = vec![(1.0, 1.0)];
        let gain = if rates.free_gain.len() > 1 {
            &rates.free_gain
        } else {
            &one
        };
        let loss = if rates.free_loss.len() > 1 {
            &rates.free_loss
        } else {
            &one
        };
        loss.iter()
            .flat_map(|&(pl, ml)| {
                gain.iter().map(move |&(pg, mg)| Class {
                    w: pl * pg,
                    rg: mg,
                    rl: ml,
                })
            })
            .collect()
    } else {
        let f = class_factors(rates);
        let w = 1.0 / f.len() as f64;
        f.into_iter().map(|(rg, rl)| Class { w, rg, rl }).collect()
    }
}

/// The (gain factor, loss factor) of every position class, equally likely. One shared factor gives `(r_k, r_k)`; separate factors give
/// the grid over loss classes (outer) and gain classes (inner); no variation gives `[(1, 1)]`.
pub fn class_factors(rates: &Rates) -> Vec<(f64, f64)> {
    if rates.is_free() {
        class_list(rates).iter().map(|c| (c.rg, c.rl)).collect()
    } else if rates.is_separate() {
        let gain = category_rates(&Rates::shared(
            rates.gain_categories.max(1),
            rates.gain_shape,
        ));
        let loss = category_rates(&Rates::shared(
            rates.loss_categories.max(1),
            rates.loss_shape,
        ));
        loss.iter()
            .flat_map(|&rl| gain.iter().map(move |&rg| (rg, rl)))
            .collect()
    } else {
        category_rates(rates).into_iter().map(|r| (r, r)).collect()
    }
}

/// One rate class: the transition matrix of every branch (indexed by node id) and the root probability.
struct Cat {
    tr: Vec<[[f64; 2]; 2]>,
    root: f64,
    /// The class probability.
    w: f64,
}

fn cats(p: &Params) -> Vec<Cat> {
    class_list(&p.rates)
        .iter()
        .map(|c| Cat {
            tr: (0..p.gain.len())
                .map(|v| transition_gl(p.gain[v], p.loss[v], c.rg, c.rl))
                .collect(),
            root: p.root_present,
            w: c.w,
        })
        .collect()
}

/// Why a branch's probabilities are not estimated (fixed at "no change").
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Fixed {
    /// Estimated.
    #[default]
    No,
    /// One branch at the root is fixed so that the root probability is identifiable (see the module documentation).
    Root,
    /// No observed tip below: the data say nothing about this branch.
    NoData,
    /// The node below has a single child with data, so only the product of this branch and that child's branch is identified;
    /// all change is put on the child's branch.
    Unary,
}

impl Fixed {
    pub fn as_str(self) -> &'static str {
        match self {
            Fixed::No => "no",
            Fixed::Root => "root",
            Fixed::NoData => "no_data",
            Fixed::Unary => "unary",
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
    /// `false` when no tip below the node has any observed cell: the branch is fixed at "no change" and its numbers are not
    /// estimates.
    pub has_data: bool,
    /// Whether, and why, the branch above the node is fixed at "no change".
    pub fixed: Fixed,
}

/// Posterior expectations for one tip.
#[derive(Clone, Debug, Default)]
pub struct TipStats {
    /// Positions observed absent (0) that the genome is expected to really have.
    pub expected_missed: f64,
    /// Positions observed present (1) that are expected to be false introns.
    pub expected_false: f64,
    /// Expected never-seen positions at which this genome is observable (not `?`) and truly lacks the intron: where a false intron
    /// could have appeared.
    pub unseen_observable: f64,
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
    /// Posterior mean rate multiplier of each character (1 without rate classes): above 1, a position that changes more readily
    /// than average. With separate gain and loss factors it is the mean of the two posterior means below.
    pub site_rate: Vec<f64>,
    /// Posterior mean multiplier of the gain rate of each character (equal to `site_rate` unless gain and loss have separate factors).
    pub site_gain_rate: Vec<f64>,
    /// Posterior mean multiplier of the loss rate of each character.
    pub site_loss_rate: Vec<f64>,
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
        MISSING => [1.0, 1.0],
        // soft evidence: the observation is "present" with probability q and "absent" otherwise
        k => {
            let q = presence_probability(k).unwrap_or(0.0);
            [q * f + (1.0 - q) * (1.0 - f), q * s + (1.0 - q) * (1.0 - s)]
        }
    }
}

/// The rate-1 transition matrix of the branch above `v`.
#[cfg(test)]
fn transition(p: &Params, v: NodeId) -> [[f64; 2]; 2] {
    transition_at(p.gain[v], p.loss[v], 1.0)
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
    never_seen: NeverSeen,
}

/// How the positions that were never seen (no intron at any observed genome) are modelled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NeverSeen {
    /// Each position is conditioned on being seen given its own pattern of missing cells: the probability that nobody shows it is
    /// computed over the genomes observable at that position. The default (SPEC §5.5).
    #[default]
    PerPattern,
    /// MALIN's correction: one probability `p0` that a position has no intron at any genome, all genomes being observable, for every
    /// position; the expected number of never-seen positions is `N p0 / (1 − p0)` and the log-likelihood is the uncorrected one minus
    /// `N ln(1 − p0)`. The same as `PerPattern` when no cell is missing.
    Global,
}

impl<'a> Ctx<'a> {
    #[cfg(test)]
    fn new(data: &'a ScoringData, err: &'a ErrorRates) -> Self {
        Self::with(data, err, NeverSeen::PerPattern)
    }

    fn with(data: &'a ScoringData, err: &'a ErrorRates, never_seen: NeverSeen) -> Self {
        let mut tip_of = vec![None; data.tree.len()];
        for (i, &t) in data.tips.iter().enumerate() {
            tip_of[t] = Some(i);
        }
        Ctx {
            tree: &data.tree,
            tip_of,
            err,
            never_seen,
        }
    }

    /// The "nobody shows an intron" column that goes with `col`.
    fn phantom_of(&self, col: &[u8]) -> Vec<u8> {
        match self.never_seen {
            NeverSeen::PerPattern => phantom(col),
            NeverSeen::Global => vec![ABSENT; col.len()],
        }
    }

    fn emit(&self, v: NodeId, col: &[u8]) -> [f64; 2] {
        match self.tip_of[v] {
            Some(i) => emission(col[i], self.err.sensitivity[i], self.err.false_rate[i]),
            None => [1.0, 1.0],
        }
    }

    /// Felsenstein pruning with per-node rescaling. Returns ln P(column) (unconditional).
    fn upward(&self, c: &Cat, col: &[u8], w: &mut Work) -> f64 {
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
                let t = c.tr[v];
                w.msg[v] = [
                    t[0][0] * l[0] + t[0][1] * l[1],
                    t[1][0] * l[0] + t[1][1] * l[1],
                ];
            }
        }
        let r = self.tree.root();
        let tot = (1.0 - c.root) * w.up[r][0] + c.root * w.up[r][1];
        log_scale + tot.ln()
    }

    /// Outside pass after [`Ctx::upward`] on the same column and class; adds `weight` × posterior expectations to `acc` (transition
    /// counts under class `k`).
    #[allow(clippy::too_many_arguments)]
    fn accumulate(
        &self,
        cat: &Cat,
        k: usize,
        col: &[u8],
        w: &mut Work,
        weight: f64,
        real: bool,
        acc: &mut Acc,
    ) {
        let tree = self.tree;
        let r = tree.root();
        w.down[r] = [1.0 - cat.root, cat.root];
        for &v in tree.postorder().iter().rev() {
            let dv = w.down[v];
            let post = norm([dv[0] * w.up[v][0], dv[1] * w.up[v][1]]);
            if real {
                acc.present[v] += weight * post[1];
                if let Some(i) = self.tip_of[v] {
                    if let Some(q) = presence_probability(col[i]) {
                        acc.missed[i] += weight * post[1] * (1.0 - q);
                        acc.false_[i] += weight * post[0] * q;
                    }
                }
            }
            if !real {
                // a never-seen position where this genome is observable and truly lacks the intron: where a false one could appear
                if let Some(i) = self.tip_of[v] {
                    if col[i] != MISSING {
                        acc.unseen_obs[i] += weight * post[0];
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
                let t = cat.tr[c];
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
                        acc.n[k][c][a][b] += weight * x;
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

impl<'a> Ctx<'a> {
    /// Felsenstein pruning that also counts independent origins (the root being present, or a 0 → 1 change on a branch), capped at
    /// 2. Returns the log of a common scale and `[P(data, 0 origins), P(data, 1 origin), P(data, 2 or more)]` divided by that scale.
    fn origin_pass(&self, c: &Cat, col: &[u8]) -> (f64, [f64; 3]) {
        type Cnt = [[f64; 3]; 2];
        let tree = self.tree;
        let mut log_scale = 0.0;
        let mut up: Vec<Cnt> = vec![[[0.0; 3]; 2]; tree.len()];
        let mut msg: Vec<Cnt> = vec![[[0.0; 3]; 2]; tree.len()];
        for &v in tree.postorder() {
            let e = self.emit(v, col);
            let mut l: Cnt = [[e[0], 0.0, 0.0], [e[1], 0.0, 0.0]];
            for &ch in tree.children(v) {
                let m = msg[ch];
                let mut out: Cnt = [[0.0; 3]; 2];
                for s in 0..2 {
                    for a in 0..3 {
                        for b in 0..3 {
                            out[s][(a + b).min(2)] += l[s][a] * m[s][b];
                        }
                    }
                }
                l = out;
            }
            let mx = l.iter().flatten().fold(0.0_f64, |a, &b| a.max(b));
            if mx <= 0.0 {
                return (f64::NEG_INFINITY, [0.0; 3]);
            }
            for row in l.iter_mut() {
                for x in row.iter_mut() {
                    *x /= mx;
                }
            }
            log_scale += mx.ln();
            up[v] = l;
            if tree.parent(v).is_some() {
                let t = c.tr[v];
                let mut m: Cnt = [[0.0; 3]; 2];
                for s in 0..2 {
                    for s2 in 0..2 {
                        let shift = usize::from(s == 0 && s2 == 1);
                        for k in 0..3 {
                            m[s][(k + shift).min(2)] += t[s][s2] * l[s2][k];
                        }
                    }
                }
                msg[v] = m;
            }
        }
        let l = up[tree.root()];
        let r = c.root;
        let p0 = (1.0 - r) * l[0][0];
        let p1 = (1.0 - r) * l[0][1] + r * l[1][0];
        let p2 = (1.0 - r) * l[0][2] + r * (l[1][1] + l[1][2]);
        (log_scale, [p0, p1, p2])
    }

    /// The class mixture of [`Ctx::origin_pass`], normalised to `P(data, n origins)` shares (they sum to 1 unless the column is
    /// impossible) together with `ln P(data)`.
    fn origin_mixture(&self, cs: &[Cat], col: &[u8]) -> (f64, [f64; 3]) {
        let parts: Vec<(f64, [f64; 3])> = cs.iter().map(|c| self.origin_pass(c, col)).collect();
        let m = parts.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
        if !m.is_finite() {
            return (f64::NEG_INFINITY, [0.0; 3]);
        }
        let mut tot = [0.0; 3];
        for (c, (ls, v)) in cs.iter().zip(&parts) {
            let f = c.w * (ls - m).exp();
            for j in 0..3 {
                tot[j] += f * v[j];
            }
        }
        let z: f64 = tot.iter().sum();
        (m + z.ln(), [tot[0] / z, tot[1] / z, tot[2] / z])
    }
}

/// How many independent origins a position has: the posterior given the data, and the evidence for convergence (SPEC §5.6).
#[derive(Clone, Copy, Debug)]
pub struct OriginPosterior {
    /// Posterior probability of no origin (the presences are annotation errors).
    pub p_none: f64,
    /// Posterior probability of exactly one origin: ancient (present at the root) or a single insertion.
    pub p_one: f64,
    /// Posterior probability of two or more independent origins (convergent insertion).
    pub p_multi: f64,
    /// `ln` of the likelihood ratio of "two or more origins" against "exactly one": the posterior odds divided by the prior odds
    /// under the fitted model. Positive: the pattern is better explained by convergence than a single origin is, per unit of prior.
    pub ln_lr_multi: f64,
}

/// Per-position evidence for convergent insertion against a single ancestral origin, under the fitted parameters `p`.
///
/// The model is the one fitted (rate classes included); an *origin* is a root presence or a 0 → 1 change on a branch. The prior odds
/// of "two or more" against "exactly one" are those of the model with every cell unobserved, so the likelihood ratio
/// `P(data | ≥ 2) / P(data | 1)` does not reward a position for being likely a priori.
pub fn origin_posteriors(data: &ScoringData, err: &ErrorRates, p: &Params) -> Vec<OriginPosterior> {
    let ctx = Ctx::new_for_origins(data, err);
    let cs = cats(p);
    let (_, prior) = ctx.origin_mixture(&cs, &vec![MISSING; data.tips.len()]);
    let prior_ln_odds = (prior[2] / prior[1]).ln();
    let pats = compress(data);
    let per: Vec<OriginPosterior> = pats
        .columns
        .iter()
        .map(|c| {
            let (_, q) = ctx.origin_mixture(&cs, c);
            OriginPosterior {
                p_none: q[0],
                p_one: q[1],
                p_multi: q[2],
                ln_lr_multi: (q[2] / q[1]).ln() - prior_ln_odds,
            }
        })
        .collect();
    pats.of_char.iter().map(|&k| per[k]).collect()
}

impl<'a> Ctx<'a> {
    fn new_for_origins(data: &'a ScoringData, err: &'a ErrorRates) -> Self {
        Self::with(data, err, NeverSeen::PerPattern)
    }
}

fn norm(x: [f64; 2]) -> [f64; 2] {
    let s = x[0] + x[1];
    [x[0] / s, x[1] / s]
}

/// The column with every observed cell turned to "absent": the event "nobody shows an intron" for this pattern of missing cells.
fn phantom(col: &[u8]) -> Vec<u8> {
    col.iter()
        .map(|&c| if c == MISSING { c } else { ABSENT })
        .collect()
}

#[derive(Clone)]
struct Acc {
    /// Expected transitions per class and branch, observed + unseen positions (for the M-step).
    n: Vec<Vec<[[f64; 2]; 2]>>,
    root: [f64; 2],
    /// Observed positions only.
    gains: Vec<f64>,
    losses: Vec<f64>,
    present: Vec<f64>,
    missed: Vec<f64>,
    false_: Vec<f64>,
    /// Never-seen positions per tip where the tip is observable.
    unseen_obs: Vec<f64>,
    /// Expected number of positions in each class, observed and unseen (for the level probabilities of free classes).
    class_w: Vec<f64>,
    log_lik: f64,
    unobserved: f64,
}

impl Acc {
    fn new(n: usize, tips: usize, classes: usize) -> Self {
        Acc {
            n: vec![vec![[[0.0; 2]; 2]; n]; classes],
            root: [0.0; 2],
            gains: vec![0.0; n],
            losses: vec![0.0; n],
            present: vec![0.0; n],
            missed: vec![0.0; tips],
            false_: vec![0.0; tips],
            unseen_obs: vec![0.0; tips],
            class_w: vec![0.0; classes],
            log_lik: 0.0,
            unobserved: 0.0,
        }
    }

    fn add(mut self, o: Acc) -> Acc {
        for (ka, kb) in self.n.iter_mut().zip(&o.n) {
            for (a, b) in ka.iter_mut().zip(kb) {
                for i in 0..2 {
                    for j in 0..2 {
                        a[i][j] += b[i][j];
                    }
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
        self.unseen_obs
            .iter_mut()
            .zip(&o.unseen_obs)
            .for_each(|(x, y)| *x += y);
        self.class_w
            .iter_mut()
            .zip(&o.class_w)
            .for_each(|(x, y)| *x += y);
        self.log_lik += o.log_lik;
        self.unobserved += o.unobserved;
        self
    }
}

/// ln of `Σ w_k exp(x_k)`.
fn ln_wsum_exp(x: &[f64], w: &[f64]) -> f64 {
    let m = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return m;
    }
    m + x
        .iter()
        .zip(w)
        .map(|(&a, &wk)| wk * (a - m).exp())
        .sum::<f64>()
        .ln()
}

/// One column over all classes: conditional ln-likelihood ln P(column) − ln(1 − Z), the no-show probability Z, and the posterior
/// mean rate multiplier.
fn column_terms(
    ctx: &Ctx,
    cs: &[Cat],
    rates: &[(f64, f64)],
    col: &[u8],
    w: &mut Work,
) -> (f64, f64, (f64, f64)) {
    let ph = ctx.phantom_of(col);
    let mut w2 = Work::new(ctx.tree.len());
    let lk: Vec<f64> = cs.iter().map(|c| ctx.upward(c, col, w)).collect();
    let z = cs
        .iter()
        .map(|c| c.w * ctx.upward(c, &ph, &mut w2).exp())
        .sum::<f64>()
        .min(1.0 - 1e-9);
    let ws: Vec<f64> = cs.iter().map(|c| c.w).collect();
    let ln_l = ln_wsum_exp(&lk, &ws);
    let (mut rg, mut rl) = (0.0, 0.0);
    for (k, (&a, &(g, l))) in lk.iter().zip(rates).enumerate() {
        let post = (a - ln_l).exp() * cs[k].w;
        rg += post * g;
        rl += post * l;
    }
    (ln_l - (1.0 - z).ln(), z, (rg, rl))
}

/// One E-step over all distinct columns. Columns are processed in fixed chunks whose partial sums are added in order, so the
/// result does not depend on how the work was split between threads (runs are reproducible bit for bit).
fn e_step(ctx: &Ctx, pats: &Patterns, p: &Params) -> Acc {
    let n = ctx.tree.len();
    let tips = ctx.err.sensitivity.len();
    let cs = cats(p);
    let kk = cs.len();
    let ws: Vec<f64> = cs.iter().map(|c| c.w).collect();
    let idx: Vec<usize> = (0..pats.columns.len()).collect();
    let parts: Vec<Acc> = idx
        .par_chunks(64)
        .map(|chunk| {
            let mut acc = Acc::new(n, tips, kk);
            let mut w: Vec<Work> = (0..kk).map(|_| Work::new(n)).collect();
            let mut wp: Vec<Work> = (0..kk).map(|_| Work::new(n)).collect();
            let (mut lk, mut zk) = (vec![0.0; kk], vec![0.0; kk]);
            for &k in chunk {
                let col = &pats.columns[k];
                let cnt = pats.counts[k];
                let ph = ctx.phantom_of(col);
                for c in 0..kk {
                    lk[c] = ctx.upward(&cs[c], col, &mut w[c]);
                    zk[c] = ctx.upward(&cs[c], &ph, &mut wp[c]).exp();
                }
                let ln_l = ln_wsum_exp(&lk, &ws);
                // Z is capped so that one nearly unobservable pattern cannot dominate the M-step (weight at most 1e9 per column)
                let zbar = zk.iter().zip(&ws).map(|(z, w)| z * w).sum::<f64>();
                let z = zbar.min(1.0 - 1e-9);
                let scale = if zbar > z { z / zbar } else { 1.0 };
                for c in 0..kk {
                    let post = (lk[c] - ln_l).exp() * ws[c];
                    acc.class_w[c] += cnt * post;
                    if post > 0.0 {
                        ctx.accumulate(&cs[c], c, col, &mut w[c], cnt * post, true, &mut acc);
                    }
                    let wt = cnt * scale * zk[c] * ws[c] / (1.0 - z);
                    acc.class_w[c] += wt;
                    if wt > 0.0 {
                        ctx.accumulate(&cs[c], c, &ph, &mut wp[c], wt, false, &mut acc);
                        acc.unobserved += wt;
                    }
                }
                acc.log_lik += cnt * (ln_l - (1.0 - z).ln());
            }
            acc
        })
        .collect();
    parts.into_iter().fold(Acc::new(n, tips, kk), Acc::add)
}

/// Expected complete-data log-likelihood of one branch with rate-1 probabilities `(g, l)`, over the classes.
fn q_branch(acc: &Acc, v: NodeId, rates: &[(f64, f64)], g: f64, l: f64) -> f64 {
    let mut q = 0.0;
    for (k, &(rg, rl)) in rates.iter().enumerate() {
        let t = transition_gl(g, l, rg, rl);
        let c = acc.n[k][v];
        for a in 0..2 {
            for b in 0..2 {
                if c[a][b] > 0.0 {
                    q += c[a][b] * t[a][b].max(1e-300).ln();
                }
            }
        }
    }
    q
}

/// Golden-section search for a maximum of `f` on `[lo, hi]`; returns the best point evaluated and its value.
fn golden_max(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> (f64, f64) {
    let gr = (5f64.sqrt() - 1.0) / 2.0;
    let mut x1 = hi - gr * (hi - lo);
    let mut x2 = lo + gr * (hi - lo);
    let (mut f1, mut f2) = (f(x1), f(x2));
    for _ in 0..100 {
        if hi - lo < 1e-10 {
            break;
        }
        if f1 < f2 {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + gr * (hi - lo);
            f2 = f(x2);
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - gr * (hi - lo);
            f1 = f(x1);
        }
    }
    if f1 >= f2 {
        (x1, f1)
    } else {
        (x2, f2)
    }
}

fn logit(x: f64) -> f64 {
    (x / (1.0 - x)).ln().clamp(-LOGIT_MAX, LOGIT_MAX)
}

fn expit(y: f64) -> f64 {
    1.0 / (1.0 + (-y).exp())
}

/// Search range of the branch coordinates (logits of `π` and `g + l`) under rate classes.
const LOGIT_MAX: f64 = 21.0;

/// The shape(s) in use, for the trace line.
fn shape_summary(r: &Rates) -> String {
    let lv = |v: &[(f64, f64)]| {
        v.iter()
            .map(|(p, m)| format!("{p:.3}@{m:.3}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    if r.is_free() {
        format!(
            "gain levels [{}]  loss levels [{}]",
            lv(&r.free_gain),
            lv(&r.free_loss)
        )
    } else if r.is_separate() {
        format!(
            "gain shape {:.4}  loss shape {:.4}",
            r.gain_shape, r.loss_shape
        )
    } else {
        format!("shape {:.4}", r.shape)
    }
}

/// Which gamma shapes are estimated.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct ShapeFit {
    shared: bool,
    gain: bool,
    loss: bool,
    /// Free levels: their probabilities and multipliers.
    free: bool,
}

/// Bounds of the multiplier of a free level.
const MULT_MIN: f64 = 1e-3;
const MULT_MAX: f64 = 1e3;

/// M-step for free levels: the level probabilities in closed form from the expected class counts (loss levels outer, gain inner), then
/// each non-reference multiplier by a one-dimensional search of the expected complete-data log-likelihood (kept only if it improves).
fn update_free_levels(acc: &Acc, p: &mut Params, free: &[NodeId]) {
    let (kg, kl) = (
        p.rates.free_gain.len().max(1),
        p.rates.free_loss.len().max(1),
    );
    let tot: f64 = acc.class_w.iter().sum();
    if tot > 0.0 {
        let normalise = |mut pr: Vec<f64>| {
            pr.iter_mut().for_each(|x| *x = x.max(1e-6));
            let t: f64 = pr.iter().sum();
            pr.iter_mut().for_each(|x| *x /= t);
            pr
        };
        if p.rates.free_loss.len() > 1 {
            let pr = normalise(
                (0..kl)
                    .map(|i| (0..kg).map(|j| acc.class_w[i * kg + j]).sum::<f64>() / tot)
                    .collect(),
            );
            for (lv, x) in p.rates.free_loss.iter_mut().zip(pr) {
                lv.0 = x;
            }
        }
        if p.rates.free_gain.len() > 1 {
            let pr = normalise(
                (0..kg)
                    .map(|j| (0..kl).map(|i| acc.class_w[i * kg + j]).sum::<f64>() / tot)
                    .collect(),
            );
            for (lv, x) in p.rates.free_gain.iter_mut().zip(pr) {
                lv.0 = x;
            }
        }
    }
    for dim in 0..2 {
        let k = if dim == 0 {
            p.rates.free_gain.len()
        } else {
            p.rates.free_loss.len()
        };
        for lvl in 1..k {
            let q = |ln_m: f64| {
                let mut r = p.rates.clone();
                let target = if dim == 0 {
                    &mut r.free_gain
                } else {
                    &mut r.free_loss
                };
                target[lvl].1 = ln_m.exp();
                let f = class_factors(&r);
                free.iter()
                    .map(|&v| q_branch(acc, v, &f, p.gain[v], p.loss[v]))
                    .sum::<f64>()
            };
            let cur_m = if dim == 0 {
                p.rates.free_gain[lvl].1
            } else {
                p.rates.free_loss[lvl].1
            };
            let cur = q(cur_m.ln());
            let (x, fx) = golden_max(q, MULT_MIN.ln(), MULT_MAX.ln());
            if fx > cur {
                let target = if dim == 0 {
                    &mut p.rates.free_gain
                } else {
                    &mut p.rates.free_loss
                };
                target[lvl].1 = x.exp();
            }
        }
    }
}

/// The new value of one shape (found by a one-dimensional search of the expected complete-data log-likelihood) when it improves on
/// the current one.
fn improve_shape(
    acc: &Acc,
    p: &Params,
    free: &[NodeId],
    get: fn(&Rates) -> f64,
    set: fn(&mut Rates, f64),
) -> Option<f64> {
    let q = |ln_a: f64| {
        let mut r = p.rates.clone();
        set(&mut r, ln_a.exp());
        let f = class_factors(&r);
        free.iter()
            .map(|&v| q_branch(acc, v, &f, p.gain[v], p.loss[v]))
            .sum::<f64>()
    };
    let cur = q(get(&p.rates).ln());
    let (x, fx) = golden_max(q, SHAPE_MIN.ln(), SHAPE_MAX.ln());
    (fx > cur).then(|| x.exp())
}

/// The M-step. Without classes it is closed form. With classes it is a generalised M-step: each estimated shape by a
/// one-dimensional search, then each branch in `(logit π, logit(g + l))` by safeguarded Newton steps (coordinate searches where
/// Newton cannot be used); a move is kept only if it raises the expected complete-data log-likelihood.
fn m_step(acc: &Acc, p: &mut Params, free: &[NodeId], fit: ShapeFit) {
    if p.rates.n_classes() <= 1 {
        for &v in free {
            let c = acc.n[0][v];
            let from0 = c[0][0] + c[0][1];
            let from1 = c[1][0] + c[1][1];
            if from0 > 0.0 {
                p.gain[v] = (c[0][1] / from0).clamp(CLAMP, 1.0 - CLAMP);
            }
            if from1 > 0.0 {
                p.loss[v] = (c[1][0] / from1).clamp(CLAMP, 1.0 - CLAMP);
            }
        }
    } else {
        if fit.shared {
            if let Some(a) = improve_shape(acc, p, free, |r| r.shape, |r, a| r.shape = a) {
                p.rates.shape = a;
            }
        }
        if fit.gain {
            if let Some(a) = improve_shape(acc, p, free, |r| r.gain_shape, |r, a| r.gain_shape = a)
            {
                p.rates.gain_shape = a;
            }
        }
        if fit.loss {
            if let Some(a) = improve_shape(acc, p, free, |r| r.loss_shape, |r, a| r.loss_shape = a)
            {
                p.rates.loss_shape = a;
            }
        }
        if fit.free {
            update_free_levels(acc, p, free);
        }
        let rates = class_factors(&p.rates);
        for &v in free {
            let to_gl = |x: f64, y: f64| {
                let (pi, s) = (expit(x), expit(y));
                (pi * s, (1.0 - pi) * s)
            };
            let q = |x: f64, y: f64| {
                let (g, l) = to_gl(x, y);
                q_branch(acc, v, &rates, g, l)
            };
            let (g0, l0) = (p.gain[v], p.loss[v]);
            let (mut x, mut y) = (logit(g0 / (g0 + l0)), logit(g0 + l0));
            let start = q(x, y);
            let mut best = start;
            // Newton steps on the 2 × 2 problem (finite differences), each kept only if it raises Q (step halving); where the
            // Hessian is not negative definite or Newton stalls, a round of one-dimensional searches instead
            for _ in 0..100 {
                let (x0, y0, f0) = (x, y, best);
                let h = 1e-4;
                let (fxp, fxm, fyp, fym) = (q(x + h, y), q(x - h, y), q(x, y + h), q(x, y - h));
                let gx = (fxp - fxm) / (2.0 * h);
                let gy = (fyp - fym) / (2.0 * h);
                let hxx = (fxp - 2.0 * best + fxm) / (h * h);
                let hyy = (fyp - 2.0 * best + fym) / (h * h);
                let hxy = (q(x + h, y + h) - q(x + h, y - h) - q(x - h, y + h) + q(x - h, y - h))
                    / (4.0 * h * h);
                let det = hxx * hyy - hxy * hxy;
                let mut moved = false;
                if hxx < 0.0 && det > 0.0 {
                    let (mut dx, mut dy) =
                        ((-hyy * gx + hxy * gy) / det, (hxy * gx - hxx * gy) / det);
                    for _ in 0..30 {
                        let (nx, ny) = (
                            (x + dx).clamp(-LOGIT_MAX, LOGIT_MAX),
                            (y + dy).clamp(-LOGIT_MAX, LOGIT_MAX),
                        );
                        let fnew = q(nx, ny);
                        if fnew > best {
                            x = nx;
                            y = ny;
                            best = fnew;
                            moved = true;
                            break;
                        }
                        dx *= 0.5;
                        dy *= 0.5;
                    }
                }
                if !moved {
                    let (nx, fx) = golden_max(|t| q(t, y), -LOGIT_MAX, LOGIT_MAX);
                    if fx > best {
                        x = nx;
                        best = fx;
                    }
                    let (ny, fy) = golden_max(|t| q(x, t), -LOGIT_MAX, LOGIT_MAX);
                    if fy > best {
                        y = ny;
                        best = fy;
                    }
                }
                if ((x - x0).abs() < 1e-10 && (y - y0).abs() < 1e-10)
                    || best - f0 < 1e-12 * best.abs().max(1.0)
                {
                    break;
                }
            }
            if best > start {
                let (g, l) = to_gl(x, y);
                p.gain[v] = g.clamp(CLAMP, 1.0 - CLAMP);
                p.loss[v] = l.clamp(CLAMP, 1.0 - CLAMP);
            }
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
    log_likelihood_with(data, err, p, NeverSeen::PerPattern)
}

/// The probability that a position shows no intron at any genome, all genomes being observable (MALIN's `p0`; the expected number
/// of never-seen positions is `N p0 / (1 − p0)`).
pub fn no_intron_probability(data: &ScoringData, err: &ErrorRates, p: &Params) -> f64 {
    let ctx = Ctx::with(data, err, NeverSeen::Global);
    let mut w = Work::new(data.tree.len());
    let ph = vec![ABSENT; data.tips.len()];
    cats(p)
        .iter()
        .map(|c| c.w * ctx.upward(c, &ph, &mut w).exp())
        .sum()
}

/// [`log_likelihood`] with the never-seen positions modelled as given.
pub fn log_likelihood_with(
    data: &ScoringData,
    err: &ErrorRates,
    p: &Params,
    never_seen: NeverSeen,
) -> f64 {
    let ctx = Ctx::with(data, err, never_seen);
    let mut w = Work::new(data.tree.len());
    let (cs, rates) = (cats(p), class_factors(&p.rates));
    data.characters
        .iter()
        .map(|ch| column_terms(&ctx, &cs, &rates, &ch.states, &mut w).0)
        .sum()
}

/// Nodes with at least one observed (non-`?`) cell among the tips below them, in any character.
///
/// A branch with nothing observed below it carries no information: its part of the likelihood sums to one whatever its
/// probabilities, so it is fixed to "no change" and reported as having no data instead of returning its starting values as estimates.
fn nodes_with_data(data: &ScoringData) -> Vec<bool> {
    let tree = &data.tree;
    let mut has = vec![false; tree.len()];
    for (i, &t) in data.tips.iter().enumerate() {
        has[t] = data.characters.iter().any(|c| c.states[i] != MISSING);
    }
    for &v in tree.postorder() {
        let below = tree.children(v).iter().any(|&c| has[c]);
        has[v] = has[v] || below;
    }
    has
}

/// Which branches are fixed at "no change", and why.
///
/// - no observed tip below the node: [`Fixed::NoData`];
/// - a non-root node with exactly one child that has data: its own branch, [`Fixed::Unary`] (only the product of the two branches is
///   identified; the change is put on the child's branch);
/// - at the root, after looking through such unary nodes: with two informative children, the larger one's branch (by tips),
///   [`Fixed::Root`]; with one, that one's branch, and the same rule again one level down (that child is then the effective root).
///   With three or more the parameters are identifiable as they are.
pub fn fixed_branches(data: &ScoringData) -> Vec<Fixed> {
    fixed_branches_with(data, true)
}

/// [`fixed_branches`] with the root rule optional. With one position class the root branches are not separately identifiable (only the
/// composite path between the root's two children matters), so the larger one is fixed and nothing is lost. With rate classes the two
/// branches are scaled differently in each class, so they are weakly identifiable and fixing one is a (small) restriction; leaving
/// them free (`root_rule` false, `FitOptions::free_root_branch`) adds a nearly flat direction and slows the fit, as MALIN's does.
pub fn fixed_branches_with(data: &ScoringData, root_rule: bool) -> Vec<Fixed> {
    let tree = &data.tree;
    let has = nodes_with_data(data);
    let mut fixed = vec![Fixed::No; tree.len()];
    let data_kids = |v: NodeId| -> Vec<NodeId> {
        tree.children(v)
            .iter()
            .copied()
            .filter(|&c| has[c])
            .collect()
    };
    for v in 0..tree.len() {
        if tree.parent(v).is_none() {
            continue;
        }
        if !has[v] {
            fixed[v] = Fixed::NoData;
        } else if !tree.children(v).is_empty() && data_kids(v).len() == 1 {
            fixed[v] = Fixed::Unary;
        }
    }
    // the root's effective children: follow each informative child down through unary nodes
    let unary: Vec<bool> = fixed.iter().map(|&f| f == Fixed::Unary).collect();
    let effective = |mut v: NodeId| {
        while unary[v] {
            v = data_kids(v)[0];
        }
        v
    };
    // walk down from the root while it has a single informative child: that child's state is then the root state, so its branch is
    // fixed and it becomes the effective root; with two informative children, fix the larger one's branch
    let mut r = tree.root();
    if root_rule {
        loop {
            let eff: Vec<NodeId> = data_kids(r).into_iter().map(effective).collect();
            match eff.len() {
                1 => {
                    fixed[eff[0]] = Fixed::Root;
                    r = eff[0];
                }
                2 => {
                    let n = |v: NodeId| tree.leaves_under(v).len();
                    fixed[if n(eff[1]) > n(eff[0]) {
                        eff[1]
                    } else {
                        eff[0]
                    }] = Fixed::Root;
                    break;
                }
                _ => break,
            }
        }
    }
    fixed
}

/// The parameters EM may change, as one vector: free gains, free losses, the root probability, and ln(shape) when it is estimated.
fn pack(p: &Params, free: &[NodeId], fit: ShapeFit) -> Vec<f64> {
    let mut x: Vec<f64> = free.iter().map(|&v| p.gain[v]).collect();
    x.extend(free.iter().map(|&v| p.loss[v]));
    x.push(p.root_present);
    if fit.shared {
        x.push(p.rates.shape.ln());
    }
    if fit.gain {
        x.push(p.rates.gain_shape.ln());
    }
    if fit.loss {
        x.push(p.rates.loss_shape.ln());
    }
    if fit.free {
        for lv in p.rates.free_gain.iter().chain(&p.rates.free_loss) {
            x.push(lv.0.ln());
        }
        for lv in p
            .rates
            .free_gain
            .iter()
            .skip(1)
            .chain(p.rates.free_loss.iter().skip(1))
        {
            x.push(lv.1.ln());
        }
    }
    x
}

/// Inverse of [`pack`], projecting into the parameter space (with classes, `g + l < 1` on every branch).
fn unpack(x: &[f64], free: &[NodeId], p: &mut Params, fit: ShapeFit) {
    let k = free.len();
    for (j, &v) in free.iter().enumerate() {
        p.gain[v] = x[j].clamp(CLAMP, 1.0 - CLAMP);
        p.loss[v] = x[k + j].clamp(CLAMP, 1.0 - CLAMP);
    }
    p.root_present = x[2 * k].clamp(CLAMP, 1.0 - CLAMP);
    let mut at = 2 * k + 1;
    let mut next = || {
        at += 1;
        x[at - 1].exp().clamp(SHAPE_MIN, SHAPE_MAX)
    };
    if fit.shared {
        p.rates.shape = next();
    }
    if fit.gain {
        p.rates.gain_shape = next();
    }
    if fit.loss {
        p.rates.loss_shape = next();
    }
    if fit.free {
        let (ng, nl) = (p.rates.free_gain.len(), p.rates.free_loss.len());
        // level probabilities were packed as logs; normalise within each dimension
        for (n, gain) in [(ng, true), (nl, false)] {
            if n == 0 {
                continue;
            }
            let raw: Vec<f64> = x[at..at + n].iter().map(|v| v.exp().max(1e-9)).collect();
            let t: f64 = raw.iter().sum();
            let lv = if gain {
                &mut p.rates.free_gain
            } else {
                &mut p.rates.free_loss
            };
            for (l, r) in lv.iter_mut().zip(raw) {
                l.0 = r / t;
            }
            at += n;
        }
        for lv in p
            .rates
            .free_gain
            .iter_mut()
            .skip(1)
            .chain(p.rates.free_loss.iter_mut().skip(1))
        {
            lv.1 = x[at].exp().clamp(MULT_MIN, MULT_MAX);
            at += 1;
        }
    }
    keep_below_one(p, free);
}

/// With rate classes, scale `(g, l)` down where `g + l` is not below 1 (the branch must be a continuous-time process).
fn keep_below_one(p: &mut Params, free: &[NodeId]) {
    if p.rates.n_classes() <= 1 {
        return;
    }
    let cap = 1.0 - 1e-6;
    for &v in free {
        let s = p.gain[v] + p.loss[v];
        if s > cap {
            p.gain[v] *= cap / s;
            p.loss[v] *= cap / s;
        }
    }
}

/// Fit the model by EM (accelerated with SQUAREM, Varadhan & Roland 2008, and never accepting a step that lowers the
/// log-likelihood) and return the parameters with the posterior expectations at the optimum.
///
/// Stops when one acceleration cycle improves the log-likelihood by less than `tol` *and* moves no parameter by more than `param_tol`
/// (default 1e-7);
/// the gain/root direction is flat with annotation error, so a log-likelihood criterion alone stops too early. `iterations` counts
/// E-steps.
/// The most E-steps one SQUAREM cycle can use: one EM step, up to six extrapolation candidates of two E-steps each, one final E-step.
const CYCLE_MAX_E_STEPS: usize = 1 + 6 * 2 + 1;

pub fn fit(data: &ScoringData, err: &ErrorRates, opt: &FitOptions) -> Result<FitResult, String> {
    let free = opt.free_gain_levels > 1 || opt.free_loss_levels > 1;
    if !free {
        return fit_once(data, err, opt, (Vec::new(), Vec::new()));
    }
    // Free levels make the likelihood multimodal (which level is slow, and how slow, trade off against the branch probabilities): start
    // from the levels carried by `init` when they match, otherwise from several spreads of the multipliers and keep the best fit.
    let spread = |k: usize, b: f64| -> Vec<(f64, f64)> {
        if k <= 1 {
            Vec::new()
        } else {
            (0..k).map(|i| (1.0 / k as f64, b.powi(i as i32))).collect()
        }
    };
    let matches = |given: &[(f64, f64)], k: usize| {
        if k <= 1 {
            given.is_empty()
        } else {
            given.len() == k
        }
    };
    let given = opt.init.as_ref().filter(|i| {
        matches(&i.rates.free_gain, opt.free_gain_levels)
            && matches(&i.rates.free_loss, opt.free_loss_levels)
    });
    let starts: Vec<LevelPair> = match given {
        Some(i) => vec![(i.rates.free_gain.clone(), i.rates.free_loss.clone())],
        None => [0.5, 0.2, 2.0, 5.0]
            .iter()
            .map(|&b| {
                (
                    spread(opt.free_gain_levels, b),
                    spread(opt.free_loss_levels, b),
                )
            })
            .collect(),
    };
    let mut best: Option<FitResult> = None;
    let mut first_err = None;
    for st in starts {
        match fit_once(data, err, opt, st) {
            Ok(r) => match &best {
                Some(b) if b.log_lik >= r.log_lik => {}
                _ => best = Some(r),
            },
            Err(e) => {
                first_err.get_or_insert(e);
            }
        }
    }
    best.ok_or_else(|| first_err.unwrap_or_default())
}

/// Free gain levels and free loss levels, each `(probability, multiplier)`.
type LevelPair = (Vec<(f64, f64)>, Vec<(f64, f64)>);

/// One EM fit from the given starting free levels (empty = none).
fn fit_once(
    data: &ScoringData,
    err: &ErrorRates,
    opt: &FitOptions,
    levels: LevelPair,
) -> Result<FitResult, String> {
    err.validate(data.tips.len())?;
    if opt.max_iter == 0 {
        return Err("max_iter must be at least 1 (the number of E-steps)".into());
    }
    if !(1..=MAX_CATEGORIES).contains(&opt.categories) {
        return Err(format!(
            "{} rate categories; use 1 (no rate variation) to {MAX_CATEGORIES}",
            opt.categories
        ));
    }
    for (what, k) in [("gain", opt.gain_categories), ("loss", opt.loss_categories)] {
        if !(1..=MAX_CATEGORIES).contains(&k) {
            return Err(format!(
                "{k} {what} rate categories; use 1 (no rate variation) to {MAX_CATEGORIES}"
            ));
        }
    }
    for (what, k) in [
        ("gain", opt.free_gain_levels),
        ("loss", opt.free_loss_levels),
    ] {
        if !(1..=MAX_CATEGORIES).contains(&k) {
            return Err(format!(
                "{k} free {what} levels; use 1 (no rate variation) to {MAX_CATEGORIES}"
            ));
        }
    }
    let free_levels = opt.free_gain_levels > 1 || opt.free_loss_levels > 1;
    if free_levels && (opt.categories > 1 || opt.gain_categories > 1 || opt.loss_categories > 1) {
        return Err(
            "free levels (--ml-gain-levels, --ml-loss-levels) cannot be combined with gamma classes (--ml-rate-categories, --ml-gain-categories, --ml-loss-categories)".into(),
        );
    }
    if free_levels && opt.free_gain_levels * opt.free_loss_levels > MAX_CATEGORIES {
        return Err(format!(
            "{} gain x {} loss levels = {} position classes; the most allowed is {MAX_CATEGORIES}",
            opt.free_gain_levels,
            opt.free_loss_levels,
            opt.free_gain_levels * opt.free_loss_levels
        ));
    }
    let separate = opt.gain_categories > 1 || opt.loss_categories > 1;
    if separate && opt.categories > 1 {
        return Err(
            "rate classes shared by gain and loss (categories) and separate gain/loss classes cannot be combined; use one or the other".into(),
        );
    }
    if separate && opt.gain_categories * opt.loss_categories > MAX_CATEGORIES {
        return Err(format!(
            "{} gain x {} loss classes = {} position classes; the most allowed is {MAX_CATEGORIES}",
            opt.gain_categories,
            opt.loss_categories,
            opt.gain_categories * opt.loss_categories
        ));
    }
    for a in [opt.shape, opt.gain_shape, opt.loss_shape]
        .into_iter()
        .flatten()
    {
        if !(SHAPE_MIN..=SHAPE_MAX).contains(&a) {
            return Err(format!(
                "gamma shape {a}; it must be in [{SHAPE_MIN}, {SHAPE_MAX}]"
            ));
        }
    }
    if data.characters.is_empty() {
        return Err("no characters to fit".into());
    }
    let tree = &data.tree;
    if let Some(v) = (0..tree.len()).find(|&v| tree.children(v).len() == 1) {
        return Err(format!(
            "node {v} has a single child; only the product of its two branches is identifiable (remove unary nodes from the tree)"
        ));
    }
    let ctx = Ctx::with(data, err, opt.never_seen);
    let pats = compress(data);
    let has_data = nodes_with_data(data);
    let informative = data.tips.iter().filter(|&&t| has_data[t]).count();
    if informative < 4 {
        return Err(format!(
            "{informative} tip(s) have data; the model needs at least 4 (conditioning on positions being seen costs one degree of \
             freedom, and with 3 tips the number of never-seen positions cannot be told from the root probability)"
        ));
    }
    let mut p = match &opt.init {
        Some(init) => {
            let ok = |x: f64| x > 0.0 && x < 1.0;
            if init.gain.len() != tree.len()
                || init.loss.len() != tree.len()
                || !ok(init.root_present)
                || (0..tree.len())
                    .any(|v| tree.parent(v).is_some() && !(ok(init.gain[v]) && ok(init.loss[v])))
            {
                return Err("starting parameters must be strictly between 0 and 1, one gain and loss per node".into());
            }
            init.clone()
        }
        None => Params::uniform(tree, 0.05, 0.2, 0.5),
    };
    p.rates = Rates {
        categories: opt.categories,
        shape: opt
            .shape
            .unwrap_or(p.rates.shape)
            .clamp(SHAPE_MIN, SHAPE_MAX),
        gain_categories: opt.gain_categories,
        gain_shape: opt
            .gain_shape
            .unwrap_or(p.rates.gain_shape)
            .clamp(SHAPE_MIN, SHAPE_MAX),
        loss_categories: opt.loss_categories,
        loss_shape: opt
            .loss_shape
            .unwrap_or(p.rates.loss_shape)
            .clamp(SHAPE_MIN, SHAPE_MAX),
        free_gain: levels.0,
        free_loss: levels.1,
    };
    let fit_shape = ShapeFit {
        shared: opt.categories > 1 && opt.shape.is_none(),
        gain: opt.gain_categories > 1 && opt.gain_shape.is_none(),
        loss: opt.loss_categories > 1 && opt.loss_shape.is_none(),
        free: free_levels,
    };
    let fixed = fixed_branches_with(data, !opt.free_root_branch);
    let free: Vec<NodeId> = (0..tree.len())
        .filter(|&v| tree.parent(v).is_some() && fixed[v] == Fixed::No)
        .collect();
    for v in 0..tree.len() {
        if tree.parent(v).is_some() && !free.contains(&v) {
            p.gain[v] = CLAMP;
            p.loss[v] = CLAMP;
        }
    }
    keep_below_one(&mut p, &free);
    // the EM map: parameters -> (parameters after one M-step, E-step at the input)
    let t_e = std::cell::Cell::new(0.0f64);
    let t_m = std::cell::Cell::new(0.0f64);
    let timed_e = |p: &Params| {
        let t = std::time::Instant::now();
        let a = e_step(&ctx, &pats, p);
        t_e.set(t_e.get() + t.elapsed().as_secs_f64());
        a
    };
    let timed_m = |acc: &Acc, q: &mut Params| {
        let t = std::time::Instant::now();
        m_step(acc, q, &free, fit_shape);
        t_m.set(t_m.get() + t.elapsed().as_secs_f64());
    };
    let em = |p: &Params| -> (Params, Acc) {
        let acc = timed_e(p);
        let mut q = p.clone();
        timed_m(&acc, &mut q);
        (q, acc)
    };
    let mut iterations = 0usize;
    let mut converged = false;
    let (mut p1, mut acc0) = em(&p);
    iterations += 1;
    while iterations < opt.max_iter {
        if opt.max_iter - iterations < CYCLE_MAX_E_STEPS {
            // Not enough budget left for a whole SQUAREM cycle: finish with plain EM steps, one E-step each, so that the limit is exact.
            let acc1 = timed_e(&p1);
            iterations += 1;
            let dlnl = acc1.log_lik - acc0.log_lik;
            let dx = pack(&p1, &free, fit_shape)
                .iter()
                .zip(pack(&p, &free, fit_shape).iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            p = p1.clone();
            acc0 = acc1;
            if opt.trace {
                eprintln!(
                    "likelihood: E-steps {iterations:>6}  lnL {:.6}  {}  dlnL {dlnl:.3e}  dx {dx:.3e}  (plain EM, iteration limit near)",
                    acc0.log_lik,
                    shape_summary(&p.rates)
                );
            }
            if dlnl < opt.tol && dx < opt.param_tol {
                converged = true;
                break;
            }
            let mut q = p.clone();
            timed_m(&acc0, &mut q);
            p1 = q;
            continue;
        }
        let (p2, acc1) = em(&p1);
        iterations += 1;
        let (x0, x1, x2) = (
            pack(&p, &free, fit_shape),
            pack(&p1, &free, fit_shape),
            pack(&p2, &free, fit_shape),
        );
        let r: Vec<f64> = x1.iter().zip(&x0).map(|(a, b)| a - b).collect();
        let v: Vec<f64> = x2
            .iter()
            .zip(&x1)
            .zip(&r)
            .map(|((a, b), c)| a - b - c)
            .collect();
        let (nr, nv) = (
            r.iter().map(|x| x * x).sum::<f64>().sqrt(),
            v.iter().map(|x| x * x).sum::<f64>().sqrt(),
        );
        // candidate: the SQUAREM extrapolation, falling back towards the plain double EM step (alpha = -1 gives x2)
        let mut best = (p2.clone(), None::<Acc>);
        if nv > 0.0 {
            let mut alpha = (-nr / nv).min(-1.0);
            for _ in 0..6 {
                let x: Vec<f64> = (0..x0.len())
                    .map(|i| x0[i] - 2.0 * alpha * r[i] + alpha * alpha * v[i])
                    .collect();
                let mut cand = p.clone();
                unpack(&x, &free, &mut cand, fit_shape);
                let (stab, _) = em(&cand); // one stabilising EM step
                iterations += 1;
                let acc_s = timed_e(&stab);
                iterations += 1;
                if acc_s.log_lik >= acc1.log_lik - 1e-9 {
                    best = (stab, Some(acc_s));
                    break;
                }
                if alpha >= -1.0 {
                    break;
                }
                alpha = (alpha - 1.0) / 2.0;
            }
        }
        let (new_p, new_acc) = match best {
            (q, Some(a)) => (q, a),
            (q, None) => {
                let a = timed_e(&q);
                iterations += 1;
                (q, a)
            }
        };
        let dlnl = new_acc.log_lik - acc0.log_lik;
        let dx = pack(&new_p, &free, fit_shape)
            .iter()
            .zip(&x0)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        p = new_p;
        acc0 = new_acc;
        if opt.trace {
            eprintln!(
                "likelihood: E-steps {iterations:>6}  lnL {:.6}  {}  dlnL {dlnl:.3e}  dx {dx:.3e}  E {:.1}s  M {:.1}s",
                acc0.log_lik,
                shape_summary(&p.rates),
                t_e.get(),
                t_m.get()
            );
        }
        if dlnl < opt.tol && dx < opt.param_tol {
            converged = true;
            break;
        }
        let mut q = p.clone();
        timed_m(&acc0, &mut q); // the next EM step comes from the E-step just done
        p1 = q;
    }
    let acc = acc0;
    let nodes = (0..tree.len())
        .map(|v| {
            let root = tree.parent(v).is_none();
            NodeStats {
                gain_prob: if root { p.root_present } else { p.gain[v] },
                loss_prob: if root { 0.0 } else { p.loss[v] },
                expected_gains: if root { acc.present[v] } else { acc.gains[v] },
                expected_losses: if root { 0.0 } else { acc.losses[v] },
                expected_present: acc.present[v],
                has_data: has_data[v],
                fixed: fixed[v],
            }
        })
        .collect();
    let tips = (0..data.tips.len())
        .map(|i| TipStats {
            expected_missed: acc.missed[i],
            expected_false: acc.false_[i],
            unseen_observable: acc.unseen_obs[i],
        })
        .collect();
    let mut w = Work::new(tree.len());
    let (cs, rates) = (cats(&p), class_factors(&p.rates));
    let per_pattern: Vec<(f64, (f64, f64))> = pats
        .columns
        .iter()
        .map(|c| {
            let (l, _, r) = column_terms(&ctx, &cs, &rates, c, &mut w);
            (l, r)
        })
        .collect();
    let site_log_lik = pats.of_char.iter().map(|&k| per_pattern[k].0).collect();
    let site_gain_rate: Vec<f64> = pats.of_char.iter().map(|&k| per_pattern[k].1 .0).collect();
    let site_loss_rate: Vec<f64> = pats.of_char.iter().map(|&k| per_pattern[k].1 .1).collect();
    let site_rate = site_gain_rate
        .iter()
        .zip(&site_loss_rate)
        .map(|(g, l)| 0.5 * (g + l))
        .collect();
    Ok(FitResult {
        params: p,
        log_lik: acc.log_lik,
        iterations,
        converged,
        nodes,
        tips,
        site_log_lik,
        site_rate,
        site_gain_rate,
        site_loss_rate,
        unobserved_positions: acc.unobserved,
    })
}

/// How a tip's false-intron rate is given.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FalseSpec {
    /// The per-position rate itself.
    Rate(f64),
    /// Annotation precision (fraction of annotated introns that are real); converted inside the fit.
    Precision(f64),
}

/// [`fit`] with false-intron rates given as precision where wanted.
///
/// Precision fixes the number of false introns in a genome, `carried·(1 − precision)`. The model applies the false-intron rate to
/// every position the genome truly lacks, including positions nobody shows (they exist but are not in the matrix), so the rate is
/// `carried·(1 − precision) / (truly absent positions in the matrix + never-seen positions where the genome is observable)`. The
/// never-seen count comes from the fit itself, so the rate and the fit are iterated to a fixed point (at most 20 rounds, until no
/// rate changes by more than 0.1 %), each fit starting from the last. Returns the fit, the rates used, and whether the rates
/// settled within the 20 rounds.
pub fn fit_with_precision(
    data: &ScoringData,
    sensitivity: &[f64],
    spec: &[FalseSpec],
    opt: &FitOptions,
) -> Result<(FitResult, ErrorRates, bool), String> {
    let n = data.tips.len();
    if sensitivity.len() != n || spec.len() != n {
        return Err(format!(
            "error specification for {} / {} tips, the tree has {n}",
            sensitivity.len(),
            spec.len()
        ));
    }
    let carried: Vec<f64> = (0..n)
        .map(|i| {
            data.characters
                .iter()
                .filter_map(|c| presence_probability(c.states[i]))
                .sum()
        })
        .collect();
    let observed: Vec<f64> = (0..n)
        .map(|i| {
            data.characters
                .iter()
                .filter(|c| c.states[i] != MISSING)
                .count() as f64
        })
        .collect();
    let rate = |i: usize, unseen: f64| -> f64 {
        match spec[i] {
            FalseSpec::Rate(f) => f,
            FalseSpec::Precision(pr) => {
                let false_ones = carried[i] * (1.0 - pr);
                let truly_absent =
                    (observed[i] - carried[i] * pr / sensitivity[i]).max(0.0) + unseen;
                // below the sensitivity, as `ErrorRates::validate` requires
                (false_ones / truly_absent.max(false_ones.max(1.0)))
                    .clamp(0.0, (0.5f64).min(0.95 * sensitivity[i]))
            }
        }
    };
    let mut err = ErrorRates {
        sensitivity: sensitivity.to_vec(),
        false_rate: (0..n).map(|i| rate(i, 0.0)).collect(),
    };
    let mut o = opt.clone();
    let mut r = fit(data, &err, &o)?;
    if spec.iter().all(|s| matches!(s, FalseSpec::Rate(_))) {
        return Ok((r, err, true));
    }
    let mut settled = false;
    let mut steps = r.iterations;
    for _ in 0..20 {
        let next: Vec<f64> = (0..n)
            .map(|i| rate(i, r.tips[i].unseen_observable))
            .collect();
        let change = next
            .iter()
            .zip(&err.false_rate)
            .map(|(a, b)| if *b > 0.0 { (a - b).abs() / b } else { a.abs() })
            .fold(0.0, f64::max);
        err.false_rate = next;
        o.init = Some(r.params.clone());
        r = fit(data, &err, &o)?;
        steps += r.iterations;
        if change < 1e-3 {
            settled = true;
            break;
        }
    }
    r.iterations = steps; // E-steps over all rounds, not only the last (warm-started) one
    Ok((r, err, settled))
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

    /// Sum over every assignment of true states: P(column, number of origins) with the origins counted directly (the root being
    /// present, plus every 0 -> 1 change on a branch), capped at 2.
    fn brute_origins(d: &ScoringData, err: &ErrorRates, p: &Params, col: &[u8]) -> [f64; 3] {
        let t = &d.tree;
        let n = t.len();
        let mut tip_of = vec![None; n];
        for (i, &v) in d.tips.iter().enumerate() {
            tip_of[v] = Some(i);
        }
        let mut out = [0.0; 3];
        for mask in 0u32..(1 << n) {
            let st = |v: usize| ((mask >> v) & 1) as usize;
            let mut prob = if st(t.root()) == 1 {
                p.root_present
            } else {
                1.0 - p.root_present
            };
            let mut origins = st(t.root());
            for v in 0..n {
                if let Some(u) = t.parent(v) {
                    prob *= transition(p, v)[st(u)][st(v)];
                    origins += usize::from(st(u) == 0 && st(v) == 1);
                }
            }
            let mut e = 1.0;
            for v in 0..n {
                if let Some(i) = tip_of[v] {
                    e *= emission(col[i], err.sensitivity[i], err.false_rate[i])[st(v)];
                }
            }
            out[origins.min(2)] += prob * e;
        }
        out
    }

    #[test]
    fn origin_posteriors_match_brute_force() {
        let cols: &[&[u8]] = &[
            &[1, 0, 1, MISSING, 0],
            &[0, 0, 0, 1, 1],
            &[1, 1, 1, 1, 1],
            &[1, 0, 0, 1, 0],
            &[MISSING, MISSING, 1, 0, MISSING],
            &[0, 0, 0, 0, 0],
        ];
        let d = data(NWK, cols);
        let err = some_errors();
        let p = some_params(&d.tree);
        let got = origin_posteriors(&d, &err, &p);
        let none = vec![MISSING; 5];
        let prior = brute_origins(&d, &err, &p, &none);
        let prior_ln = (prior[2] / prior[1]).ln();
        for (col, g) in cols.iter().zip(&got) {
            let b = brute_origins(&d, &err, &p, col);
            let z: f64 = b.iter().sum();
            let want = [b[0] / z, b[1] / z, b[2] / z];
            for (x, y) in [
                (g.p_none, want[0]),
                (g.p_one, want[1]),
                (g.p_multi, want[2]),
            ] {
                assert!((x - y).abs() < 1e-12, "{col:?}: {x} vs {y}");
            }
            let lr = (b[2] / b[1]).ln() - prior_ln;
            assert!((g.ln_lr_multi - lr).abs() < 1e-9, "{col:?}");
            assert!((g.p_none + g.p_one + g.p_multi - 1.0).abs() < 1e-12);
        }
        // sanity on the direction: a pattern shared by distant tips only (A and E) is more convergent than one shared by a clade
        let d2 = data(NWK, &[&[1, 1, 0, 0, 0], &[1, 0, 0, 0, 1]]);
        let g2 = origin_posteriors(&d2, &err, &p);
        assert!(g2[1].ln_lr_multi > g2[0].ln_lr_multi);
    }

    #[test]
    fn origin_posteriors_with_rate_classes_match_a_brute_force_mixture() {
        let t = tree(NWK);
        let mut p = some_params(&t);
        p.rates = Rates::shared(3, 0.6);
        let err = some_errors();
        let cols: &[&[u8]] = &[&[1, 0, 1, 0, 0], &[1, 1, 1, 1, 0], &[0, 1, 0, 0, 1]];
        let d = data(NWK, cols);
        let got = origin_posteriors(&d, &err, &p);
        let classes = class_list(&p.rates);
        let mix = |col: &[u8]| {
            let mut tot = [0.0; 3];
            for c in &classes {
                let b = brute_origins(&d, &err, &at_factors(&p, c.rg, c.rl), col);
                for j in 0..3 {
                    tot[j] += c.w * b[j];
                }
            }
            tot
        };
        let prior = mix(&[MISSING; 5]);
        for (col, g) in cols.iter().zip(&got) {
            let b = mix(col);
            let z: f64 = b.iter().sum();
            assert!((g.p_multi - b[2] / z).abs() < 1e-12, "{col:?}");
            assert!((g.p_one - b[1] / z).abs() < 1e-12, "{col:?}");
            let lr = (b[2] / b[1]).ln() - (prior[2] / prior[1]).ln();
            assert!((g.ln_lr_multi - lr).abs() < 1e-9, "{col:?}");
        }
    }

    #[test]
    fn soft_evidence_is_the_mixture_of_the_two_hard_cells() {
        use crate::score::{apply_soft_evidence, soft_state};
        let err = some_errors();
        let t = tree(NWK);
        let p = some_params(&t);
        let base: &[u8] = &[1, 0, 1, 0, 0];
        let mut d = data(NWK, &[base]);
        // tips 0 and 2 become soft; the tip order of `data` is the tree's tip order
        let names = d.tip_names.clone();
        let (q0, q2) = (0.25, 0.8);
        let n = apply_soft_evidence(
            &mut d,
            &[
                ("0:0".to_string(), names[0].clone(), q0),
                ("0:0".to_string(), names[2].clone(), q2),
            ],
        )
        .unwrap();
        assert_eq!(n, 2);
        let (r0, r2) = (
            presence_probability(soft_state(q0)).unwrap(),
            presence_probability(soft_state(q2)).unwrap(),
        );
        // quantised to 252 steps
        assert!((r0 - q0).abs() < 0.003 && (r2 - q2).abs() < 0.003);
        // P(soft cell) = q P(cell seen present) + (1 - q) P(cell seen absent), for each soft cell (multilinear)
        // compare the uncorrected likelihoods (the global correction's denominator is the same everywhere)
        let un = |dd: &ScoringData| {
            log_likelihood_with(dd, &err, &p, NeverSeen::Global)
                + (1.0 - no_intron_probability(dd, &err, &p)).ln()
        };
        let got_un = un(&d).exp();
        let want_un = r0 * r2 * un(&data(NWK, &[&[1, 0, 1, 0, 0]])).exp()
            + r0 * (1.0 - r2) * un(&data(NWK, &[&[1, 0, 0, 0, 0]])).exp()
            + (1.0 - r0) * r2 * un(&data(NWK, &[&[0, 0, 1, 0, 0]])).exp()
            + (1.0 - r0) * (1.0 - r2) * un(&data(NWK, &[&[0, 0, 0, 0, 0]])).exp();
        assert!(
            (got_un - want_un).abs() < 1e-12 * want_un.max(1.0),
            "{got_un} vs {want_un}"
        );
        // p = 0 and p = 1 are the hard cells; MISSING becomes an observation; unknown site or taxon is an error
        let mut d2 = data(NWK, &[&[MISSING, 0, 1, 0, 0]]);
        let nm = d2.tip_names.clone();
        apply_soft_evidence(&mut d2, &[("0:0".into(), nm[0].clone(), 1.0)]).unwrap();
        assert_eq!(d2.characters[0].states[0], PRESENT);
        apply_soft_evidence(&mut d2, &[("0:0".into(), nm[0].clone(), 0.0)]).unwrap();
        assert_eq!(d2.characters[0].states[0], ABSENT);
        assert!(apply_soft_evidence(&mut d2, &[("9:9".into(), nm[0].clone(), 0.5)]).is_err());
        assert!(apply_soft_evidence(&mut d2, &[("0:0".into(), "nobody".into(), 0.5)]).is_err());
        assert!(apply_soft_evidence(&mut d2, &[("0:0".into(), nm[0].clone(), 1.5)]).is_err());
    }

    #[test]
    fn soft_evidence_posteriors_and_the_fit_stay_consistent() {
        use crate::score::apply_soft_evidence;
        // posterior presence at a soft tip lies between its two hard cases; EM with soft cells never decreases the likelihood
        let err = some_errors();
        let cols: &[&[u8]] = &[
            &[1, 1, 1, 0, 0],
            &[0, 0, 1, 1, 1],
            &[1, 0, 0, 1, 0],
            &[1, 1, 0, 0, 0],
        ];
        let mut d = data(NWK, cols);
        let names = d.tip_names.clone();
        apply_soft_evidence(
            &mut d,
            &[
                ("0:0".into(), names[3].clone(), 0.5),
                ("2:0".into(), names[1].clone(), 0.3),
                ("3:0".into(), names[4].clone(), 0.9),
            ],
        )
        .unwrap();
        let r = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!(r.log_lik.is_finite() && r.converged);
        // soft evidence of 0.5 gives a likelihood between the two hard versions at the fitted parameters
        let ll = |c3: u8| {
            let mut cc: Vec<Vec<u8>> = cols.iter().map(|c| c.to_vec()).collect();
            cc[0][3] = c3;
            let refs: Vec<&[u8]> = cc.iter().map(|c| c.as_slice()).collect();
            let dd = data(NWK, &refs);
            log_likelihood(&dd, &err, &r.params)
        };
        let (a, b) = (ll(0), ll(1));
        let mid = {
            let mut dd = data(NWK, cols);
            apply_soft_evidence(&mut dd, &[("0:0".into(), names[3].clone(), 0.5)]).unwrap();
            log_likelihood(&dd, &err, &r.params)
        };
        assert!(
            mid > a.min(b) - 1e-9 && mid < a.max(b) + 1e-9,
            "{a} {mid} {b}"
        );
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
            let ln = ctx.upward(&cats(&p)[0], col, &mut w);
            assert!((ln - lik.ln()).abs() < 1e-10, "{ln} vs {}", lik.ln());
            let (cond, z2, _) = column_terms(&ctx, &cats(&p), &[(1.0, 1.0)], col, &mut w);
            assert!((z - z2).abs() < 1e-12 && (cond - (lik.ln() - (1.0 - z).ln())).abs() < 1e-10);
            let mut acc = Acc::new(d.tree.len(), 5, 1);
            ctx.upward(&cats(&p)[0], col, &mut w);
            ctx.accumulate(&cats(&p)[0], 0, col, &mut w, 1.0, true, &mut acc);
            for v in 0..d.tree.len() {
                assert!((acc.present[v] - present[v]).abs() < 1e-10, "node {v}");
                if d.tree.parent(v).is_some() {
                    for a in 0..2 {
                        for b in 0..2 {
                            assert!(
                                (acc.n[0][v][a][b] - pairs[v][a][b]).abs() < 1e-10,
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
            total += column_terms(&ctx, &cats(&p), &[(1.0, 1.0)], &col, &mut w)
                .0
                .exp();
        }
        assert!((total - 1.0).abs() < 1e-12, "{total}");
        // with a missing cell, the patterns over the other four cells sum to one as well
        let mut total = 0.0;
        for mask in 1u32..16 {
            let mut col: Vec<u8> = (0..4).map(|i| ((mask >> i) & 1) as u8).collect();
            col.insert(2, MISSING);
            total += column_terms(&ctx, &cats(&p), &[(1.0, 1.0)], &col, &mut w)
                .0
                .exp();
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
    fn max_iter_is_an_exact_limit_on_e_steps() {
        // The limit used to be checked once per SQUAREM cycle, so a fit could report up to 13 E-steps more than asked for.
        let t = tree(NWK);
        let truth = Params::uniform(&t, 0.05, 0.15, 0.7);
        let err = some_errors();
        let d = sim_data(NWK, simulate(&t, &truth, &err, 3_000, 11));
        let mut last = f64::NEG_INFINITY;
        for it in [1usize, 2, 3, 5, 13, 14, 15, 20, 27, 28, 29, 40, 41] {
            let r = fit(
                &d,
                &err,
                &FitOptions {
                    max_iter: it,
                    tol: 0.0, // never converges on the log-likelihood criterion, so the whole budget is used
                    param_tol: 0.0,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(r.iterations, it, "max_iter {it}");
            assert!(!r.converged);
            assert!(r.log_lik >= last - 1e-6, "{it}: {} < {last}", r.log_lik);
            last = r.log_lik;
        }
        // with the default stopping rule the count never exceeds the limit
        for it in [5usize, 30, 100] {
            let r = fit(
                &d,
                &err,
                &FitOptions {
                    max_iter: it,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(r.iterations <= it, "{} > {it}", r.iterations);
        }
        assert!(fit(
            &d,
            &err,
            &FitOptions {
                max_iter: 0,
                ..Default::default()
            }
        )
        .is_err());
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
                    ..Default::default()
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

    #[test]
    fn polytomy_posteriors_match_brute_force() {
        // an internal node with three children, a three-child root, missing cells (from the independent review)
        let nwk = "(A,(B,C,D),E);";
        let cols: &[&[u8]] = &[
            &[1, 0, 1, MISSING, 0],
            &[0, 1, 1, 1, 0],
            &[MISSING, 1, 0, 0, 1],
        ];
        let d = data(nwk, cols);
        let err = some_errors();
        let p = some_params(&d.tree);
        let ctx = Ctx::new(&d, &err);
        for col in cols {
            let (lik, _, present, pairs) = brute(&d, &err, &p, col);
            let mut w = Work::new(d.tree.len());
            assert!((ctx.upward(&cats(&p)[0], col, &mut w) - lik.ln()).abs() < 1e-10);
            let mut acc = Acc::new(d.tree.len(), 5, 1);
            ctx.accumulate(&cats(&p)[0], 0, col, &mut w, 1.0, true, &mut acc);
            for v in 0..d.tree.len() {
                assert!((acc.present[v] - present[v]).abs() < 1e-10, "present {v}");
                if d.tree.parent(v).is_some() {
                    for a in 0..2 {
                        for b in 0..2 {
                            assert!(
                                (acc.n[0][v][a][b] - pairs[v][a][b]).abs() < 1e-10,
                                "pair {v}"
                            );
                        }
                    }
                }
            }
        }
    }

    fn truth_with_fixed_root(t: &Tree, gain: f64, loss: f64, root: f64) -> Params {
        let mut p = Params::uniform(t, gain, loss, root);
        let f = fixed_root_branch(t).unwrap();
        p.gain[f] = CLAMP;
        p.loss[f] = CLAMP;
        p
    }

    #[test]
    fn a_tip_without_data_is_marked_and_adds_nothing() {
        let t = tree(NWK);
        let err = ErrorRates::none(5);
        let mut cols = simulate(
            &t,
            &truth_with_fixed_root(&t, 0.05, 0.15, 0.7),
            &err,
            5000,
            3,
        );
        for c in cols.iter_mut() {
            c[4] = MISSING; // E is never observed
        }
        cols.retain(|c| c.contains(&PRESENT));
        let d = sim_data(NWK, cols);
        let r = fit(&d, &err, &FitOptions::default()).unwrap();
        let e = t.tip_index("E").unwrap();
        assert!(!r.nodes[e].has_data && r.nodes.iter().filter(|n| !n.has_data).count() == 1);
        assert_eq!((r.params.gain[e], r.params.loss[e]), (CLAMP, CLAMP));
        assert!(
            r.nodes[e].expected_gains < 1e-3 && r.nodes[e].expected_losses < 1e-3,
            "{:?}",
            r.nodes[e]
        );
        // the same fit as with E left out of the data (its branch integrates out)
        let keep: Vec<NodeId> = d.tips[..4].to_vec();
        let d4 = ScoringData {
            tips: keep,
            tip_names: d.tip_names[..4].to_vec(),
            annotated_tip: vec![true; 4],
            characters: d
                .characters
                .iter()
                .map(|c| Character {
                    states: c.states[..4].to_vec(),
                    ..c.clone()
                })
                .collect(),
            tree: d.tree.clone(),
        };
        let r4 = fit(&d4, &ErrorRates::none(4), &FitOptions::default()).unwrap();
        assert!(
            (r.log_lik - r4.log_lik).abs() < 1e-4,
            "{} vs {}",
            r.log_lik,
            r4.log_lik
        );
    }

    #[test]
    fn precision_is_converted_with_the_never_seen_positions() {
        // tip A has 1 % false introns; its precision is measured from the simulated truth, then converted inside the fit
        let t = tree(NWK);
        let truth = truth_with_fixed_root(&t, 0.02, 0.15, 0.7);
        let err = ErrorRates {
            sensitivity: vec![0.9, 1.0, 0.9, 0.9, 0.9],
            false_rate: vec![0.01, 0.0, 0.0, 0.0, 0.0],
        };
        let tips = t.tips();
        let order: Vec<NodeId> = t.postorder().iter().rev().copied().collect();
        let mut rng = Rng(9);
        let (mut carried, mut real) = (0.0, 0.0);
        let mut cols = Vec::new();
        while cols.len() < 100_000 {
            let mut st = vec![0u8; t.len()];
            for &v in &order {
                st[v] = match t.parent(v) {
                    None => (rng.f() < truth.root_present) as u8,
                    Some(u) if st[u] == 0 => (rng.f() < truth.gain[v]) as u8,
                    Some(_) => (rng.f() >= truth.loss[v]) as u8,
                };
            }
            let col: Vec<u8> = tips
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    if st[v] == 1 {
                        (rng.f() < err.sensitivity[i]) as u8
                    } else {
                        (rng.f() < err.false_rate[i]) as u8
                    }
                })
                .collect();
            if col.contains(&PRESENT) {
                if col[0] == 1 {
                    carried += 1.0;
                    real += (st[tips[0]] == 1) as u8 as f64;
                }
                cols.push(col);
            }
        }
        let d = sim_data(NWK, cols);
        let mut spec = vec![FalseSpec::Rate(0.0); 5];
        spec[0] = FalseSpec::Precision(real / carried);
        let (r, used, settled) =
            fit_with_precision(&d, &err.sensitivity, &spec, &FitOptions::default()).unwrap();
        assert!(settled);
        let single = fit(&d, &used, &FitOptions::default()).unwrap();
        assert!(r.iterations > 0 && (r.log_lik - single.log_lik).abs() < 1e-3);
        let f = used.false_rate[0];
        assert!(
            (f - 0.01).abs() < 0.002,
            "converted false rate {f}, true 0.01 (the seen-only formula gave ~0.023)"
        );
        let false_true = carried - real;
        assert!(
            (r.tips[0].expected_false - false_true).abs() < 0.15 * false_true,
            "{} vs {false_true}",
            r.tips[0].expected_false
        );
        // rates given directly are passed through unchanged
        let (_, used, _) = fit_with_precision(
            &d,
            &err.sensitivity,
            &[FalseSpec::Rate(0.01); 5],
            &FitOptions::default(),
        )
        .unwrap();
        assert_eq!(used.false_rate, vec![0.01; 5]);
    }

    #[test]
    fn the_fit_is_a_stationary_point() {
        let t = tree(NWK);
        let err = some_errors();
        let d = sim_data(
            NWK,
            simulate(
                &t,
                &truth_with_fixed_root(&t, 0.05, 0.15, 0.7),
                &err,
                20_000,
                5,
            ),
        );
        let r = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!(r.converged, "{} iterations", r.iterations);
        let base = log_likelihood(&d, &err, &r.params);
        assert!((base - r.log_lik).abs() < 1e-6 * base.abs());
        let f = fixed_root_branch(&t).unwrap();
        for v in 0..t.len() {
            if t.parent(v).is_none() || v == f {
                continue;
            }
            for (which, h) in [(0, 1e-4), (0, -1e-4), (1, 1e-4), (1, -1e-4)] {
                let mut q = r.params.clone();
                if which == 0 {
                    q.gain[v] += h
                } else {
                    q.loss[v] += h
                }
                assert!(
                    log_likelihood(&d, &err, &q) <= base + 1e-6,
                    "moving node {v} ({which}) by {h} improves the fit"
                );
            }
        }
        for h in [1e-4, -1e-4] {
            let mut q = r.params.clone();
            q.root_present += h;
            assert!(
                log_likelihood(&d, &err, &q) <= base + 1e-6,
                "moving the root probability by {h} improves the fit"
            );
        }
        // a far stricter run (tolerance 1e-12) ends at the same log-likelihood
        let slow = {
            let o = FitOptions {
                max_iter: 400_000,
                tol: 1e-12,
                init: None,
                ..Default::default()
            };
            fit(&d, &err, &o).unwrap()
        };
        assert!(
            (slow.log_lik - r.log_lik).abs() < 1e-3,
            "{} vs {}",
            slow.log_lik,
            r.log_lik
        );
    }

    #[test]
    fn unary_nodes_and_bad_starting_values_are_refused() {
        let d = data("((A,B),((C)));", &[&[1, 0, 1]]);
        assert!(fit(&d, &ErrorRates::none(3), &FitOptions::default())
            .unwrap_err()
            .contains("single child"));
        let d = data(NWK, &[&[1, 0, 1, 0, 0]]);
        let mut p = Params::uniform(&d.tree, 0.05, 0.2, 0.5);
        p.loss[2] = 0.0;
        let o = FitOptions {
            init: Some(p),
            ..Default::default()
        };
        assert!(fit(&d, &ErrorRates::none(5), &o)
            .unwrap_err()
            .contains("strictly between"));
    }

    #[test]
    fn effectively_unary_nodes_are_fixed_and_the_split_does_not_depend_on_the_start() {
        // E never observed: (D,E) is left with one informative child, D. Found by the second review: without fixing one of the
        // two branches, their split moved with the starting values at the same log-likelihood.
        let t = tree(NWK);
        let err = some_errors();
        let mut cols = simulate(
            &t,
            &truth_with_fixed_root(&t, 0.05, 0.15, 0.7),
            &err,
            8000,
            13,
        );
        for c in cols.iter_mut() {
            c[4] = MISSING;
        }
        cols.retain(|c| c.contains(&PRESENT));
        let d = sim_data(NWK, cols);
        let fx = fixed_branches(&d);
        let (e, dn) = (t.tip_index("E").unwrap(), t.tip_index("D").unwrap());
        let de = t.parent(e).unwrap();
        assert_eq!(
            (fx[e], fx[de], fx[dn]),
            (Fixed::NoData, Fixed::Unary, Fixed::No)
        );
        assert_eq!(fx.iter().filter(|&&f| f == Fixed::Root).count(), 1);
        let fits: Vec<FitResult> = [(0.05, 0.2, 0.5), (0.2, 0.05, 0.8), (0.01, 0.4, 0.3)]
            .iter()
            .map(|&(g, l, r0)| {
                fit(
                    &d,
                    &err,
                    &FitOptions {
                        init: Some(Params::uniform(&t, g, l, r0)),
                        ..Default::default()
                    },
                )
                .unwrap()
            })
            .collect();
        for f in &fits[1..] {
            assert!((f.log_lik - fits[0].log_lik).abs() < 1e-4);
            for v in [dn, de] {
                let (a, b) = (f.nodes[v].expected_losses, fits[0].nodes[v].expected_losses);
                assert!(
                    (a - b).abs() < 1e-3 * b.max(1.0) + 0.5,
                    "node {v}: {a} vs {b}"
                );
            }
        }
        assert!(
            fits[0].nodes[de].expected_losses < 1e-2 && fits[0].nodes[de].fixed == Fixed::Unary
        );
    }

    #[test]
    fn a_root_with_one_informative_child_fixes_that_branch() {
        // C, D and E never observed: the root's only informative child is ((A,B),(F,G)), which then acts as the root
        let nwk = "(((A,B),(F,G)),(C,(D,E)));";
        let t = tree(nwk);
        let truth = truth_with_fixed_root(&t, 0.05, 0.15, 0.7);
        let mut cols = simulate(&t, &truth, &ErrorRates::none(7), 30_000, 17);
        let names: Vec<&str> = t.tips().iter().map(|&v| t.label(v).unwrap()).collect();
        for c in cols.iter_mut() {
            for (i, n) in names.iter().enumerate() {
                if ["C", "D", "E"].contains(n) {
                    c[i] = MISSING;
                }
            }
        }
        cols.retain(|c| c.contains(&PRESENT));
        let d = sim_data(nwk, cols);
        let fx = fixed_branches(&d);
        let top = t
            .parent(t.parent(t.tip_index("A").unwrap()).unwrap())
            .unwrap();
        assert_eq!(fx[top], Fixed::Root);
        assert_eq!(fx.iter().filter(|&&f| f == Fixed::Root).count(), 2); // and one of (A,B), (F,G)
        assert_eq!(fx.iter().filter(|&&f| f == Fixed::NoData).count(), 5); // C, D, E, (D,E), (C,(D,E))
        let a = fit(&d, &ErrorRates::none(7), &FitOptions::default()).unwrap();
        let b = fit(
            &d,
            &ErrorRates::none(7),
            &FitOptions {
                init: Some(Params::uniform(&t, 0.3, 0.02, 0.2)),
                ..Default::default()
            },
        )
        .unwrap();
        assert!((a.log_lik - b.log_lik).abs() < 1e-4);
        assert!(
            (a.params.root_present - b.params.root_present).abs() < 1e-3,
            "{} vs {}",
            a.params.root_present,
            b.params.root_present
        );
        // too few informative tips is refused
        let d3 = data("((A,B),(C,D));", &[&[1, 0, 1, MISSING]]);
        assert!(fit(&d3, &ErrorRates::none(4), &FitOptions::default())
            .unwrap_err()
            .contains("at least 4"));
    }

    #[test]
    fn e_steps_are_reproducible() {
        let t = tree(NWK);
        let d = sim_data(
            NWK,
            simulate(
                &t,
                &truth_with_fixed_root(&t, 0.05, 0.15, 0.7),
                &some_errors(),
                20_000,
                23,
            ),
        );
        let a = fit(&d, &some_errors(), &FitOptions::default()).unwrap();
        let b = fit(&d, &some_errors(), &FitOptions::default()).unwrap();
        assert_eq!(
            (a.iterations, a.log_lik.to_bits()),
            (b.iterations, b.log_lik.to_bits())
        );
    }

    // ---- rate variation across positions ----

    #[test]
    fn incomplete_gamma_and_discrete_gamma_match_known_values() {
        for x in [0.01, 0.3, 1.0, 2.5, 9.0] {
            assert!(
                (gamma_p(1.0, x) - (1.0 - (-x).exp())).abs() < 1e-13,
                "P(1, {x})"
            );
        }
        // P(2, x) = 1 − (1 + x) e^{−x}
        for x in [0.2, 1.7, 6.0] {
            assert!((gamma_p(2.0, x) - (1.0 - (1.0 + x) * (-x).exp())).abs() < 1e-13);
        }
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
        assert!((ln_gamma(10.0) - 362_880f64.ln()).abs() < 1e-11);
        for (a, q) in [(0.3, 0.25), (1.0, 0.5), (4.0, 0.75)] {
            assert!((gamma_p(a, gamma_quantile(a, q)) - q).abs() < 1e-10);
        }
        // Yang (1994), Table 1: mean rates of four categories
        let check = |shape: f64, want: [f64; 4]| {
            let r = category_rates(&Rates::shared(4, shape));
            for (a, b) in r.iter().zip(want) {
                assert!((a - b).abs() < 5e-4, "shape {shape}: {r:?}");
            }
            assert!((r.iter().sum::<f64>() / 4.0 - 1.0).abs() < 1e-12);
        };
        check(0.5, [0.0334, 0.2519, 0.8203, 2.8944]);
        check(1.0, [0.1369, 0.4768, 1.0000, 2.3863]);
        assert_eq!(category_rates(&Rates::shared(1, 0.3)), vec![1.0]);
        // extreme shapes stay finite and ordered
        for a in [SHAPE_MIN, SHAPE_MAX] {
            let r = category_rates(&Rates::shared(8, a));
            assert!(
                r.windows(2).all(|w| w[0] <= w[1]) && r.iter().all(|x| x.is_finite() && *x > 0.0),
                "{a}: {r:?}"
            );
        }
    }

    #[test]
    fn class_transitions_form_a_continuous_time_process() {
        let (g, l) = (0.12, 0.3);
        assert_eq!(transition_at(g, l, 1.0), [[1.0 - g, g], [l, 1.0 - l]]);
        let mul = |a: [[f64; 2]; 2], b: [[f64; 2]; 2]| {
            let mut c = [[0.0; 2]; 2];
            for i in 0..2 {
                for j in 0..2 {
                    c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j];
                }
            }
            c
        };
        // T(r1) T(r2) = T(r1 + r2), and rate 1 reproduces (g, l)
        for (r1, r2) in [(0.3, 0.7), (0.05, 2.2), (1.5, 1.5)] {
            let a = mul(transition_at(g, l, r1), transition_at(g, l, r2));
            let b = transition_at(g, l, r1 + r2);
            for i in 0..2 {
                for j in 0..2 {
                    assert!((a[i][j] - b[i][j]).abs() < 1e-12);
                }
            }
        }
        let t = transition_at(g, l, 1.0 + 1e-12);
        assert!((t[0][1] - g).abs() < 1e-10 && (t[1][0] - l).abs() < 1e-10);
    }

    /// `p` with four rate classes of shape `shape`.
    fn with_classes(mut p: Params, k: usize, shape: f64) -> Params {
        p.rates = Rates::shared(k, shape);
        p
    }

    /// `p` with separate gain and loss classes.
    fn with_separate(mut p: Params, kg: usize, ag: f64, kl: usize, al: f64) -> Params {
        p.rates = Rates::separate(kg, ag, kl, al);
        p
    }

    /// `p` with the gain rate multiplied by `rg` and the loss rate by `rl` on every branch, without classes (for brute force).
    fn at_factors(p: &Params, rg: f64, rl: f64) -> Params {
        let mut q = Params::uniform(&tree(NWK), 0.1, 0.1, p.root_present);
        q.gain = vec![0.0; p.gain.len()];
        q.loss = vec![0.0; p.gain.len()];
        for v in 0..p.gain.len() {
            let t = transition_gl(p.gain[v], p.loss[v], rg, rl);
            q.gain[v] = t[0][1];
            q.loss[v] = t[1][0];
        }
        q
    }

    /// The class mixture (likelihood, no-show probability, posterior mean factors, posterior presence) of `p` against brute-force
    /// enumeration under each class, and conditional probabilities of all seen patterns summing to one.
    fn check_class_mixture(p: &Params) {
        let cols: &[&[u8]] = &[&[1, 0, 1, MISSING, 0], &[0, 0, 0, 1, 1], &[1, 1, 1, 1, 1]];
        let d = data(NWK, cols);
        let err = some_errors();
        let classes = class_list(&p.rates);
        let factors: Vec<(f64, f64)> = classes.iter().map(|c| (c.rg, c.rl)).collect();
        let cs = cats(p);
        let ctx = Ctx::new(&d, &err);
        let mut w = Work::new(d.tree.len());
        for col in cols {
            let (mut lik, mut z, mut post_g, mut post_l) = (0.0, 0.0, 0.0, 0.0);
            for c in &classes {
                let (l, zz, _, _) = brute(&d, &err, &at_factors(p, c.rg, c.rl), col);
                lik += c.w * l;
                z += c.w * zz;
                post_g += c.rg * c.w * l;
                post_l += c.rl * c.w * l;
            }
            let (cond, z2, (rg, rl)) = column_terms(&ctx, &cs, &factors, col, &mut w);
            assert!((z - z2).abs() < 1e-12);
            assert!((cond - (lik.ln() - (1.0 - z).ln())).abs() < 1e-10);
            assert!((rg - post_g / lik).abs() < 1e-10, "gain factor");
            assert!((rl - post_l / lik).abs() < 1e-10, "loss factor");
        }
        // the E-step's posterior presence is the class mixture of the brute-force posteriors
        let d1 = data(NWK, &[cols[0]]);
        let acc = e_step(&ctx, &compress(&d1), p);
        let mut present = vec![0.0; d.tree.len()];
        let mut tot = 0.0;
        for c in &classes {
            let (l, _, pr, _) = brute(&d, &err, &at_factors(p, c.rg, c.rl), cols[0]);
            tot += c.w * l;
            for v in 0..d.tree.len() {
                present[v] += c.w * l * pr[v];
            }
        }
        for v in 0..d.tree.len() {
            assert!(
                (acc.present[v] - present[v] / tot).abs() < 1e-10,
                "node {v}"
            );
        }
        // conditional probabilities of all seen patterns still sum to one
        let mut total = 0.0;
        for mask in 1u32..32 {
            let col: Vec<u8> = (0..5).map(|i| ((mask >> i) & 1) as u8).collect();
            total += column_terms(&ctx, &cs, &factors, &col, &mut w).0.exp();
        }
        assert!((total - 1.0).abs() < 1e-12, "{total}");
    }

    #[test]
    fn class_mixture_matches_brute_force_and_sums_to_one() {
        let t = tree(NWK);
        check_class_mixture(&with_classes(some_params(&t), 4, 0.6));
    }

    #[test]
    fn separate_gain_and_loss_classes_match_brute_force_and_sum_to_one() {
        let t = tree(NWK);
        // loss classes only (the structure of MALIN's published rate-variation runs), gain classes only, and the grid of both
        check_class_mixture(&with_separate(some_params(&t), 1, 1.0, 3, 0.5));
        check_class_mixture(&with_separate(some_params(&t), 3, 0.7, 1, 1.0));
        check_class_mixture(&with_separate(some_params(&t), 2, 0.4, 3, 1.5));
    }

    #[test]
    fn global_never_seen_correction_matches_brute_force() {
        // MALIN's correction: log-likelihood = sum of ln P(column) minus N ln(1 - p0), p0 = P(no intron at any genome, all observable)
        let t = tree(NWK);
        let p = some_params(&t);
        let err = some_errors();
        let with_missing: &[&[u8]] = &[
            &[1, 0, 1, MISSING, 0],
            &[0, 0, 0, 1, 1],
            &[1, MISSING, 1, 1, MISSING],
            &[1, 1, 1, 1, 1],
        ];
        let d = data(NWK, with_missing);
        let (_, p0, _, _) = brute(&d, &err, &p, &[ABSENT; 5]);
        let want: f64 = with_missing
            .iter()
            .map(|c| brute(&d, &err, &p, c).0.ln() - (1.0 - p0).ln())
            .sum();
        let got = log_likelihood_with(&d, &err, &p, NeverSeen::Global);
        assert!((got - want).abs() < 1e-10, "{got} vs {want}");
        // with missing cells the per-pattern conditioning differs; without them the two are the same function
        let per = log_likelihood_with(&d, &err, &p, NeverSeen::PerPattern);
        assert!((per - got).abs() > 1e-4, "{per} vs {got}");
        let full: &[&[u8]] = &[&[1, 0, 1, 0, 0], &[0, 0, 0, 1, 1], &[1, 1, 1, 1, 1]];
        let df = data(NWK, full);
        let (a, b) = (
            log_likelihood_with(&df, &err, &p, NeverSeen::Global),
            log_likelihood_with(&df, &err, &p, NeverSeen::PerPattern),
        );
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        // the expected number of never-seen positions is N p0 / (1 - p0), as MALIN's "predicted missing"
        let ctx = Ctx::with(&d, &err, NeverSeen::Global);
        let acc = e_step(&ctx, &compress(&d), &p);
        let n = with_missing.len() as f64;
        assert!(
            (acc.unobserved - n * p0 / (1.0 - p0)).abs() < 1e-9 * n,
            "{} vs {}",
            acc.unobserved,
            n * p0 / (1.0 - p0)
        );
        assert!((no_intron_probability(&d, &err, &p) - p0).abs() < 1e-12);
    }

    #[test]
    fn em_with_the_global_correction_never_decreases_and_reaches_a_stationary_point() {
        // data with missing cells: simulate, then hide 15 % of the cells
        let t = tree(NWK8);
        let err = ErrorRates {
            sensitivity: vec![0.95; 8],
            false_rate: vec![0.005; 8],
        };
        let truth = truth8(0.5);
        let mut rng = Rng(31);
        let cols: Vec<Vec<u8>> = simulate(&t, &truth, &err, 20_000, 31)
            .into_iter()
            .map(|c| {
                c.into_iter()
                    .map(|x| if rng.f() < 0.15 { MISSING } else { x })
                    .collect::<Vec<u8>>()
            })
            .filter(|c| c.contains(&PRESENT))
            .collect();
        let d = sim_data(NWK8, cols);
        let opt = |it: usize| FitOptions {
            max_iter: it,
            tol: 0.0,
            never_seen: NeverSeen::Global,
            ..Default::default()
        };
        let mut last = f64::NEG_INFINITY;
        for it in [1, 2, 4, 8, 16, 32] {
            let r = fit(&d, &err, &opt(it)).unwrap();
            assert!(
                r.log_lik >= last - 1e-6,
                "iteration {it}: {} < {last}",
                r.log_lik
            );
            last = r.log_lik;
        }
        let r = fit(
            &d,
            &err,
            &FitOptions {
                never_seen: NeverSeen::Global,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(r.converged);
        let base = log_likelihood_with(&d, &err, &r.params, NeverSeen::Global);
        assert!((base - r.log_lik).abs() < 1e-6 * base.abs());
        for v in 0..t.len() {
            if r.nodes[v].fixed != Fixed::No || t.parent(v).is_none() {
                continue;
            }
            for (which, h) in [(0, 1e-4), (0, -1e-4), (1, 1e-4), (1, -1e-4)] {
                let mut q = r.params.clone();
                if which == 0 {
                    q.gain[v] += h
                } else {
                    q.loss[v] += h
                }
                assert!(
                    log_likelihood_with(&d, &err, &q, NeverSeen::Global) <= base + 1e-4,
                    "node {v} ({which}) by {h}"
                );
            }
        }
        // the two corrections differ on data with missing cells, so their fits do too
        let per = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!((per.log_lik - r.log_lik).abs() > 1e-3);
    }

    #[test]
    fn free_levels_match_brute_force_and_their_probabilities_multiply() {
        let t = tree(NWK);
        let mut p = some_params(&t);
        // loss levels like MALIN's published runs (probability 0.2 at multiplier 0.4, 0.8 at the reference), and a grid with gain levels
        p.rates = Rates::free(vec![], vec![(0.8, 1.0), (0.2, 0.4)]);
        check_class_mixture(&p);
        p.rates = Rates::free(
            vec![(0.6, 1.0), (0.4, 2.5)],
            vec![(0.5, 1.0), (0.3, 0.5), (0.2, 0.1)],
        );
        check_class_mixture(&p);
        let c = class_list(&p.rates);
        assert_eq!(c.len(), 6);
        assert!((c.iter().map(|x| x.w).sum::<f64>() - 1.0).abs() < 1e-12);
        // loss level outer, gain level inner; probability = product
        assert_eq!((c[1].rg, c[1].rl), (2.5, 1.0));
        assert!((c[1].w - 0.4 * 0.5).abs() < 1e-12);
        assert_eq!((c[5].rg, c[5].rl), (2.5, 0.1));
        assert!((c[5].w - 0.4 * 0.2).abs() < 1e-12);
        assert_eq!(p.rates.n_classes(), 6);
        // one free dimension leaves the other at the single reference level
        let q = Rates::free(vec![], vec![(0.8, 1.0), (0.2, 0.4)]);
        assert!(class_list(&q).iter().all(|x| x.rg == 1.0));
        assert!(class_list(&Rates::default()).len() == 1);
    }

    #[test]
    fn em_with_free_loss_levels_recovers_probability_and_multiplier() {
        // MALIN's published structure: a reference loss level (probability 0.8) and a slow one (0.2 at multiplier 0.4)
        let t = tree(NWK8);
        let err = ErrorRates {
            sensitivity: vec![0.95; 8],
            false_rate: vec![0.005; 8],
        };
        let mut truth = truth8(0.5);
        truth.rates = Rates::free(vec![], vec![(0.8, 1.0), (0.2, 0.4)]);
        let d = sim_data(NWK8, simulate_classes(&t, &truth, &err, 40_000, 29));
        let opt = |it: usize| FitOptions {
            max_iter: it,
            tol: 0.0,
            free_loss_levels: 2,
            ..Default::default()
        };
        let mut last = f64::NEG_INFINITY;
        for it in [1, 2, 4, 8, 16, 32] {
            let r = fit(&d, &err, &opt(it)).unwrap();
            assert!(
                r.log_lik >= last - 1e-6,
                "iteration {it}: {} < {last}",
                r.log_lik
            );
            last = r.log_lik;
        }
        let r = fit(
            &d,
            &err,
            &FitOptions {
                free_loss_levels: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(r.converged, "{} iterations", r.iterations);
        let lv = &r.params.rates.free_loss;
        assert_eq!(lv.len(), 2);
        assert!((lv[0].0 + lv[1].0 - 1.0).abs() < 1e-9 && lv[0].1 == 1.0);
        assert!(
            (0.05..0.4).contains(&lv[1].0),
            "slow level probability {}",
            lv[1].0
        );
        assert!(
            (0.15..0.8).contains(&lv[1].1),
            "slow level multiplier {}",
            lv[1].1
        );
        assert!(r.log_lik >= log_likelihood(&d, &err, &truth) - 1e-6);
        let flat = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!(
            r.log_lik > flat.log_lik + 50.0,
            "{} vs {}",
            r.log_lik,
            flat.log_lik
        );
        // stationary in the level probability, the multiplier and the branches
        let base = log_likelihood(&d, &err, &r.params);
        assert!((base - r.log_lik).abs() < 1e-6 * base.abs());
        for h in [0.01, -0.01] {
            let mut q = r.params.clone();
            q.rates.free_loss[1].1 *= 1.0 + h;
            assert!(
                log_likelihood(&d, &err, &q) <= base + 1e-4,
                "multiplier × {h}"
            );
            let mut q = r.params.clone();
            q.rates.free_loss[1].0 += h * 0.1;
            q.rates.free_loss[0].0 -= h * 0.1;
            assert!(
                log_likelihood(&d, &err, &q) <= base + 1e-4,
                "probability + {h}"
            );
        }
        // option checks
        let bad = |o: FitOptions| fit(&d, &err, &o).unwrap_err();
        assert!(bad(FitOptions {
            free_loss_levels: 2,
            categories: 2,
            ..Default::default()
        })
        .contains("cannot be combined"));
        assert!(bad(FitOptions {
            free_loss_levels: 2,
            loss_categories: 2,
            ..Default::default()
        })
        .contains("cannot be combined"));
        assert!(bad(FitOptions {
            free_gain_levels: 5,
            free_loss_levels: 5,
            ..Default::default()
        })
        .contains("position classes"));
        assert!(bad(FitOptions {
            free_gain_levels: 0,
            ..Default::default()
        })
        .contains("free gain levels"));
    }

    #[test]
    fn separate_factors_define_a_two_state_process_with_scaled_rates() {
        // The branch (g, l) at rate 1 is a two-state process with gain rate a and loss rate b over unit time. With factors (rg, rl) its
        // transition matrix must be the matrix exponential of the generator with rates (a·rg, b·rl): the closed form of a two-state chain.
        let two_state = |a: f64, b: f64| -> [[f64; 2]; 2] {
            let (s, e) = (a + b, (-(a + b)).exp());
            [
                [(b + a * e) / s, a * (1.0 - e) / s],
                [b * (1.0 - e) / s, (a + b * e) / s],
            ]
        };
        for (g, l) in [(0.12, 0.3), (0.02, 0.5), (0.4, 0.05)] {
            let (s, pi) = (g + l, g / (g + l));
            let tau = -(1.0f64 - s).ln();
            let (a, b) = (pi * tau, (1.0 - pi) * tau);
            // rate 1 reproduces (g, l) themselves
            let t1 = two_state(a, b);
            assert!((t1[0][1] - g).abs() < 1e-12 && (t1[1][0] - l).abs() < 1e-12);
            for (rg, rl) in [
                (1.0, 1.0),
                (0.3, 1.0),
                (1.0, 0.3),
                (2.5, 0.4),
                (0.05, 3.0),
                (1.7, 1.7),
            ] {
                let want = two_state(a * rg, b * rl);
                let got = transition_gl(g, l, rg, rl);
                for i in 0..2 {
                    for j in 0..2 {
                        assert!(
                            (want[i][j] - got[i][j]).abs() < 1e-12,
                            "g {g} l {l} rg {rg} rl {rl}: {want:?} vs {got:?}"
                        );
                    }
                }
            }
            // equal factors are the shared-factor form
            assert_eq!(transition_gl(g, l, 0.8, 0.8), transition_at(g, l, 0.8));
        }
    }

    #[test]
    fn class_factors_form_the_expected_grid_and_have_mean_one() {
        assert_eq!(class_factors(&Rates::default()), vec![(1.0, 1.0)]);
        let shared = class_factors(&Rates::shared(4, 0.6));
        assert!(shared.iter().all(|&(g, l)| g == l) && shared.len() == 4);
        let sep = Rates::separate(2, 0.4, 3, 1.5);
        let f = class_factors(&sep);
        assert_eq!(f.len(), 6);
        assert_eq!(sep.n_classes(), 6);
        // loss classes are the outer index, gain classes the inner one
        let (g2, l3) = (
            category_rates(&Rates::shared(2, 0.4)),
            category_rates(&Rates::shared(3, 1.5)),
        );
        for (i, &(rg, rl)) in f.iter().enumerate() {
            assert_eq!((rg, rl), (g2[i % 2], l3[i / 2]));
        }
        assert!((f.iter().map(|x| x.0).sum::<f64>() / 6.0 - 1.0).abs() < 1e-12);
        assert!((f.iter().map(|x| x.1).sum::<f64>() / 6.0 - 1.0).abs() < 1e-12);
        // one separate dimension with a single class leaves the other factor at 1
        assert!(class_factors(&Rates::separate(1, 1.0, 3, 1.5))
            .iter()
            .all(|x| x.0 == 1.0));
    }

    const NWK8: &str = "(((A,B),(C,D)),((E,F),(G,H)));";

    /// Simulate seen columns with a gamma rate class drawn per position (no missing cells).
    fn simulate_classes(
        t: &Tree,
        p: &Params,
        err: &ErrorRates,
        n: usize,
        seed: u64,
    ) -> Vec<Vec<u8>> {
        let cs = cats(p);
        let tips = t.tips();
        let mut rng = Rng(seed);
        let order: Vec<NodeId> = t.postorder().iter().rev().copied().collect();
        let mut out = Vec::new();
        while out.len() < n {
            let u = rng.f();
            let mut acc_w = 0.0;
            let mut pick = cs.len() - 1;
            for (k, c) in cs.iter().enumerate() {
                acc_w += c.w;
                if u < acc_w {
                    pick = k;
                    break;
                }
            }
            let c = &cs[pick];
            let mut st = vec![0u8; t.len()];
            for &v in &order {
                st[v] = match t.parent(v) {
                    None => (rng.f() < c.root) as u8,
                    Some(u) => {
                        let tr = c.tr[v][st[u] as usize];
                        (rng.f() < tr[1]) as u8
                    }
                };
            }
            let col: Vec<u8> = tips
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    if st[v] == 1 {
                        (rng.f() < err.sensitivity[i]) as u8
                    } else {
                        (rng.f() < err.false_rate[i]) as u8
                    }
                })
                .collect();
            if col.contains(&PRESENT) {
                out.push(col);
            }
        }
        out
    }

    fn truth8(shape: f64) -> Params {
        let t = tree(NWK8);
        let mut p = with_classes(Params::uniform(&t, 0.08, 0.25, 0.6), 4, shape);
        for v in 0..t.len() {
            p.gain[v] = 0.04 + 0.01 * (v % 5) as f64;
            p.loss[v] = 0.12 + 0.05 * (v % 4) as f64;
        }
        let f = fixed_root_branch(&t).unwrap();
        p.gain[f] = CLAMP;
        p.loss[f] = CLAMP;
        p
    }

    #[test]
    fn em_with_separate_loss_classes_recovers_the_shape_and_never_decreases() {
        // The structure of MALIN's published rate-variation runs: several loss-rate classes, one gain-rate class.
        let t = tree(NWK8);
        let err = ErrorRates {
            sensitivity: vec![0.95; 8],
            false_rate: vec![0.005; 8],
        };
        let mut truth = truth8(0.5);
        truth.rates = Rates::separate(1, 1.0, 3, 0.5);
        let d = sim_data(NWK8, simulate_classes(&t, &truth, &err, 40_000, 17));
        let opt = |it: usize| FitOptions {
            max_iter: it,
            tol: 0.0,
            loss_categories: 3,
            ..Default::default()
        };
        let mut last = f64::NEG_INFINITY;
        for it in [1, 2, 4, 8, 16, 32] {
            let r = fit(&d, &err, &opt(it)).unwrap();
            assert!(
                r.log_lik >= last - 1e-6,
                "iteration {it}: {} < {last}",
                r.log_lik
            );
            last = r.log_lik;
        }
        let r = fit(
            &d,
            &err,
            &FitOptions {
                loss_categories: 3,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(r.converged, "{} iterations", r.iterations);
        let a = r.params.rates.loss_shape;
        assert!((0.25..1.2).contains(&a), "loss shape {a}");
        assert_eq!(
            (r.params.rates.gain_categories, r.params.rates.categories),
            (1, 1)
        );
        assert!(r.log_lik >= log_likelihood(&d, &err, &truth) - 1e-6);
        // far better than without classes, and better than three classes shared by gain and loss (the truth has no gain variation)
        let flat = fit(&d, &err, &FitOptions::default()).unwrap();
        let shared = fit(
            &d,
            &err,
            &FitOptions {
                categories: 3,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            r.log_lik > flat.log_lik + 50.0,
            "{} vs {}",
            r.log_lik,
            flat.log_lik
        );
        assert!(
            r.log_lik > shared.log_lik + 5.0,
            "{} vs shared {}",
            r.log_lik,
            shared.log_lik
        );
        // the posterior loss factor varies between positions and the gain factor does not (one gain class)
        let spread = |x: &[f64]| {
            let m = x.iter().sum::<f64>() / x.len() as f64;
            x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / x.len() as f64
        };
        assert!(
            spread(&r.site_loss_rate) > 0.01,
            "loss factor spread {}",
            spread(&r.site_loss_rate)
        );
        assert!(r.site_gain_rate.iter().all(|&g| (g - 1.0).abs() < 1e-12));
        // stationary in the loss shape and in the branches
        let base = log_likelihood(&d, &err, &r.params);
        assert!((base - r.log_lik).abs() < 1e-6 * base.abs());
        for h in [1.01, 0.99] {
            let mut q = r.params.clone();
            q.rates.loss_shape *= h;
            assert!(
                log_likelihood(&d, &err, &q) <= base + 1e-4,
                "loss shape × {h}"
            );
        }
        for v in 0..t.len() {
            if r.nodes[v].fixed != Fixed::No || t.parent(v).is_none() {
                continue;
            }
            for (which, h) in [(0, 1e-4), (0, -1e-4), (1, 1e-4), (1, -1e-4)] {
                let mut q = r.params.clone();
                if which == 0 {
                    q.gain[v] += h
                } else {
                    q.loss[v] += h
                }
                assert!(
                    log_likelihood(&d, &err, &q) <= base + 1e-4,
                    "node {v} ({which}) by {h}"
                );
            }
        }
    }

    #[test]
    fn em_with_gain_and_loss_classes_estimates_both_shapes_and_the_grid_option_is_checked() {
        let t = tree(NWK8);
        let err = ErrorRates {
            sensitivity: vec![0.95; 8],
            false_rate: vec![0.005; 8],
        };
        let mut truth = truth8(0.5);
        truth.rates = Rates::separate(2, 0.6, 2, 0.6);
        let d = sim_data(NWK8, simulate_classes(&t, &truth, &err, 40_000, 23));
        let r = fit(
            &d,
            &err,
            &FitOptions {
                gain_categories: 2,
                loss_categories: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(r.converged);
        assert!(r.log_lik >= log_likelihood(&d, &err, &truth) - 1e-6);
        assert!(r.params.rates.gain_shape > SHAPE_MIN && r.params.rates.loss_shape > SHAPE_MIN);
        // fixing a shape keeps it, and the other one is still estimated
        let f = fit(
            &d,
            &err,
            &FitOptions {
                gain_categories: 2,
                gain_shape: Some(0.6),
                loss_categories: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(f.params.rates.gain_shape, 0.6);
        assert!(f.params.rates.loss_shape != 1.0);
        // option checks: shared and separate classes are not combined, the grid is capped, shapes are in range
        let bad = |o: FitOptions| fit(&d, &err, &o).unwrap_err();
        assert!(bad(FitOptions {
            categories: 2,
            loss_categories: 2,
            ..Default::default()
        })
        .contains("cannot be combined"));
        assert!(bad(FitOptions {
            gain_categories: 5,
            loss_categories: 5,
            ..Default::default()
        })
        .contains("position classes"));
        assert!(bad(FitOptions {
            loss_categories: 2,
            loss_shape: Some(0.0),
            ..Default::default()
        })
        .contains("shape"));
        assert!(bad(FitOptions {
            gain_categories: 0,
            ..Default::default()
        })
        .contains("gain rate categories"));
    }

    #[test]
    fn em_with_classes_recovers_the_shape_and_never_decreases() {
        let t = tree(NWK8);
        let err = ErrorRates {
            sensitivity: vec![0.95; 8],
            false_rate: vec![0.005; 8],
        };
        let truth = truth8(0.5);
        let d = sim_data(NWK8, simulate_classes(&t, &truth, &err, 40_000, 11));
        let opt = |it: usize| FitOptions {
            max_iter: it,
            tol: 0.0,
            categories: 4,
            ..Default::default()
        };
        let mut last = f64::NEG_INFINITY;
        for it in [1, 2, 4, 8, 16, 32] {
            let r = fit(&d, &err, &opt(it)).unwrap();
            assert!(
                r.log_lik >= last - 1e-6,
                "iteration {it}: {} < {last}",
                r.log_lik
            );
            last = r.log_lik;
        }
        let r = fit(
            &d,
            &err,
            &FitOptions {
                categories: 4,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(r.converged, "{} iterations", r.iterations);
        let a = r.params.rates.shape;
        assert!((0.3..0.8).contains(&a), "shape {a}");
        assert!(r.log_lik >= log_likelihood(&d, &err, &truth) - 1e-6);
        // much better than without classes on heterogeneous data
        let flat = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!(
            r.log_lik > flat.log_lik + 50.0,
            "{} vs {}",
            r.log_lik,
            flat.log_lik
        );
        // fitted rate-1 losses close to the truth; every branch a continuous-time process
        for v in 0..t.len() {
            if t.parent(v).is_some() && r.nodes[v].fixed == Fixed::No {
                assert!(
                    (r.params.loss[v] - truth.loss[v]).abs() < 0.06,
                    "loss {v}: {} vs {}",
                    r.params.loss[v],
                    truth.loss[v]
                );
                assert!(r.params.gain[v] + r.params.loss[v] < 1.0);
            }
        }
        // positions shown by many genomes but not all (changing ones) get higher posterior rates than all-present ones
        let rate_of = |pred: &dyn Fn(&[u8]) -> bool| {
            let (mut s, mut n) = (0.0, 0.0);
            for (c, &x) in d.characters.iter().zip(&r.site_rate) {
                if pred(&c.states) {
                    s += x;
                    n += 1.0;
                }
            }
            s / n
        };
        let all = rate_of(&|c| c.iter().all(|&x| x == PRESENT));
        let mixed = rate_of(&|c| c.iter().filter(|&&x| x == PRESENT).count() == 4);
        assert!(mixed > all, "{mixed} vs {all}");
        // the fit is a stationary point in the shape and the branches
        let base = log_likelihood(&d, &err, &r.params);
        assert!((base - r.log_lik).abs() < 1e-6 * base.abs());
        for h in [1.01, 0.99] {
            let mut q = r.params.clone();
            q.rates.shape *= h;
            assert!(log_likelihood(&d, &err, &q) <= base + 1e-4, "shape × {h}");
        }
        for v in 0..t.len() {
            if r.nodes[v].fixed != Fixed::No || t.parent(v).is_none() {
                continue;
            }
            for (which, h) in [(0, 1e-4), (0, -1e-4), (1, 1e-4), (1, -1e-4)] {
                let mut q = r.params.clone();
                if which == 0 {
                    q.gain[v] += h
                } else {
                    q.loss[v] += h
                }
                assert!(
                    log_likelihood(&d, &err, &q) <= base + 1e-4,
                    "node {v} ({which}) by {h}"
                );
            }
        }
    }

    #[test]
    fn homogeneous_data_give_a_large_shape_and_a_fixed_shape_is_kept() {
        let t = tree(NWK8);
        let err = ErrorRates::none(8);
        let mut truth = truth8(1.0);
        truth.rates.categories = 1;
        let d = sim_data(NWK8, simulate_classes(&t, &truth, &err, 30_000, 3));
        let r = fit(
            &d,
            &err,
            &FitOptions {
                categories: 4,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(r.params.rates.shape > 5.0, "shape {}", r.params.rates.shape);
        let flat = fit(&d, &err, &FitOptions::default()).unwrap();
        assert!(
            (r.log_lik - flat.log_lik).abs() < 5.0,
            "{} vs {}",
            r.log_lik,
            flat.log_lik
        );
        let fixed = fit(
            &d,
            &err,
            &FitOptions {
                categories: 4,
                shape: Some(0.7),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(fixed.params.rates.shape, 0.7);
        assert!(fixed.log_lik < r.log_lik);
    }

    #[test]
    fn bad_rate_options_are_refused() {
        let t = tree(NWK);
        let err = ErrorRates::none(5);
        let d = sim_data(
            NWK,
            simulate(
                &t,
                &truth_with_fixed_root(&t, 0.05, 0.15, 0.7),
                &err,
                500,
                2,
            ),
        );
        for (k, a) in [
            (0, None),
            (MAX_CATEGORIES + 1, None),
            (4, Some(0.01)),
            (4, Some(1e4)),
        ] {
            let o = FitOptions {
                categories: k,
                shape: a,
                ..Default::default()
            };
            assert!(fit(&d, &err, &o).is_err(), "{k} {a:?}");
        }
    }
}
