//! Scoring shared intron positions: Dollo parsimony, maximum parsimony, and counts.
//!
//! Each intron position becomes a discrete binary character ("is there an intron at alignment
//! column *c*, phase *p*?") whose states are observed in the aligned taxa. What differs between
//! methods is the weighting:
//!
//! * **Dollo parsimony** — a position may be gained many times but, once lost, never regained in
//!   that lineage. Gains cost 1, losses cost 0, so the score is the minimum number of independent
//!   origins. A site with one gain is *compatible* with Dollo's law; a site needing two is evidence
//!   either for convergent insertion or against the topology. Implemented as an exact three-state
//!   DP over the tree (`U` never yet gained, `P` present, `L` lost after a gain), which is what
//!   makes the "never regained" constraint bite — a plain Fitch pass cannot express it.
//! * **Maximum parsimony** — symmetric 0↔1 changes cost 1 each (Fitch), the conventional
//!   unrestricted comparison.
//! * **Counts** — no tree at all: how many taxa carry each position, which is the input to
//!   sharing-frequency analyses.
//!
//! Missing data is first-class: a taxon that was never annotated gets `?`, not `0`, because treating
//! "no annotation" as "intron lost" fabricates losses.

use rayon::prelude::*;

use crate::error::{Error, Result};
use crate::matrix::{Cell, SiteMatrix};
use crate::tree::{NodeId, Tree};

/// Tip state codes.
pub const ABSENT: u8 = 0;
pub const PRESENT: u8 = 1;
pub const MISSING: u8 = 2;

/// Internal Dollo states.
const D_NEVER: u8 = 0;
const D_PRESENT: u8 = 1;
const D_LOST: u8 = 2;

/// Which scoring model to apply to the site matrix.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Method {
    Dollo,
    MaxParsimony,
    Counts,
}

impl std::fmt::Display for Method {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Method::Dollo => "dollo",
            Method::MaxParsimony => "mp",
            Method::Counts => "counts",
        })
    }
}

impl Method {
    pub fn parse(s: &str) -> Option<Method> {
        match s.to_ascii_lowercase().as_str() {
            "dollo" | "d" => Some(Method::Dollo),
            "mp" | "parsimony" | "maxparsimony" | "fitch" => Some(Method::MaxParsimony),
            "count" | "counts" | "frequency" => Some(Method::Counts),
            _ => None,
        }
    }
}

/// One scored position.
#[derive(Clone, Debug)]
pub struct Character {
    /// `<column>:<phase>` label from the matrix.
    pub label: String,
    pub col: u32,
    pub phase: u8,
    /// State per tip, in [`ScoringData::tips`] order.
    pub states: Vec<u8>,
}

/// A tree plus the characters prepared against it.
#[derive(Clone, Debug)]
pub struct ScoringData {
    pub tree: Tree,
    /// Tree tips in left-to-right order.
    pub tips: Vec<NodeId>,
    pub tip_names: Vec<String>,
    pub characters: Vec<Character>,
    /// Per-tip note: `None` when the taxon had no row in the matrix (scored as missing).
    pub annotated_tip: Vec<bool>,
}

impl ScoringData {
    /// Position of a tip label in the per-character state vectors.
    pub fn tip_index_of(&self, name: &str) -> usize {
        self.tip_names
            .iter()
            .position(|n| n == name)
            .expect("tip exists in the scoring data")
    }
}

/// What [`prepare`] dropped or renamed, for reporting.
#[derive(Clone, Debug, Default)]
pub struct PrepareReport {
    pub used: Vec<String>,
    /// Matrix taxa absent from the tree: their introns cannot be scored.
    pub taxa_not_in_tree: Vec<String>,
    /// Tree tips absent from the matrix: scored as missing data.
    pub tree_tip_without_data: Vec<String>,
    pub uninformative: usize,
    /// Sites dropped by `min_observed`: too few taxa had an observation.
    pub too_sparse: usize,
}

/// Align matrix taxa with tree tips.
///
/// Names must match exactly — a silent fuzzy match here would place introns on the wrong species,
/// which is precisely the failure mode that produces a convincing-looking but false result.
pub fn prepare(
    m: &SiteMatrix,
    tree: &Tree,
    outgroups: &[String],
    min_taxa: usize,
) -> Result<(ScoringData, PrepareReport)> {
    prepare_with(
        m,
        tree,
        outgroups,
        &PrepareOptions {
            min_taxa,
            min_observed: 0,
        },
    )
}

/// Which sites are scored.
#[derive(Clone, Copy, Debug, Default)]
pub struct PrepareOptions {
    /// Skip sites carried by fewer than this many scored taxa.
    pub min_taxa: usize,
    /// Skip sites observed (present or absent, i.e. not `?`) in fewer than this many scored taxa.
    /// A site that is mostly missing says little about loss; MALIN applies the same kind of filter.
    pub min_observed: usize,
}

/// [`prepare`] with the full set of site filters.
pub fn prepare_with(
    m: &SiteMatrix,
    tree: &Tree,
    outgroups: &[String],
    opts: &PrepareOptions,
) -> Result<(ScoringData, PrepareReport)> {
    let min_taxa = opts.min_taxa;
    let rooted = if outgroups.is_empty() {
        tree.clone()
    } else {
        tree.with_outgroups(outgroups)?
    };
    let tips = rooted.tips();
    let tip_names: Vec<String> = tips
        .iter()
        .map(|&t| rooted.label(t).unwrap_or("?").to_string())
        .collect();

    let mut col_of_taxon: Vec<Option<usize>> = vec![None; tips.len()];
    let mut report = PrepareReport::default();
    for (t, name) in tip_names.iter().enumerate() {
        if let Some(i) = m.taxa.iter().position(|x| x == name) {
            col_of_taxon[t] = Some(i);
            if m.annotated[i] {
                report.used.push(name.clone());
            }
        } else {
            report.tree_tip_without_data.push(name.clone());
        }
    }
    for (i, name) in m.taxa.iter().enumerate() {
        if col_of_taxon.iter().all(|&c| c != Some(i)) && m.annotated[i] {
            report.taxa_not_in_tree.push(name.clone());
        }
    }
    if report.used.is_empty() {
        return Err(Error::NotFound {
            kind: "taxon",
            name: m.taxa.first().cloned().unwrap_or_default(),
            context: format!(
                "no matrix taxon matches a tree tip ({} matrix rows, {} tips) — check naming, or map alignment row ids to tree tips with --taxon-map",
                m.taxa.len(),
                tips.len()
            ),
        });
    }

    let mut characters = Vec::with_capacity(m.sites.len());
    for (si, site) in m.sites.iter().enumerate() {
        let mut states = vec![ABSENT; tips.len()];
        for (t, &maybe) in col_of_taxon.iter().enumerate() {
            states[t] = match maybe {
                Some(mi) => match site.cell(mi) {
                    Cell::Present => PRESENT,
                    Cell::Absent => ABSENT,
                    // Unannotated, or not aligned around this site: no observation.
                    Cell::Missing => MISSING,
                },
                // A tip with no matrix row at all: missing, not absent.
                None => MISSING,
            };
        }
        // `min_taxa` counts carriers that are actually scored: a carrier absent from the tree, or in an
        // unannotated row, must not make the site look shared.
        let n_present = states.iter().filter(|&&s| s == PRESENT).count();
        let n_observed = states.iter().filter(|&&s| s != MISSING).count();
        if n_present == 0 || n_present < min_taxa {
            report.uninformative += 1;
            continue;
        }
        if n_observed < opts.min_observed {
            report.too_sparse += 1;
            continue;
        }
        characters.push(Character {
            label: m.site_label(si),
            col: m.site_column(si),
            phase: site.phase,
            states,
        });
    }
    let annotated_tip = col_of_taxon
        .iter()
        .map(|&c| c.is_some_and(|i| m.annotated[i]))
        .collect();
    Ok((
        ScoringData {
            tree: rooted,
            tips,
            tip_names,
            characters,
            annotated_tip,
        },
        report,
    ))
}

/// Per-site Dollo outcome.
///
/// Classic Dollo parsimony assumes a complex character (here: an intron at an exact
/// position-and-phase) arises **once**, at the most recent common ancestor of the taxa that have it,
/// and is then lost independently any number of times. The informative statistic is therefore the
/// minimum number of independent *losses*; a site needing none is perfectly compatible, and a site
/// needing many is evidence either for convergent insertion or against the topology.
#[derive(Clone, Debug)]
pub struct DolloSite {
    pub label: String,
    pub col: u32,
    pub phase: u8,
    pub n_present: u32,
    pub n_missing: u32,
    /// Always one gain under classic Dollo; zero for an all-absent character.
    pub gains: u32,
    /// Minimum independent losses given that single gain.
    pub losses: u32,
    /// The taxa carrying the site are exactly the leaves of one clade.
    pub monophyletic: bool,
    /// `losses == 0`: one origin, no reversal needed.
    pub compatible: bool,
    /// Where the single origin was placed.
    pub gain_node: Option<String>,
    pub loss_edges: Vec<String>,
    /// Tree node the single origin sits at, and the nodes whose branch above carries a loss. Node ids
    /// are unambiguous where the readable names above are not (two internal branches can print alike).
    pub gain_id: Option<NodeId>,
    pub loss_ids: Vec<NodeId>,
    /// Node states: `0` never gained, `1` present, `2` lost.
    pub reconstruction: Option<Vec<u8>>,
}

/// Dollo gains and losses accumulated on one branch of the tree, over all sites.
#[derive(Clone, Debug)]
pub struct BranchTally {
    pub node: NodeId,
    /// Tip name, non-numeric node label, or `node<id>` for an unlabelled internal branch.
    pub name: String,
    /// Every tip below the branch.
    pub tips: Vec<String>,
    /// Origins placed at this node.
    pub gains: u64,
    /// Losses on the branch leading to this node.
    pub losses: u64,
}

/// Aggregate of a Dollo run.
#[derive(Clone, Debug, Default)]
pub struct DolloResult {
    pub sites: Vec<DolloSite>,
    pub total_gains: u64,
    pub total_losses: u64,
    pub compatible: u64,
    pub incompatible: u64,
    pub monophyletic: u64,
    pub n_tips: usize,
    pub n_characters: usize,
    /// `true` when the Newick looked unrooted and no `--outgroup` was given, so the reported losses
    /// depend on an arbitrary rooting.
    pub rooted_arbitrarily: bool,
    /// Per-branch totals, in tree-node order.
    pub branches: Vec<BranchTally>,
}

/// Score every character under classic Dollo parsimony (one origin, minimum losses).
pub fn score_dollo(data: &ScoringData, reconstructions: bool) -> DolloResult {
    let names = data.tree.node_names();
    let sites: Vec<DolloSite> = data
        .characters
        .par_iter()
        .map(|ch| {
            let r = dollo_one(&data.tree, &data.tips, &ch.states, reconstructions);
            DolloSite {
                label: ch.label.clone(),
                col: ch.col,
                phase: ch.phase,
                n_present: ch.states.iter().filter(|&&s| s == PRESENT).count() as u32,
                n_missing: ch.states.iter().filter(|&&s| s == MISSING).count() as u32,
                gains: if r.gain_node.is_some() { 1 } else { 0 },
                losses: r.losses,
                monophyletic: r.monophyletic,
                compatible: r.losses == 0,
                gain_node: r.gain_id.map(|g| names[g].clone()),
                loss_edges: {
                    let mut v: Vec<String> = r.loss_ids.iter().map(|&l| names[l].clone()).collect();
                    v.sort();
                    v
                },
                gain_id: r.gain_id,
                loss_ids: r.loss_ids,
                reconstruction: r.recon,
            }
        })
        .collect();
    let mut out = DolloResult {
        rooted_arbitrarily: data.tree.looks_unrooted(),
        n_tips: data.tips.len(),
        n_characters: data.characters.len(),
        sites,
        ..Default::default()
    };
    for s in &out.sites {
        out.total_gains += u64::from(s.gains);
        out.total_losses += u64::from(s.losses);
        out.compatible += s.compatible as u64;
        out.monophyletic += s.monophyletic as u64;
        if !s.compatible {
            out.incompatible += 1;
        }
    }
    // Per-branch tallies: the quantity that can be compared with another tool's per-lineage losses.
    let tree = &data.tree;
    let mut gains = vec![0u64; tree.len()];
    let mut losses = vec![0u64; tree.len()];
    for s in &out.sites {
        if let Some(g) = s.gain_id {
            gains[g] += 1;
        }
        for &l in &s.loss_ids {
            losses[l] += 1;
        }
    }
    out.branches = (0..tree.len())
        .map(|v| {
            let tips: Vec<String> = tree
                .leaves_under(v)
                .iter()
                .map(|&t| tree.label(t).unwrap_or("?").to_string())
                .collect();
            // Readable and unique: a tip name, a label, or the tip list (see `Tree::node_names`).
            let name = names[v].clone();
            BranchTally {
                node: v,
                name,
                tips,
                gains: gains[v],
                losses: losses[v],
            }
        })
        .collect();
    out
}

#[derive(Debug)]
struct DolloOne {
    losses: u32,
    monophyletic: bool,
    gain_node: Option<String>,
    gain_id: Option<NodeId>,
    loss_ids: Vec<NodeId>,
    recon: Option<Vec<u8>>,
}

const INF: u32 = 1 << 28;

/// One character under classic Dollo.
///
/// The gain is pinned at the MRCA of the observed carriers — any deeper placement can only add
/// losses, and any shallower one cannot cover all carriers. Below the MRCA a two-state irreversible
/// DP (`present` → `lost`, never back) gives the minimum number of losses.
fn dollo_one(tree: &Tree, tips: &[NodeId], obs: &[u8], keep_recon: bool) -> DolloOne {
    let mut present_tips: Vec<NodeId> = Vec::new();
    for (i, &t) in tips.iter().enumerate() {
        if obs[i] == PRESENT {
            present_tips.push(t);
        }
    }
    if present_tips.is_empty() {
        return DolloOne {
            losses: 0,
            monophyletic: false,
            gain_node: None,
            gain_id: None,
            loss_ids: Vec::new(),
            recon: None,
        };
    }
    let n = tree.len();
    let mut tip_state = vec![MISSING; n];
    for (i, &t) in tips.iter().enumerate() {
        tip_state[t] = obs[i];
    }
    // One postorder pass finds the MRCA (the first node, children before parents, that holds all k
    // carriers is the lowest such node) and how many observed tips lie under it. Monophyly is judged
    // over *observed* tips: a lineage with no data inside the carriers' clade says nothing against the
    // carriers forming that clade.
    let k = present_tips.len() as u32;
    let (mut carriers, mut observed) = (vec![0u32; n], vec![0u32; n]);
    let mut mrca_found: Option<NodeId> = None;
    for &v in tree.postorder() {
        let kids = tree.children(v);
        if kids.is_empty() {
            carriers[v] = u32::from(tip_state[v] == PRESENT);
            observed[v] = u32::from(tip_state[v] != MISSING);
        } else {
            for &c in kids {
                carriers[v] += carriers[c];
                observed[v] += observed[c];
            }
        }
        if mrca_found.is_none() && carriers[v] == k {
            mrca_found = Some(v);
        }
    }
    let mrca = mrca_found.expect("the root holds every carrier");
    let monophyletic = observed[mrca] == k;

    // Irreversible two-state DP over the subtree below the gain.
    let mut cp = vec![0u32; n];
    let mut cl = vec![0u32; n];
    let mut in_subtree = vec![false; n];
    {
        let mut stack = vec![mrca];
        while let Some(v) = stack.pop() {
            in_subtree[v] = true;
            stack.extend(tree.children(v).iter().copied());
        }
    }
    for &v in tree.postorder() {
        if !in_subtree[v] {
            continue;
        }
        let kids = tree.children(v);
        if kids.is_empty() {
            match tip_state[v] {
                PRESENT => {
                    cp[v] = 0;
                    cl[v] = INF;
                }
                ABSENT => {
                    cp[v] = INF;
                    cl[v] = 0;
                }
                _ => {
                    cp[v] = 0;
                    cl[v] = 0;
                }
            }
            continue;
        }
        let (mut sp, mut sl) = (0u32, 0u32);
        for &c in kids {
            sp = sp.saturating_add(cp[c].min(cl[c].saturating_add(1)));
            sl = sl.saturating_add(cl[c]);
        }
        cp[v] = sp;
        cl[v] = sl;
    }

    // Traceback from the gain node downward.
    let mut recon = keep_recon.then(|| vec![D_NEVER; n]);
    let mut loss_ids: Vec<NodeId> = Vec::new();
    let mut stack = vec![(mrca, true)]; // (node, is present)
    while let Some((v, present)) = stack.pop() {
        if let Some(r) = recon.as_mut() {
            r[v] = if present { D_PRESENT } else { D_LOST };
        }
        for &c in tree.children(v) {
            let child_present = present && cp[c] <= cl[c].saturating_add(1) && cp[c] < INF;
            if present && !child_present {
                loss_ids.push(c);
            }
            stack.push((c, child_present));
        }
    }
    loss_ids.sort_unstable();
    DolloOne {
        losses: cp[mrca].min(INF),
        monophyletic,
        // Readable names are filled in once per tree by `score_dollo` (from ids), not per site.
        gain_node: Some(String::new()),
        gain_id: Some(mrca),
        loss_ids,
        recon,
    }
}

/// Sensitivity analysis: the weighted irreversible model, where a gain costs `gain_weight` and a
/// loss `loss_weight`, and a lost position can never come back.
///
/// With `gain_weight = 1, loss_weight = 0` every site is explained by one gain at the root, so the
/// interesting region is `loss_weight > 0` (e.g. 1:1 trades an extra origin against the losses it
/// would avoid). Returns `(site_label, gains, losses)` per site.
pub fn dollo_weighted(
    data: &ScoringData,
    gain_weight: u32,
    loss_weight: u32,
) -> Vec<(String, u32, u32)> {
    data.characters
        .par_iter()
        .map(|ch| {
            let (g, l) = weighted_one(&data.tree, &data.tips, &ch.states, gain_weight, loss_weight);
            (ch.label.clone(), g, l)
        })
        .collect()
}

/// Three-state weighted DP: `U` not yet gained, `P` present, `L` lost (terminal).
fn weighted_one(tree: &Tree, tips: &[NodeId], obs: &[u8], gw: u32, lw: u32) -> (u32, u32) {
    let n = tree.len();
    let mut u = vec![0u64; n];
    let mut p = vec![0u64; n];
    let mut l = vec![0u64; n];
    let mut tip_state = vec![MISSING; n];
    for (i, &t) in tips.iter().enumerate() {
        tip_state[t] = obs[i];
    }
    for &v in tree.postorder() {
        let kids = tree.children(v);
        if kids.is_empty() {
            match tip_state[v] {
                PRESENT => {
                    u[v] = INF_W;
                    p[v] = 0;
                    l[v] = INF_W;
                }
                ABSENT => {
                    u[v] = 0;
                    p[v] = INF_W;
                    l[v] = 0;
                }
                _ => {}
            }
            continue;
        }
        let (mut su, mut sp, mut sl) = (0u64, 0u64, 0u64);
        for &c in kids {
            su = su.saturating_add(u[c].min(wadd(p[c], gw)));
            sp = sp.saturating_add(p[c].min(wadd(l[c], lw)));
            sl = sl.saturating_add(l[c]);
        }
        u[v] = su;
        p[v] = sp;
        l[v] = sl;
    }
    let root = tree.root();
    // The cheaper root state decides whether a gain sits on the root itself; either way the
    // traceback below re-walks the same recurrences, so gains*gw + losses*lw equals the DP cost.
    let use_present = wadd(p[root], gw) < u[root];
    let mut gains = 0u32;
    let mut losses = 0u32;
    let mut state = vec![D_NEVER; n];
    if use_present {
        state[root] = D_PRESENT;
        gains += 1;
    }
    let mut stack = vec![root];
    while let Some(v) = stack.pop() {
        for &c in tree.children(v) {
            let st = match state[v] {
                D_NEVER => {
                    if wadd(p[c], gw) < u[c] {
                        gains += 1;
                        D_PRESENT
                    } else {
                        D_NEVER
                    }
                }
                D_PRESENT => {
                    if wadd(l[c], lw) < p[c] {
                        losses += 1;
                        D_LOST
                    } else {
                        D_PRESENT
                    }
                }
                _ => D_LOST,
            };
            state[c] = st;
            stack.push(c);
        }
    }
    (gains, losses)
}

/// "Impossible" cost for the weighted DP. 64-bit so that any `u32` weight, summed over a tree, stays
/// far below it: a gain costing more than the sentinel would silently read as impossible.
const INF_W: u64 = 1 << 60;

fn wadd(x: u64, w: u32) -> u64 {
    if x >= INF_W {
        INF_W
    } else {
        (x + u64::from(w)).min(INF_W - 1)
    }
}

/// Per-site unrestricted maximum-parsimony outcome.
#[derive(Clone, Debug)]
pub struct MpSite {
    pub label: String,
    pub col: u32,
    pub phase: u8,
    pub steps: u32,
    pub n_present: u32,
    pub n_missing: u32,
    /// Informative for the tree: at least two taxa on each side.
    pub informative: bool,
    /// Ambiguous sites have more than one equally short reconstruction.
    pub ambiguous: bool,
    pub reconstruction: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Default)]
pub struct MpResult {
    pub sites: Vec<MpSite>,
    pub total_steps: u64,
    pub informative: u64,
    pub ambiguous: u64,
    /// Consistency index over the scored characters: informative sites / total steps.
    pub ci: f64,
}

/// Fitch parsimony: unit cost for every state change, gains and losses symmetric.
pub fn score_mp(data: &ScoringData, reconstructions: bool) -> MpResult {
    let sites: Vec<MpSite> = data
        .characters
        .par_iter()
        .map(|ch| {
            let (steps, ambiguous, recon) =
                fitch_one(&data.tree, &data.tips, &ch.states, reconstructions);
            let n_present = ch.states.iter().filter(|&&s| s == PRESENT).count() as u32;
            let n_absent = ch.states.iter().filter(|&&s| s == ABSENT).count() as u32;
            MpSite {
                label: ch.label.clone(),
                col: ch.col,
                phase: ch.phase,
                steps,
                n_present,
                n_missing: ch.states.iter().filter(|&&s| s == MISSING).count() as u32,
                informative: n_present >= 2 && n_absent >= 2,
                ambiguous,
                reconstruction: recon,
            }
        })
        .collect();
    let mut out = MpResult {
        sites,
        ..Default::default()
    };
    for s in &out.sites {
        out.total_steps += u64::from(s.steps);
        out.informative += s.informative as u64;
        out.ambiguous += s.ambiguous as u64;
    }
    out.ci = if out.total_steps == 0 {
        0.0
    } else {
        out.informative as f64 / out.total_steps as f64
    };
    out
}

/// Fitch down-pass/up-pass on one binary character.
///
/// Tip states are `0`, `1` or `2` (missing). A node whose whole subtree is missing carries no
/// information and contributes nothing to the step count.
fn fitch_one(
    tree: &Tree,
    tips: &[NodeId],
    obs: &[u8],
    keep_recon: bool,
) -> (u32, bool, Option<Vec<u8>>) {
    let n = tree.len();
    // State sets are bitmasks: bit0 = "can be absent", bit1 = "can be present".
    let mut set = vec![0b11u8; n];
    let mut no_info = vec![false; n];
    let mut tip_state = vec![MISSING; n];
    for (i, &t) in tips.iter().enumerate() {
        tip_state[t] = obs[i];
    }

    let mut steps: u32 = 0;
    for &v in tree.postorder() {
        let kids = tree.children(v);
        if kids.is_empty() {
            match tip_state[v] {
                PRESENT => set[v] = 0b10,
                ABSENT => set[v] = 0b01,
                _ => {
                    set[v] = 0b11;
                    no_info[v] = true;
                }
            }
            continue;
        }
        // Hartigan's n-ary generalisation of Fitch: for each state count the informed children whose
        // set allows it; the node takes the state(s) with the maximum count, and every informed child
        // outside that count costs one step. Pairwise intersection is only right for binary nodes —
        // children {0},{0},{1},{1},{1} need 2 steps, not 1 — and unrooted trees have a polytomous root.
        let mut cnt = [0u32; 2];
        let mut informed = 0u32;
        for &c in kids {
            if no_info[c] {
                continue;
            }
            informed += 1;
            cnt[0] += u32::from(set[c] & 0b01 != 0);
            cnt[1] += u32::from(set[c] & 0b10 != 0);
        }
        no_info[v] = informed == 0;
        if informed == 0 {
            set[v] = 0b11;
        } else {
            let best = cnt[0].max(cnt[1]);
            steps += informed - best;
            set[v] = u8::from(cnt[0] == best) | (u8::from(cnt[1] == best) << 1);
        }
    }

    let ambiguous = has_several_optimal_labellings(tree, tips, obs);
    if !keep_recon {
        return (steps, ambiguous, None);
    }

    // Up-pass: keep the parent's state where the down-pass allows it, so the result is minimal.
    let mut recon = vec![MISSING; n];
    let root = tree.root();
    if !no_info[root] {
        // Prefer absence on a tie: a single-tip present state stays present, and the total is unchanged.
        recon[root] = if set[root] & 0b01 != 0 {
            ABSENT
        } else {
            PRESENT
        };
    }
    let mut stack = vec![root];
    while let Some(v) = stack.pop() {
        for &c in tree.children(v) {
            recon[c] = if no_info[c] {
                MISSING
            } else if recon[v] != MISSING && set[c] & bit(recon[v]) != 0 {
                recon[v]
            } else if set[c] & 0b01 != 0 {
                ABSENT
            } else {
                PRESENT
            };
            stack.push(c);
        }
    }
    (steps, ambiguous, Some(recon))
}

/// `true` when the site has more than one most-parsimonious labelling of the internal nodes.
///
/// A node whose Fitch set is `{0,1}` is *not* enough: `((A,B),(C,D))` with `1110` has such a node but a
/// unique optimal reconstruction (everything present, one loss on D). So the optimal labellings are
/// counted exactly: per node and state, the cheapest cost below it and the number of labellings that
/// achieve it (saturating at 2, since only "one or more than one" matters). Unobserved tips and
/// subtrees without any observed tip are free and add no labelling choice.
fn has_several_optimal_labellings(tree: &Tree, tips: &[NodeId], obs: &[u8]) -> bool {
    const BIG: u32 = 1 << 20;
    let n = tree.len();
    let mut tip_state = vec![MISSING; n];
    for (i, &t) in tips.iter().enumerate() {
        tip_state[t] = obs[i];
    }
    // cost[v][s], ways[v][s]; `informed[v]` = some observed tip lies below v.
    let mut cost = vec![[0u32; 2]; n];
    let mut ways = vec![[1u32; 2]; n];
    let mut informed = vec![false; n];
    for &v in tree.postorder() {
        let kids = tree.children(v);
        if kids.is_empty() {
            match tip_state[v] {
                PRESENT => {
                    informed[v] = true;
                    cost[v] = [BIG, 0];
                    ways[v] = [0, 1];
                }
                ABSENT => {
                    informed[v] = true;
                    cost[v] = [0, BIG];
                    ways[v] = [1, 0];
                }
                _ => {}
            }
            continue;
        }
        for s in 0..2usize {
            let (mut c, mut w) = (0u32, 1u32);
            for &k in kids.iter().filter(|&&k| informed[k]) {
                let via = |y: usize| cost[k][y].saturating_add(u32::from(y != s));
                let best = via(0).min(via(1));
                let ww = (0..2usize)
                    .filter(|&y| via(y) == best)
                    .fold(0u32, |acc, y| (acc + ways[k][y]).min(2));
                c = c.saturating_add(best);
                w = (w * ww).min(2);
            }
            cost[v][s] = c;
            ways[v][s] = w;
        }
        informed[v] = kids.iter().any(|&k| informed[k]);
    }
    let root = tree.root();
    if !informed[root] {
        return false;
    }
    let best = cost[root][0].min(cost[root][1]);
    let total: u32 = (0..2)
        .filter(|&s| cost[root][s] == best)
        .map(|s| ways[root][s])
        .sum();
    total > 1
}

fn bit(state: u8) -> u8 {
    1u8 << state
}

/// Tree-free frequency table for each position.
#[derive(Clone, Debug)]
pub struct CountsSite {
    pub label: String,
    pub col: u32,
    pub phase: u8,
    pub n_present: u32,
    pub n_absent: u32,
    /// Taxa with no observation at this site (unannotated, or not aligned around it).
    pub n_missing: u32,
    /// `n_present / (n_present + n_absent)`: missing cells are excluded, not counted as absent.
    pub frequency: f64,
    pub taxa: Vec<String>,
}

pub fn score_counts(m: &SiteMatrix) -> Vec<CountsSite> {
    (0..m.sites.len())
        .map(|i| {
            let site = &m.sites[i];
            let n_present = site.n_taxa();
            // Presence is a subset of observability by construction, so this cannot underflow.
            let n_absent = site.n_observable() - n_present;
            let n_missing = m.taxa.len() as u32 - site.n_observable();
            CountsSite {
                label: m.site_label(i),
                col: m.site_column(i),
                phase: site.phase,
                n_present,
                n_absent,
                n_missing,
                frequency: if n_present + n_absent == 0 {
                    0.0
                } else {
                    n_present as f64 / (n_present + n_absent) as f64
                },
                taxa: m.taxa_names_at(i).iter().map(|s| s.to_string()).collect(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::ProteinSite;
    use crate::project::{ColumnSite, ProjectionReport, TaxonProjection};

    fn tree(s: &str) -> Tree {
        Tree::from_newick(s).unwrap()
    }

    /// Build a matrix from `(label, column, phase, taxa_present)` and prepare it against a tree.
    fn data(newick: &str, sites: &[(&str, u32, u8, &[&str])], annotated: &[&str]) -> ScoringData {
        let m = build_matrix(&tree(newick), sites, annotated);
        let (d, _) = prepare(&m, &tree(newick), &[], 1).unwrap();
        d
    }

    fn build_matrix(
        t: &Tree,
        sites: &[(&str, u32, u8, &[&str])],
        annotated: &[&str],
    ) -> SiteMatrix {
        let names: Vec<String> = t.tip_labels();
        let mut ap = crate::project::AlignmentProjection {
            alignment: "t".into(),
            alen: 40,
            ..Default::default()
        };
        for (i, n) in names.iter().enumerate() {
            ap.taxa.push(TaxonProjection {
                taxon: i,
                row_id: n.clone(),
                taxon_name: n.clone(),
                protein_idx: Some(i),
                sites: Vec::new(),
                report: ProjectionReport::default(),
                residues_aligned: 10,
                has_annotation: annotated.is_empty() || annotated.contains(&n.as_str()),
                unannotated_reason: None,
                in_sequence_set: true,
                zero_introns: false,
            });
        }
        for &(_, col, phase, taxa) in sites {
            for tx in taxa {
                let ti = names.iter().position(|x| x == tx).unwrap();
                ap.sites.push(ColumnSite {
                    taxon: ti as u32,
                    col,
                    phase,
                    site: ProteinSite::from_cds_offset(col * 3 + phase as u32),
                    clamped: false,
                });
            }
        }
        ap.sites.sort_unstable_by_key(|c| (c.col, c.phase, c.taxon));
        SiteMatrix::build(&ap)
    }

    /// Exhaustive check of classic Dollo: try every possible origin node, and for each enumerate all
    /// irreversible present/lost labellings below it, minimising losses.
    fn brute_force_losses(tree: &Tree, tips: &[NodeId], obs: &[u8]) -> u32 {
        let present: Vec<NodeId> = tips
            .iter()
            .copied()
            .zip(obs)
            .filter(|&(_, &o)| o == PRESENT)
            .map(|(t, _)| t)
            .collect();
        if present.is_empty() {
            return 0;
        }
        let mut best = u32::MAX;
        for candidate in 0..tree.len() {
            if !present.iter().all(|&t| is_under(tree, t, candidate)) {
                continue; // an origin here cannot explain every carrier
            }
            let mut assign = vec![0u8; tree.len()];
            assign[candidate] = 1;
            best = best.min(bf_losses(tree, candidate, &mut assign, tips, obs));
        }
        best
    }

    fn is_under(tree: &Tree, mut n: NodeId, ancestor: NodeId) -> bool {
        loop {
            if n == ancestor {
                return true;
            }
            match tree.parent(n) {
                Some(p) => n = p,
                None => return false,
            }
        }
    }

    fn bf_losses(tree: &Tree, v: NodeId, assign: &mut [u8], tips: &[NodeId], obs: &[u8]) -> u32 {
        if tree.children(v).is_empty() {
            let observed = tips.iter().position(|&t| t == v).map(|i| obs[i]);
            return match observed {
                Some(PRESENT) if assign[v] != 1 => u32::MAX,
                Some(ABSENT) if assign[v] == 1 => u32::MAX,
                _ => 0,
            };
        }
        let mut total = 0u32;
        let mut options: Vec<u32> = Vec::new();
        for &c in tree.children(v) {
            // Child present (no loss on this edge) or child lost (one loss, and then all its
            // descendants are lost too).
            let mut keep: Vec<u8> = assign.to_vec();
            keep[c] = assign[v];
            let cost_keep = if assign[v] == 1 {
                bf_losses(tree, c, &mut keep, tips, obs)
            } else {
                u32::MAX
            };
            let mut lost: Vec<u8> = assign.to_vec();
            lost[c] = 0;
            let cost_lost = if assign[v] == 1 {
                1u32.saturating_add(bf_absent(tree, c, &mut lost, tips, obs))
            } else {
                u32::MAX
            };
            options.push(cost_keep.min(cost_lost));
        }
        total += options.iter().sum::<u32>();
        if total == u32::MAX {
            return u32::MAX;
        }
        total
    }

    /// Once a lineage is lost, nothing below it can be present (Dollo irreversibility).
    fn bf_absent(tree: &Tree, v: NodeId, assign: &mut [u8], tips: &[NodeId], obs: &[u8]) -> u32 {
        assign[v] = 0;
        if let Some(i) = tips.iter().position(|&t| t == v) {
            if obs[i] == PRESENT {
                return u32::MAX;
            }
        }
        let mut total = 0;
        for &c in tree.children(v) {
            let sub = bf_absent(tree, c, assign, tips, obs);
            if sub == u32::MAX {
                return u32::MAX;
            }
            total += sub;
        }
        total
    }

    #[test]
    fn a_missing_tip_inside_the_carriers_clade_does_not_break_monophyly() {
        let t = tree("(((A,B),C),D);");
        let tips = t.tips();
        // A and C carry; B (between them) has no data; D observed absent.
        let obs = [PRESENT, MISSING, PRESENT, ABSENT];
        let r = dollo_one(&t, &tips, &obs, false);
        assert!(
            r.monophyletic,
            "B is unobserved, so {{A,C}} is the observed clade"
        );
        assert_eq!(r.losses, 0);
    }

    #[test]
    fn monophyletic_site_costs_one_gain_and_no_losses() {
        let d = data(
            "(((A,B),(C,D)),E);",
            &[("10:0", 10, 0, &["A", "B"])],
            &["A", "B", "C", "D", "E"],
        );
        let r = score_dollo(&d, false);
        assert_eq!(r.sites.len(), 1);
        let s = &r.sites[0];
        assert_eq!((s.gains, s.losses), (1, 0));
        assert!(s.compatible && s.monophyletic);
        assert_eq!(s.gain_node.as_deref(), Some("A,B"));
    }

    #[test]
    fn non_sister_carriers_imply_one_deep_gain_and_two_losses() {
        // A and C are in different subclades: one origin at their MRCA, then independent losses on
        // the branches to B and D. Under classic Dollo this costs 1 gain + 2 losses.
        let d = data(
            "(((A,B),(C,D)),E);",
            &[("7:0", 7, 0, &["A", "C"])],
            &["A", "B", "C", "D", "E"],
        );
        let s = &score_dollo(&d, false).sites[0];
        assert_eq!((s.gains, s.losses), (1, 2));
        assert!(!s.monophyletic && !s.compatible);
        assert_eq!(s.gain_node.as_deref(), Some("A,B,C,D"));
        assert_eq!(s.loss_edges, vec!["B".to_string(), "D".to_string()]);
        // Unrestricted MP explains the same pattern with two symmetric changes and no losses.
        let mp = score_mp(&d, false);
        assert_eq!(
            mp.sites[0].steps, 2,
            "two independent origins cost the same as one gain plus two losses"
        );
    }

    #[test]
    fn a_widespread_site_originates_at_the_root() {
        let d = data(
            "(((A,B),(C,D)),E);",
            &[("3:1", 3, 1, &["A", "B", "C", "D", "E"])],
            &["A", "B", "C", "D", "E"],
        );
        let s = &score_dollo(&d, false).sites[0];
        assert_eq!((s.gains, s.losses), (1, 0));
        assert!(
            s.gain_node.as_deref().unwrap().contains("root"),
            "{:?}",
            s.gain_node
        );
    }

    #[test]
    fn three_islands_share_one_deep_gain_and_cost_three_losses() {
        // Losses are free under Dollo, so a single early origin still wins over three origins; this
        // is the property that makes the "one gain, many losses" score meaningful at all.
        let d = data(
            "(((A,B),(C,D)),(E,F));",
            &[("9:0", 9, 0, &["A", "C", "E"])],
            &["A", "B", "C", "D", "E", "F"],
        );
        let s = &score_dollo(&d, false).sites[0];
        assert_eq!((s.gains, s.losses), (1, 3));
        assert_eq!(
            s.loss_edges,
            vec!["B".to_string(), "D".to_string(), "F".to_string()]
        );
    }

    #[test]
    fn missing_data_does_not_create_losses() {
        // D sits inside the gained clade. Observed absent it costs a loss; unannotated it costs none.
        let sites = &[("4:0", 4u32, 0u8, &["A", "B", "C"][..])];
        let with_d = data("(((A,B),(C,D)),E);", sites, &["A", "B", "C", "D", "E"]);
        let without_d = data("(((A,B),(C,D)),E);", sites, &["A", "B", "C", "E"]);
        assert_eq!(
            score_dollo(&with_d, false).sites[0].losses,
            1,
            "D lost the intron"
        );
        let s = &score_dollo(&without_d, false).sites[0];
        assert_eq!(s.losses, 0, "D is '?', so no loss is inferred");
        assert_eq!(s.n_missing, 1);
        assert!(s.compatible);
    }

    #[test]
    fn dp_matches_brute_force_over_every_pattern() {
        for nw in [
            "((A,B),C);",
            "(((A,B),C),D);",
            "(((A,B),(C,D)),E);",
            "((((A,B),C),(D,E)),F);",
            "((A,(B,(C,D))),E);",
        ] {
            let t = tree(nw);
            let tips = t.tips();
            for mask in 0u32..(1 << tips.len()) {
                let obs: Vec<u8> = (0..tips.len())
                    .map(|i| if mask >> i & 1 == 1 { PRESENT } else { ABSENT })
                    .collect();
                let got = dollo_one(&t, &tips, &obs, false).losses;
                let want = brute_force_losses(&t, &tips, &obs);
                assert_eq!(
                    got, want,
                    "{nw} obs {obs:?}: DP {got} losses vs brute force {want}"
                );
            }
        }
    }

    /// Exhaustive minimum number of state changes over every labelling of the internal nodes (and of
    /// missing tips, which are free), for a symmetric 0/1 character.
    fn brute_force_mp(tree: &Tree, tips: &[NodeId], obs: &[u8]) -> u32 {
        let n = tree.len();
        let mut free: Vec<NodeId> = (0..n).filter(|&v| !tree.children(v).is_empty()).collect();
        for (i, &t) in tips.iter().enumerate() {
            if obs[i] == MISSING {
                free.push(t);
            }
        }
        let mut best = u32::MAX;
        for mask in 0u32..(1 << free.len()) {
            let mut st = vec![0u8; n];
            for (i, &t) in tips.iter().enumerate() {
                if obs[i] != MISSING {
                    st[t] = obs[i];
                }
            }
            for (k, &v) in free.iter().enumerate() {
                st[v] = (mask >> k & 1) as u8;
            }
            let mut cost = 0u32;
            for v in 0..n {
                for &c in tree.children(v) {
                    cost += u32::from(st[v] != st[c]);
                }
            }
            best = best.min(cost);
        }
        best
    }

    const SHAPES: [&str; 8] = [
        "((A,B),C);",
        "(A,B,C,D,E);",
        "((A,B,C),D);",
        "((A,B,C),(D,E));",
        "(((A,B),C,D),E,F);",
        "((A,B),(C,D),(E,F));",
        "((A,(B,(C,D))),E);",
        "(A,(B,C,D),(E,F));",
    ];

    fn patterns(n: usize) -> impl Iterator<Item = Vec<u8>> {
        (0..3u32.pow(n as u32)).map(move |mut m| {
            (0..n)
                .map(|_| {
                    let d = (m % 3) as u8;
                    m /= 3;
                    d
                })
                .collect()
        })
    }

    #[test]
    fn fitch_matches_brute_force_with_polytomies_and_missing_tips() {
        for nw in SHAPES {
            let t = tree(nw);
            let tips = t.tips();
            for obs in patterns(tips.len()) {
                let (steps, _, recon) = fitch_one(&t, &tips, &obs, true);
                let want = brute_force_mp(&t, &tips, &obs);
                assert_eq!(
                    steps, want,
                    "{nw} obs {obs:?}: Fitch {steps} vs brute force {want}"
                );
                // The reconstruction must realise exactly the reported number of steps.
                let recon = recon.unwrap();
                let mut realised = 0;
                for v in 0..t.len() {
                    for &c in t.children(v) {
                        if recon[v] != MISSING && recon[c] != MISSING && recon[v] != recon[c] {
                            realised += 1;
                        }
                    }
                }
                assert_eq!(realised, steps, "{nw} obs {obs:?}: reconstruction steps");
            }
        }
    }

    #[test]
    fn dollo_matches_brute_force_with_polytomies_and_missing_tips() {
        for nw in SHAPES {
            let t = tree(nw);
            let tips = t.tips();
            for obs in patterns(tips.len()) {
                let got = dollo_one(&t, &tips, &obs, false).losses;
                let want = brute_force_losses(&t, &tips, &obs);
                assert_eq!(
                    got, want,
                    "{nw} obs {obs:?}: DP {got} vs brute force {want}"
                );
            }
        }
    }

    #[test]
    fn branch_tallies_add_up_to_the_site_totals_and_name_internal_branches_unambiguously() {
        let d = data(
            "(((A,B),C),((D,E),F));",
            &[
                ("1:0", 1, 0, &["A", "B", "D"]),
                ("2:0", 2, 0, &["A", "B", "C", "D", "E"]),
                ("3:1", 3, 1, &["D", "E"]),
            ],
            &["A", "B", "C", "D", "E", "F"],
        );
        let r = score_dollo(&d, false);
        assert_eq!(
            r.branches.iter().map(|b| b.gains).sum::<u64>(),
            r.total_gains
        );
        assert_eq!(
            r.branches.iter().map(|b| b.losses).sum::<u64>(),
            r.total_losses
        );
        // Two internal branches can print the same truncated tip list; ids and full tip lists cannot.
        let mut names: Vec<&str> = r.branches.iter().map(|b| b.name.as_str()).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), r.branches.len(), "branch names must be unique");
        let ab = r.branches.iter().find(|b| b.tips == ["A", "B"]).unwrap();
        assert_eq!(ab.losses, 0);
    }

    #[test]
    fn unrooted_polytomy_is_scored_with_the_true_minimum() {
        // (A,B,C,D,E) with 0,0,1,1,1 needs two steps, not one.
        let t = tree("(A,B,C,D,E);");
        let tips = t.tips();
        let (steps, _, _) = fitch_one(
            &t,
            &tips,
            &[ABSENT, ABSENT, PRESENT, PRESENT, PRESENT],
            false,
        );
        assert_eq!(steps, 2);
    }

    #[test]
    fn weighted_dollo_does_not_lose_a_gain_at_huge_weights() {
        let t = tree("((A,B),C);");
        let tips = t.tips();
        let (g, l) = weighted_one(&t, &tips, &[PRESENT, ABSENT, ABSENT], 3_000_000_000, 1);
        assert_eq!(
            (g, l),
            (1, 0),
            "one gain on the carrier's own branch; it must not read as impossible"
        );
    }

    #[test]
    fn reconstruction_is_a_legal_dollo_scenario() {
        let d = data(
            "((((A,B),(C,D)),(E,F)),G);",
            &[
                ("1:0", 1, 0, &["A", "C", "F"]),
                ("2:1", 2, 1, &["A", "B", "C", "D"]),
            ],
            &["A", "B", "C", "D", "E", "F", "G"],
        );
        let r = score_dollo(&d, true);
        for site in &r.sites {
            let recon = site.reconstruction.as_ref().unwrap();
            let mut gains = 0;
            let mut losses = 0;
            for (p, c) in d.tree.edges() {
                match (recon[p], recon[c]) {
                    (D_NEVER, D_PRESENT) => gains += 1,
                    (D_PRESENT, D_LOST) => losses += 1,
                    (D_LOST, D_PRESENT) => panic!("a lost intron came back on site {}", site.label),
                    (D_NEVER, D_LOST) => panic!("loss without a prior gain on site {}", site.label),
                    _ => {}
                }
            }
            assert_eq!(gains, site.gains, "site {}", site.label);
            assert_eq!(losses, site.losses, "site {}", site.label);
            for (i, &t) in d.tips.iter().enumerate() {
                let ch = &d.characters[d
                    .characters
                    .iter()
                    .position(|c| c.label == site.label)
                    .unwrap()];
                match ch.states[i] {
                    PRESENT => {
                        assert_eq!(recon[t], D_PRESENT, "{} at {}", site.label, d.tip_names[i])
                    }
                    ABSENT => {
                        assert_ne!(recon[t], D_PRESENT, "{} at {}", site.label, d.tip_names[i])
                    }
                    _ => {}
                }
            }
        }
    }

    #[test]
    fn weighted_dollo_trades_origins_against_losses() {
        let d = data(
            "(((A,B),(C,D)),E);",
            &[("7:0", 7, 0, &["A", "C"])],
            &["A", "B", "C", "D", "E"],
        );
        // Equal weights: two cheap origins beat one origin plus two losses.
        assert_eq!(dollo_weighted(&d, 1, 1), vec![("7:0".to_string(), 2, 0)]);
        // Expensive origins: the deep single gain with two losses now wins.
        assert_eq!(dollo_weighted(&d, 3, 1), vec![("7:0".to_string(), 1, 2)]);
        // A monophyletic site never needs more than one origin, whatever the weights.
        let mono = data(
            "(((A,B),(C,D)),E);",
            &[("7:0", 7, 0, &["A", "B"])],
            &["A", "B", "C", "D", "E"],
        );
        assert_eq!(dollo_weighted(&mono, 1, 1), vec![("7:0".to_string(), 1, 0)]);
        assert_eq!(dollo_weighted(&mono, 5, 1), vec![("7:0".to_string(), 1, 0)]);
    }

    #[test]
    fn mp_counts_and_totals() {
        let d = data(
            "(((A,B),(C,D)),E);",
            &[
                ("1:0", 1, 0, &["A", "B"]),
                ("2:0", 2, 0, &["A", "C"]),
                ("3:2", 3, 2, &["A", "B", "C", "D"]),
            ],
            &["A", "B", "C", "D", "E"],
        );
        let mp = score_mp(&d, false);
        assert_eq!(mp.sites[0].steps, 1, "one origin for a sister pair");
        assert_eq!(mp.sites[1].steps, 2);
        assert_eq!(mp.sites[2].steps, 1, "present everywhere but E: one change");
        assert_eq!(mp.total_steps, 4);
        assert_eq!(
            mp.informative, 2,
            "the last site is 4-1, not parsimony-informative"
        );
        let dollo = score_dollo(&d, false);
        assert_eq!((dollo.total_gains, dollo.total_losses), (3, 2));
        assert_eq!(dollo.compatible, 2);
        assert_eq!(dollo.incompatible, 1);
        assert_eq!(dollo.n_tips, 5);
    }

    #[test]
    fn mp_reconstruction_realises_its_own_length() {
        let d = data(
            "(((A,B),(C,D)),E);",
            &[("1:0", 1, 0, &["A", "D"])],
            &["A", "B", "C", "D", "E"],
        );
        let site = &score_mp(&d, true).sites[0];
        let recon = site.reconstruction.as_ref().unwrap();
        let steps = d
            .tree
            .edges()
            .filter(|&(p, c)| recon[p] != recon[c])
            .count() as u32;
        assert_eq!(
            steps, site.steps,
            "reconstruction must realise the reported length"
        );
        for (i, &t) in d.tips.iter().enumerate() {
            if d.characters[0].states[i] != MISSING {
                assert_eq!(
                    recon[t], d.characters[0].states[i],
                    "tips keep their observed state"
                );
            }
        }
    }

    #[test]
    fn min_taxa_filters_autapomorphies() {
        let t = tree("((A,B),C);");
        let mut ap = crate::project::AlignmentProjection {
            alignment: "t".into(),
            ..Default::default()
        };
        for (i, n) in ["A", "B", "C"].iter().enumerate() {
            ap.taxa.push(TaxonProjection {
                taxon: i,
                row_id: n.to_string(),
                taxon_name: n.to_string(),
                protein_idx: Some(i),
                sites: Vec::new(),
                report: ProjectionReport::default(),
                residues_aligned: 5,
                has_annotation: true,
                unannotated_reason: None,
                in_sequence_set: true,
                zero_introns: false,
            });
        }
        ap.sites.push(ColumnSite {
            taxon: 0,
            col: 1,
            phase: 0,
            site: ProteinSite::from_cds_offset(3),
            clamped: false,
        });
        ap.sites.push(ColumnSite {
            taxon: 0,
            col: 2,
            phase: 0,
            site: ProteinSite::from_cds_offset(6),
            clamped: false,
        });
        ap.sites.push(ColumnSite {
            taxon: 1,
            col: 2,
            phase: 0,
            site: ProteinSite::from_cds_offset(6),
            clamped: false,
        });
        let m = SiteMatrix::build(&ap);
        let (d, rep) = prepare(&m, &t, &[], 2).unwrap();
        assert_eq!(
            d.characters.len(),
            1,
            "only the shared site survives min_taxa=2"
        );
        assert_eq!(rep.uninformative, 1);
        assert_eq!(prepare(&m, &t, &[], 1).unwrap().0.characters.len(), 2);
    }

    #[test]
    fn rooting_is_explicit_and_validated() {
        let t = tree("(((A,B),(C,D)),E);");
        let m = build_matrix(&t, &[("1:0", 1, 0, &["A", "B"])], &[]);
        let (d, _) = prepare(&m, &t, &["E".into()], 1).unwrap();
        assert_eq!(
            d.tree.children(d.tree.root()).len(),
            2,
            "a new binary root replaced the trifurcation"
        );
        assert_eq!(d.tree.n_tips(), 5, "no taxon was lost by re-rooting");
        assert!(!d.tree.looks_unrooted());
        assert!(prepare(&m, &t, &["nope".into()], 1).is_err());
        assert!(
            prepare(&m, &t, &["A".into(), "C".into()], 1).is_err(),
            "non-monophyletic outgroup"
        );
    }

    #[test]
    fn disjoint_taxa_sets_are_reported() {
        let m = build_matrix(&tree("((X,Y),Z);"), &[("1:0", 1, 0, &["X"])], &[]);
        let e = prepare(&m, &tree("((A,B),C);"), &[], 1).unwrap_err();
        assert!(
            e.to_string().contains("no matrix taxon matches a tree tip"),
            "{e}"
        );
    }

    #[test]
    fn unannotated_matrix_taxa_are_not_scored_as_absent() {
        let t = tree("((A,B),C);");
        let m = build_matrix(&t, &[("1:0", 1, 0, &["A"])], &["A"]);
        let (d, rep) = prepare(&m, &t, &[], 1).unwrap();
        assert_eq!(d.characters[0].states[d.tip_index_of("B")], MISSING);
        assert_eq!(rep.used, vec!["A".to_string()]);
    }

    #[test]
    fn counts_are_tree_free() {
        let t = tree("((A,B),(C,D));");
        let m = build_matrix(
            &t,
            &[
                ("1:0", 1, 0, &["A", "B"]),
                ("2:0", 2, 0, &["A"]),
                ("3:1", 3, 1, &["A", "B", "C", "D"]),
            ],
            &[],
        );
        let c = score_counts(&m);
        assert_eq!(c.len(), 3);
        assert_eq!((c[0].n_present, c[0].n_absent), (2, 2));
        assert!((c[0].frequency - 0.5).abs() < 1e-12);
        assert!((c[2].frequency - 1.0).abs() < 1e-12);
        assert_eq!(c[1].taxa, vec!["A"]);
    }

    #[test]
    fn method_names_parse() {
        assert!(matches!(Method::parse("Dollo"), Some(Method::Dollo)));
        assert!(matches!(Method::parse("fitch"), Some(Method::MaxParsimony)));
        assert!(matches!(Method::parse("count"), Some(Method::Counts)));
        assert_eq!(Method::parse("nonsense"), None);
        assert_eq!(Method::Dollo.to_string(), "dollo");
    }

    // ---- proposed: weighted Dollo vs brute force, random trees, ambiguity ----
    fn brute_force_weighted(tree: &Tree, tips: &[NodeId], obs: &[u8], gw: u64, lw: u64) -> u64 {
        let n = tree.len();
        let internal: Vec<NodeId> = (0..n).filter(|&v| !tree.children(v).is_empty()).collect();
        let mut tip_obs = vec![MISSING; n];
        for (i, &t) in tips.iter().enumerate() {
            tip_obs[t] = obs[i];
        }
        // 0 = not yet gained, 1 = present, 2 = lost
        let edge = |a: u8, b: u8| -> Option<u64> {
            match (a, b) {
                (0, 0) | (1, 1) | (2, 2) => Some(0),
                (0, 1) => Some(gw),
                (1, 2) => Some(lw),
                _ => None,
            }
        };
        let allowed = |v: NodeId, s: u8| match tip_obs[v] {
            PRESENT => s == 1,
            ABSENT => s != 1,
            _ => true,
        };
        let mut best = u64::MAX;
        'outer: for code0 in 0..3u64.pow(internal.len() as u32) {
            let mut code = code0;
            let mut st = vec![0u8; n];
            for &v in &internal {
                st[v] = (code % 3) as u8;
                code /= 3;
            }
            let mut cost = match st[tree.root()] {
                0 => 0,
                1 => gw,
                _ => continue,
            };
            for &v in &internal {
                for &c in tree.children(v) {
                    let x = if tree.children(c).is_empty() {
                        (0..3u8)
                            .filter(|&s| allowed(c, s))
                            .filter_map(|s| edge(st[v], s))
                            .min()
                    } else {
                        edge(st[v], st[c])
                    };
                    match x {
                        Some(x) => cost += x,
                        None => continue 'outer,
                    }
                }
            }
            best = best.min(cost);
        }
        best
    }

    #[test]
    fn weighted_dollo_matches_brute_force_with_missing_tips() {
        for nw in SHAPES {
            let t = tree(nw);
            let tips = t.tips();
            for (gw, lw) in [(1u32, 1u32), (2, 1), (3, 1), (1, 2), (5, 2)] {
                for obs in patterns(tips.len()) {
                    let (g, l) = weighted_one(&t, &tips, &obs, gw, lw);
                    let got = u64::from(g) * u64::from(gw) + u64::from(l) * u64::from(lw);
                    let want = brute_force_weighted(&t, &tips, &obs, gw.into(), lw.into());
                    assert_eq!(got, want, "{nw} obs {obs:?} gw {gw} lw {lw}: {g}g {l}l");
                }
            }
        }
    }

    struct Rng(u64);
    impl Rng {
        fn below(&mut self, n: u64) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 % n
        }
    }

    fn random_newick(rng: &mut Rng, n_tips: usize) -> String {
        let mut parts: Vec<String> = (0..n_tips).map(|i| format!("T{i}")).collect();
        while parts.len() > 1 {
            let k = if parts.len() >= 3 && rng.below(4) == 0 {
                3
            } else {
                2
            };
            let mut group = Vec::new();
            for _ in 0..k {
                let i = rng.below(parts.len() as u64) as usize;
                group.push(parts.swap_remove(i));
            }
            parts.push(format!("({})", group.join(",")));
        }
        format!("{};", parts[0])
    }

    #[test]
    fn dollo_and_fitch_match_brute_force_on_random_trees() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for _ in 0..40 {
            let n = 3 + rng.below(5) as usize; // 3..=7 tips
            let nw = random_newick(&mut rng, n);
            let t = tree(&nw);
            let tips = t.tips();
            for _ in 0..80 {
                let obs: Vec<u8> = (0..n)
                    .map(|_| [PRESENT, ABSENT, ABSENT, MISSING][rng.below(4) as usize])
                    .collect();
                let d = dollo_one(&t, &tips, &obs, true);
                assert_eq!(
                    d.losses,
                    brute_force_losses(&t, &tips, &obs),
                    "{nw} {obs:?}"
                );
                // The traceback must realise exactly the reported losses.
                if let Some(r) = &d.recon {
                    let realised = t
                        .edges()
                        .filter(|&(p, c)| r[p] == D_PRESENT && r[c] == D_LOST)
                        .count();
                    assert_eq!(realised as u32, d.losses, "{nw} {obs:?}");
                    assert_eq!(d.loss_ids.len() as u32, d.losses);
                }
                assert_eq!(
                    fitch_one(&t, &tips, &obs, false).0,
                    brute_force_mp(&t, &tips, &obs),
                    "{nw} {obs:?}"
                );
            }
        }
    }

    #[test]
    fn a_three_to_one_split_on_four_taxa_has_a_unique_reconstruction() {
        // Only most-parsimonious labelling: every internal node present, one change on D.
        let t = tree("((A,B),(C,D));");
        let tips = t.tips();
        let (steps, ambiguous, _) =
            fitch_one(&t, &tips, &[PRESENT, PRESENT, PRESENT, ABSENT], false);
        assert_eq!(steps, 1);
        assert!(!ambiguous, "the reconstruction is unique");
    }

    /// Number of optimal labellings of the informed internal nodes, by exhaustive search.
    fn brute_force_optimal_labellings(tree: &Tree, tips: &[NodeId], obs: &[u8]) -> u32 {
        let n = tree.len();
        let mut tip_obs = vec![MISSING; n];
        for (i, &t) in tips.iter().enumerate() {
            tip_obs[t] = obs[i];
        }
        let mut informed = vec![false; n];
        for &v in tree.postorder() {
            informed[v] = if tree.children(v).is_empty() {
                tip_obs[v] != MISSING
            } else {
                tree.children(v).iter().any(|&c| informed[c])
            };
        }
        let internal: Vec<NodeId> = (0..n)
            .filter(|&v| !tree.children(v).is_empty() && informed[v])
            .collect();
        let (mut best, mut count) = (u32::MAX, 0u32);
        for mask in 0u32..(1 << internal.len()) {
            let mut st = vec![0u8; n];
            for (k, &v) in internal.iter().enumerate() {
                st[v] = (mask >> k & 1) as u8;
            }
            let mut cost = 0u32;
            for &v in &internal {
                for &c in tree.children(v) {
                    if !informed[c] {
                        continue;
                    }
                    let cs = if tree.children(c).is_empty() {
                        tip_obs[c]
                    } else {
                        st[c]
                    };
                    cost += u32::from(st[v] != cs);
                }
            }
            if cost < best {
                best = cost;
                count = 1;
            } else if cost == best {
                count += 1;
            }
        }
        count
    }

    #[test]
    fn the_ambiguous_flag_means_several_optimal_reconstructions() {
        for nw in SHAPES {
            let t = tree(nw);
            let tips = t.tips();
            for obs in patterns(tips.len()) {
                if obs.iter().all(|&o| o == MISSING) {
                    continue;
                }
                let want = brute_force_optimal_labellings(&t, &tips, &obs) > 1;
                let (_, got, _) = fitch_one(&t, &tips, &obs, false);
                assert_eq!(got, want, "{nw} obs {obs:?}");
            }
        }
    }
}
