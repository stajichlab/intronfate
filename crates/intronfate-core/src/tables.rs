//! Tabular output.
//!
//! Everything is plain TSV with a header row: positionable with `awk`/`datamash`, loadable into R
//! and pandas without a reader library, and diffable between runs. Every table repeats the full
//! coordinate stack (`cds_offset`, `phase`, `codon`, `residue_before`, `residue_after`, `column`) so
//! a slice of one table is never ambiguous about which space a number lives in.

use std::fmt::Write as _;
use std::path::Path;

use crate::error::Result;
use crate::matrix::{PairStat, SiteMatrix};
use crate::project::{AlignmentProjection, Projection};
use crate::report::Diagnostics;
use crate::score::{CountsSite, DolloResult, MpResult, ScoringData};
use crate::seqio::ProteinSet;

/// An in-memory TSV table.
#[derive(Clone, Debug, Default)]
pub struct Table {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
    /// Optional leading `# key=value` line: format metadata that TSV has nowhere else for.
    pub comment: Option<String>,
}

pub trait Cells {
    fn cells(&self, out: &mut Vec<String>);
}

impl<F: Fn(&mut Vec<String>)> Cells for F {
    fn cells(&self, out: &mut Vec<String>) {
        self(out)
    }
}

impl Table {
    pub fn new(header: &[&str]) -> Self {
        Table {
            header: header.iter().map(|s| s.to_string()).collect(),
            rows: Vec::new(),
            comment: None,
        }
    }

    pub fn with_capacity(header: &[&str], rows: usize) -> Self {
        let mut t = Table::new(header);
        t.rows = Vec::with_capacity(rows);
        t
    }

    pub fn push(&mut self, cells: Vec<String>) {
        debug_assert!(
            cells.len() == self.header.len(),
            "row has {} cells but the header has {}: {cells:?}",
            cells.len(),
            self.header.len()
        );
        self.rows.push(cells);
    }

    /// Build a row through a closure, keeping column order next to the header definition.
    pub fn push_row(&mut self, f: impl Fn(&mut Vec<String>)) {
        let mut cells = Vec::with_capacity(self.header.len());
        f(&mut cells);
        self.push(cells);
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    pub fn to_tsv(&self) -> String {
        let mut s = String::new();
        if let Some(c) = &self.comment {
            s.push_str("# ");
            s.push_str(c);
            s.push('\n');
        }
        s.push_str(&self.header.join("\t"));
        s.push('\n');
        for r in &self.rows {
            s.push_str(&r.join("\t"));
            s.push('\n');
        }
        s
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<usize> {
        let path = path.as_ref();
        crate::util::write_atomic(path, self.to_tsv().as_bytes())?;
        Ok(self.rows.len())
    }
}

fn cell(v: impl std::fmt::Display) -> String {
    v.to_string()
}

fn opt(o: Option<impl std::fmt::Display>) -> String {
    o.map(|v| v.to_string()).unwrap_or_default()
}

fn list(items: &[impl std::fmt::Display]) -> String {
    items
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// Write `introns.tsv` straight to disk, one row at a time. Byte-identical to
/// `introns_table(proj).write(path)`, without the one-`String`-per-cell table (a 5·10⁵-intron genome made
/// ~9·10⁶ small strings and held the whole file twice).
pub fn write_introns_tsv(proj: &Projection, path: impl AsRef<Path>) -> Result<usize> {
    crate::util::write_atomic_with(path.as_ref(), |w| {
        writeln!(w, "{}", INTRON_HEADER.join("\t"))?;
        for r in &proj.rows {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                proj.protein_id(r),
                proj.gene_id(r).unwrap_or(""),
                proj.transcript_id(r),
                proj.seqid(r),
                proj.strand(r),
                r.intron_rank,
                r.coding_rank,
                r.genomic_start,
                r.genomic_end,
                r.intron_len,
                r.site.cds_offset,
                r.site.phase,
                r.site.codon(),
                r.site.residue_before(),
                r.site.residue_after(),
                r.site.cds_offset,
                r.protein_len,
                if r.in_range { "yes" } else { "no" },
            )?;
        }
        Ok(())
    })?;
    Ok(proj.rows.len())
}

const INTRON_HEADER: [&str; 18] = [
    "protein_id",
    "gene_id",
    "transcript_id",
    "seqid",
    "strand",
    "intron_rank",
    "coding_rank",
    "genomic_start",
    "genomic_end",
    "intron_len",
    "cds_offset",
    "phase",
    "codon",
    "residue_before",
    "residue_after",
    "site_key",
    "protein_len",
    "in_range",
];

/// Every coding intron in annotation space: genome coordinates plus all protein-space views.
pub fn introns_table(proj: &Projection) -> Table {
    let mut t = Table::new(&INTRON_HEADER);
    t.rows.reserve(proj.rows.len());
    for r in &proj.rows {
        t.push_row(|c| {
            c.push(proj.protein_id(r).to_string());
            c.push(opt(proj.gene_id(r)));
            c.push(proj.transcript_id(r).to_string());
            c.push(proj.seqid(r).to_string());
            c.push(cell(proj.strand(r)));
            c.push(cell(r.intron_rank));
            c.push(cell(r.coding_rank));
            c.push(cell(r.genomic_start));
            c.push(cell(r.genomic_end));
            c.push(cell(r.intron_len));
            c.push(cell(r.site.cds_offset));
            c.push(cell(r.site.phase));
            c.push(cell(r.site.codon()));
            c.push(cell(r.site.residue_before()));
            c.push(cell(r.site.residue_after()));
            c.push(r.site.site_key());
            c.push(cell(r.protein_len));
            c.push(cell(if r.in_range { "yes" } else { "no" }));
        });
    }
    t
}

/// Projected sites for one alignment: one row per taxon per intron position.
pub fn sites_table(ap: &AlignmentProjection, proteins: &ProteinSet) -> Table {
    let mut t = Table::new(&[
        "alignment",
        "taxon",
        "protein_id",
        "column",
        "phase",
        "cds_offset",
        "codon",
        "residue_before",
        "residue_after",
        "site_label",
        "protein_len",
        "residues_aligned",
    ]);
    t.rows.reserve(ap.sites.len());
    for taxon in &ap.taxa {
        let Some(pi) = taxon.protein_idx else {
            continue;
        };
        if taxon.sites.is_empty() {
            continue;
        }
        let pid = proteins.id(pi).unwrap_or("").to_string();
        let plen = proteins.get(pi).map_or(0, |s| s.len);
        for s in &taxon.sites {
            t.push_row(|c| {
                c.push(ap.alignment.clone());
                c.push(taxon.taxon_name.clone());
                c.push(pid.clone());
                c.push(cell(s.col));
                c.push(cell(s.site.phase));
                c.push(cell(s.site.cds_offset));
                c.push(cell(s.site.codon()));
                c.push(cell(s.site.residue_before()));
                c.push(cell(s.site.residue_after()));
                c.push(format!("{}:{}", s.col, s.site.phase));
                c.push(cell(plen));
                c.push(cell(taxon.residues_aligned));
            });
        }
    }
    t
}

/// Per-taxon participation and quality: what was projected, what was dropped and why.
pub fn taxa_table(ap: &AlignmentProjection, proteins: &ProteinSet, m: &SiteMatrix) -> Table {
    let mut t = Table::new(&[
        "alignment",
        "taxon",
        "protein_id",
        "protein_len",
        "residues_aligned",
        "n_sites",
        "n_intron_records",
        "tail_sites",
        "clamped_sites",
        "annotated",
        "reason_if_not",
    ]);
    for taxon in &ap.taxa {
        let pid = taxon
            .protein_idx
            .and_then(|i| proteins.id(i))
            .unwrap_or("")
            .to_string();
        let plen = taxon
            .protein_idx
            .and_then(|i| proteins.residues(i))
            .map_or(0, |r| r.len());
        let n_sites = m.taxon_sites.get(taxon.taxon).map_or(0, |v| v.len());
        t.push_row(|c| {
            c.push(ap.alignment.clone());
            c.push(taxon.taxon_name.clone());
            c.push(pid.clone());
            c.push(cell(plen));
            c.push(cell(taxon.residues_aligned));
            c.push(cell(n_sites));
            c.push(cell(taxon.report.n_sites()));
            c.push(cell(taxon.report.tail_sites));
            c.push(cell(taxon.report.clamped));
            c.push(cell(if taxon.has_annotation { "yes" } else { "no" }));
            c.push(opt(taxon.unannotated_reason));
        });
    }
    t
}

/// Write `matrix.tsv` straight to disk, row by row. Byte-identical to
/// `matrix_table(m).write(path)`, but peak memory is one row instead of one `String` per cell (a
/// 1e5-site × 1e3-taxon matrix is 1e8 cells).
pub fn write_matrix_tsv(m: &SiteMatrix, path: impl AsRef<Path>) -> Result<usize> {
    crate::util::write_atomic_with(path.as_ref(), |w| {
        writeln!(w, "# {}", matrix_comment(m))?;
        write!(w, "site\tcolumn\tphase")?;
        if m.is_combined() {
            write!(w, "\tfamily\tfamily_column")?;
        }
        write!(w, "\tn_taxa\tn_observable\tpattern")?;
        for name in &m.taxa {
            write!(w, "\t{name}")?;
        }
        writeln!(w)?;
        let mut cells = String::with_capacity(m.taxa.len() * 2);
        for (i, site) in m.sites.iter().enumerate() {
            cells.clear();
            for t in 0..m.taxa.len() {
                cells.push('\t');
                cells.push(site.cell(t).symbol());
            }
            write!(w, "{}\t{}\t{}\t", m.site_label(i), site.col, site.phase)?;
            if let Some(&f) = m.site_family.get(i) {
                write!(w, "{}\t{}\t", m.families[f as usize].name, m.site_column(i))?;
            }
            write!(w, "{}\t{}\t", site.n_taxa(), site.n_observable())?;
            for t in 0..m.taxa.len() {
                write!(w, "{}", site.cell(t).symbol())?;
            }
            w.write_all(cells.as_bytes())?;
            writeln!(w)?;
        }
        Ok(())
    })?;
    Ok(m.sites.len())
}

/// The binary `site × taxon` matrix, with the pattern string first for quick reading.
pub fn matrix_table(m: &SiteMatrix) -> Table {
    let mut header = vec!["site", "column", "phase"];
    if m.is_combined() {
        header.extend(["family", "family_column"]);
    }
    header.extend(["n_taxa", "n_observable", "pattern"]);
    header.extend(m.taxa.iter().map(|s| s.as_str()));
    let mut t = Table::new(&header);
    for (i, site) in m.sites.iter().enumerate() {
        t.push_row(|c| {
            c.push(m.site_label(i));
            c.push(cell(site.col));
            c.push(cell(site.phase));
            if let Some(&f) = m.site_family.get(i) {
                c.push(m.families[f as usize].name.clone());
                c.push(cell(m.site_column(i)));
            }
            c.push(cell(site.n_taxa()));
            c.push(cell(site.n_observable()));
            c.push(m.pattern(i));
            for t_i in 0..m.taxa.len() {
                c.push(site.cell(t_i).symbol().to_string());
            }
        });
    }
    t.comment = Some(matrix_comment(m));
    t
}

/// The `# key=value ...` metadata line of `matrix.tsv`.
fn matrix_comment(m: &SiteMatrix) -> String {
    let unannotated: Vec<&str> = m
        .taxa
        .iter()
        .zip(&m.annotated)
        .filter(|(_, a)| !*a)
        .map(|(n, _)| n.as_str())
        .collect();
    format!(
        "format=2 alignment={} alen={} phase_mode={} flank={} masked_present={} taxa={} sites={} unannotated={} residues={}{}",
        m.alignment.replace(' ', "_"),
        m.alen,
        if matches!(m.phase_mode, crate::project::PhaseMode::Ignore) {
            "ignore"
        } else {
            "strict"
        },
        m.flank,
        m.masked_present,
        m.taxa.len(),
        m.sites.len(),
        if unannotated.is_empty() {
            "-".to_string()
        } else {
            unannotated
                .iter()
                .map(|n| crate::util::encode_list_item(n))
                .collect::<Vec<_>>()
                .join(",")
        },
        if m.aligned_residues.is_empty() {
            "-".to_string()
        } else {
            m.aligned_residues
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",")
        },
        if m.is_combined() {
            format!(" families={} family_pad={}", m.families.len(), m.family_pad)
        } else {
            String::new()
        },
    )
}

/// Long form of the matrix: one row per (site, taxon) *presence* only.
pub fn presence_table(m: &SiteMatrix) -> Table {
    let mut t = Table::new(&["site", "column", "phase", "taxon"]);
    for (i, site) in m.sites.iter().enumerate() {
        for tx in site.taxa() {
            if let Some(name) = m.taxa.get(tx) {
                t.push(vec![
                    m.site_label(i),
                    cell(m.site_column(i)),
                    cell(site.phase),
                    name.clone(),
                ]);
            }
        }
    }
    t
}

pub fn pairwise_table(stats: &[PairStat]) -> Table {
    let mut t = Table::new(&[
        "taxon_a",
        "taxon_b",
        "n_a",
        "n_b",
        "shared",
        "jaccard",
        "fraction_of_smaller",
        "both_annotated",
        "n_both_observable",
    ]);
    for p in stats {
        t.push_row(|c| {
            c.push(p.taxon_a.clone());
            c.push(p.taxon_b.clone());
            c.push(cell(p.n_a));
            c.push(cell(p.n_b));
            c.push(cell(p.shared));
            c.push(format!("{:.6}", p.jaccard));
            c.push(format!("{:.6}", p.fraction_of_smaller));
            c.push(cell(if p.both_annotated { "yes" } else { "no" }));
            c.push(cell(p.n_both_observable));
        });
    }
    t
}

/// Dense symmetric sharing matrix (PHYLIP-ish square matrix with named rows and columns).
pub fn sharing_matrix_table(m: &SiteMatrix, tolerance: u32) -> Table {
    sharing_matrix_table_from(m, &m.pairwise(tolerance, true))
}

/// [`sharing_matrix_table`] from pairs that were already computed (at the same tolerance and phase rule).
pub fn sharing_matrix_table_from(m: &SiteMatrix, pairs: &[PairStat]) -> Table {
    let dense = m.matrix_from_pairs(pairs);
    let mut header = vec!["taxon".to_string(), "n_sites".to_string()];
    header.extend(m.taxa.iter().cloned());
    let mut t = Table {
        header,
        rows: Vec::with_capacity(m.taxa.len()),
        comment: None,
    };
    let counts = m.per_taxon_counts();
    for (i, name) in m.taxa.iter().enumerate() {
        let mut row = vec![name.clone(), cell(counts[i])];
        for j in 0..m.taxa.len() {
            row.push(cell(dense[i * m.taxa.len() + j]));
        }
        t.rows.push(row);
    }
    t
}

/// How many sites are shared by exactly k taxa.
pub fn spectrum_table(m: &SiteMatrix) -> Table {
    let h = m.conservation_spectrum();
    let mut t = Table::new(&["n_taxa_with_site", "n_sites", "share_of_sites"]);
    let total: u64 = h.iter().sum();
    for (k, &n) in h.iter().enumerate() {
        if n == 0 {
            continue;
        }
        t.push(vec![
            cell(k),
            cell(n),
            if total == 0 {
                "0".into()
            } else {
                format!("{:.6}", n as f64 / total as f64)
            },
        ]);
    }
    t
}

// Node ids index the tree, the names and the results together.
#[allow(clippy::needless_range_loop)]
/// Per node: fitted transition probabilities and posterior expectations (with bootstrap intervals when given).
pub fn likelihood_branch_table(
    data: &ScoringData,
    r: &crate::likelihood::FitResult,
    boot: Option<&crate::likelihood::Bootstrap>,
) -> Table {
    let mut head = vec![
        "branch",
        "node",
        "n_tips",
        "gain_prob",
        "loss_prob",
        "expected_gains",
        "expected_losses",
        "expected_present",
        "has_data",
        "fixed",
    ];
    if boot.is_some() {
        head.extend([
            "gains_lo",
            "gains_hi",
            "losses_lo",
            "losses_hi",
            "present_lo",
            "present_hi",
        ]);
    }
    head.push("tips");
    let mut t = Table::new(&head);
    let names = data.tree.node_names();
    let f = |x: f64| format!("{x:.2}");
    let pr = |x: f64| format!("{x:.6}");
    for v in 0..data.tree.len() {
        let n = &r.nodes[v];
        let tips: Vec<String> = data
            .tree
            .leaves_under(v)
            .iter()
            .map(|&l| data.tree.label(l).unwrap_or("?").to_string())
            .collect();
        t.push_row(|c| {
            c.push(names[v].clone());
            c.push(cell(v));
            c.push(cell(tips.len()));
            c.push(pr(n.gain_prob));
            c.push(pr(n.loss_prob));
            c.push(f(n.expected_gains));
            c.push(f(n.expected_losses));
            c.push(f(n.expected_present));
            c.push(cell(if n.has_data { "yes" } else { "no" }));
            c.push(cell(n.fixed.as_str()));
            if let Some(b) = boot {
                for (lo, hi) in [b.gains[v], b.losses[v], b.present[v]] {
                    c.push(f(lo));
                    c.push(f(hi));
                }
            }
            c.push(tips.join(","));
        });
    }
    t
}

/// Per tip: the error rates used and the expected misses and false introns.
pub fn likelihood_tip_table(
    data: &ScoringData,
    err: &crate::likelihood::ErrorRates,
    r: &crate::likelihood::FitResult,
) -> Table {
    let mut t = Table::new(&[
        "taxon",
        "sensitivity",
        "false_rate",
        "observed_present",
        "observed_absent",
        "expected_missed",
        "expected_false",
    ]);
    for (i, name) in data.tip_names.iter().enumerate() {
        let count = |s: u8| data.characters.iter().filter(|c| c.states[i] == s).count();
        t.push_row(|c| {
            c.push(name.clone());
            c.push(format!("{:.4}", err.sensitivity[i]));
            c.push(format!("{:.6}", err.false_rate[i]));
            c.push(cell(count(crate::score::PRESENT)));
            c.push(cell(count(crate::score::ABSENT)));
            c.push(format!("{:.2}", r.tips[i].expected_missed));
            c.push(format!("{:.2}", r.tips[i].expected_false));
        });
    }
    t
}

/// Conditional log-likelihood per position (for site-wise tree comparisons, e.g. CONSEL/AU).
pub fn likelihood_site_table(data: &ScoringData, r: &crate::likelihood::FitResult) -> Table {
    let mut t = Table::with_capacity(&["site", "log_lik"], data.characters.len());
    for (ch, l) in data.characters.iter().zip(&r.site_log_lik) {
        t.push(vec![ch.label.clone(), format!("{l:.6}")]);
    }
    t
}

/// One weighting of a weighted-Dollo sweep: `(gain_weight, loss_weight, [(site, gains, losses)])`.
pub type WeightedDolloRun = (u32, u32, Vec<(String, u32, u32)>);

/// Per-site weighted-Dollo results, long form: one row per site and weighting.
pub fn dollo_weighted_table(runs: &[WeightedDolloRun]) -> Table {
    let n = runs.iter().map(|r| r.2.len()).sum();
    let mut t = Table::with_capacity(
        &[
            "site",
            "gain_weight",
            "loss_weight",
            "gains",
            "losses",
            "cost",
        ],
        n,
    );
    for (gw, lw, sites) in runs {
        for (label, g, l) in sites {
            t.push_row(|c| {
                c.push(label.clone());
                c.push(cell(gw));
                c.push(cell(lw));
                c.push(cell(g));
                c.push(cell(l));
                c.push(cell(
                    u64::from(*gw) * u64::from(*g) + u64::from(*lw) * u64::from(*l),
                ));
            });
        }
    }
    t
}

/// Totals per weighting of a weighted-Dollo sweep.
pub fn dollo_weighted_summary_table(runs: &[WeightedDolloRun]) -> Table {
    let mut t = Table::new(&[
        "gain_weight",
        "loss_weight",
        "sites",
        "gains",
        "losses",
        "sites_one_gain",
        "sites_several_gains",
        "sites_with_losses",
    ]);
    for (gw, lw, sites) in runs {
        let gains: u64 = sites.iter().map(|s| u64::from(s.1)).sum();
        let losses: u64 = sites.iter().map(|s| u64::from(s.2)).sum();
        t.push_row(|c| {
            c.push(cell(gw));
            c.push(cell(lw));
            c.push(cell(sites.len()));
            c.push(cell(gains));
            c.push(cell(losses));
            c.push(cell(sites.iter().filter(|s| s.1 == 1).count()));
            c.push(cell(sites.iter().filter(|s| s.1 > 1).count()));
            c.push(cell(sites.iter().filter(|s| s.2 > 0).count()));
        });
    }
    t
}

pub fn dollo_table(r: &DolloResult) -> Table {
    let mut t = Table::new(&[
        "site",
        "column",
        "phase",
        "n_present",
        "n_missing",
        "gains",
        "losses",
        "monophyletic",
        "dollo_compatible",
        "gained_at",
        "lost_in",
    ]);
    for s in &r.sites {
        t.push_row(|c| {
            c.push(s.label.clone());
            c.push(cell(s.col));
            c.push(cell(s.phase));
            c.push(cell(s.n_present));
            c.push(cell(s.n_missing));
            c.push(cell(s.gains));
            c.push(cell(s.losses));
            c.push(cell(if s.monophyletic { "yes" } else { "no" }));
            c.push(cell(if s.compatible { "yes" } else { "no" }));
            c.push(opt(s.gain_node.as_deref()));
            // `;` between branches: a branch is itself named by a comma-separated tip list.
            c.push(s.loss_edges.join(";"));
        });
    }
    t
}

/// Dollo gains and losses per branch of the tree, summed over sites. Branch `losses` is the count of
/// sites lost on the branch leading *to* that node; `gains` the origins placed at it.
pub fn dollo_branch_table(r: &DolloResult) -> Table {
    let mut t = Table::new(&["branch", "node", "n_tips", "gains", "losses", "tips"]);
    for b in &r.branches {
        t.push_row(|c| {
            c.push(b.name.clone());
            c.push(cell(b.node));
            c.push(cell(b.tips.len()));
            c.push(cell(b.gains));
            c.push(cell(b.losses));
            c.push(list(&b.tips));
        });
    }
    t
}

pub fn mp_table(r: &MpResult) -> Table {
    let mut t = Table::new(&[
        "site",
        "column",
        "phase",
        "n_present",
        "n_missing",
        "steps",
        "informative",
        "ambiguous",
    ]);
    for s in &r.sites {
        t.push_row(|c| {
            c.push(s.label.clone());
            c.push(cell(s.col));
            c.push(cell(s.phase));
            c.push(cell(s.n_present));
            c.push(cell(s.n_missing));
            c.push(cell(s.steps));
            c.push(cell(if s.informative { "yes" } else { "no" }));
            c.push(cell(if s.ambiguous { "yes" } else { "no" }));
        });
    }
    t
}

pub fn counts_table(r: &[CountsSite]) -> Table {
    let mut t = Table::new(&[
        "site",
        "column",
        "phase",
        "n_present",
        "n_absent",
        "n_missing",
        "frequency",
        "taxa",
    ]);
    for s in r {
        t.push_row(|c| {
            c.push(s.label.clone());
            c.push(cell(s.col));
            c.push(cell(s.phase));
            c.push(cell(s.n_present));
            c.push(cell(s.n_absent));
            c.push(cell(s.n_missing));
            c.push(format!("{:.6}", s.frequency));
            c.push(list(&s.taxa));
        });
    }
    t
}

/// Per-node ancestral states, one row per (method, site, node).
pub fn reconstruction_table(data: &ScoringData, sites: &[(&str, &str, Option<Vec<u8>>)]) -> Table {
    let mut t = Table::new(&["method", "site", "node", "is_tip", "state"]);
    let node_names: Vec<String> = (0..data.tree.len())
        .map(|i| {
            data.tree.label(i).map(str::to_string).unwrap_or_else(|| {
                if data.tree.node(i).is_leaf() {
                    format!("tip{i}")
                } else {
                    format!("node{}", data.tree.edge_name(i).replace(',', "-"))
                }
            })
        })
        .collect();
    for (method, label, recon) in sites {
        let Some(recon) = recon else { continue };
        for (i, &st) in recon.iter().enumerate() {
            let name = if *method == "dollo" {
                match st {
                    1 => "present",
                    2 => "lost",
                    _ => "absent",
                }
            } else {
                match st {
                    0 => "absent",
                    1 => "present",
                    _ => "missing",
                }
            };
            t.push(vec![
                method.to_string(),
                label.to_string(),
                node_names[i].clone(),
                cell(data.tree.node(i).is_leaf()),
                name.to_string(),
            ]);
        }
    }
    t
}

/// Unresolved annotation records, so nothing silently disappears.
pub fn unmatched_table(proj: &Projection) -> Table {
    let mut t = Table::new(&[
        "transcript_id",
        "gene_id",
        "attempted_protein_id",
        "reason",
        "source",
    ]);
    for u in &proj.unmatched {
        t.push(vec![
            u.transcript_id.clone(),
            opt(u.gene_id.as_deref()),
            u.attempted_id.clone(),
            u.reason.to_string(),
            u.source.clone(),
        ]);
    }
    for p in &proj.proteins_without_annotation {
        t.push(vec![
            String::new(),
            String::new(),
            p.clone(),
            "protein_without_annotation".into(),
            String::new(),
        ]);
    }
    t
}

/// Flattened run diagnostics: `kind` can be grepped, `subject` traced back to an input row.
pub fn diagnostics_table(d: &Diagnostics) -> Table {
    let mut t = Table::new(&["level", "kind", "subject", "message"]);
    for x in d.retained() {
        t.push(vec![
            x.level.to_string(),
            x.kind.to_string(),
            x.subject.clone(),
            x.message.clone(),
        ]);
    }
    t
}

pub fn diagnostics_tally_table(d: &Diagnostics) -> Table {
    let mut t = Table::new(&["kind", "count"]);
    for (k, n) in d.tally() {
        t.push(vec![k.to_string(), cell(n)]);
    }
    t
}

/// Key/value run summary.
pub fn summary_table(entries: &[(&str, String)]) -> Table {
    let mut t = Table::new(&["key", "value"]);
    for (k, v) in entries {
        t.push(vec![k.to_string(), v.clone()]);
    }
    t
}

/// Newick with each branch annotated with the number of Dollo gains placed on it, for FigTree-style
/// tools that read comments.
pub fn tree_with_gain_counts(data: &ScoringData, dollo: &DolloResult) -> String {
    // Gains are indexed by node id: matching on printed names would put a gain on the true root
    // whenever an internal node happens to be labelled `root`.
    let mut per_edge = vec![0u32; data.tree.len()];
    for b in &dollo.branches {
        per_edge[b.node] = b.gains as u32;
    }
    let mut out = String::new();
    write_annotated(&data.tree, data.tree.root(), &per_edge, &mut out);
    out.push(';');
    out
}

fn write_annotated(tree: &crate::tree::Tree, n: usize, per_edge: &[u32], out: &mut String) {
    let node = tree.node(n);
    if !node.children.is_empty() {
        out.push('(');
        for (i, &c) in node.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            write_annotated(tree, c, per_edge, out);
        }
        out.push(')');
    }
    if let Some(l) = &node.label {
        let _ = write!(out, "{l}");
    }
    if per_edge[n] > 0 {
        let _ = write!(out, "[gains={}]", per_edge[n]);
    }
    if let Some(bl) = node.branch_length {
        let _ = write!(out, ":{bl}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{AnnotatedIntron, ProteinSite, Strand, TranscriptAnnotation};
    use crate::project::{join_annotations, project_into_alignment, JoinOptions, ProjectOptions};
    use crate::score::{prepare, score_counts, score_dollo, score_mp};
    use crate::seqio::ProteinSet;
    use crate::tree::Tree;
    use std::io::Cursor;

    fn fixture() -> (SiteMatrix, AlignmentProjection, Projection, ProteinSet) {
        let mut ps = ProteinSet::new();
        ps.read_one(
            Cursor::new(b">A\nMKVLGAWEQP\n>B\nAWEQP\n>C\nMKVLGAWEQP\n"),
            std::path::Path::new("m"),
        )
        .unwrap();
        let mk = |p: &str, off: u32| TranscriptAnnotation {
            transcript_id: format!("{p}.t"),
            gene_id: Some(format!("{p}.g")),
            protein_id: Some(p.into()),
            seqid: "chr1".into(),
            strand: Strand::Plus,
            cds_len: u64::from(off) + 3,
            cds_segments: 2,
            introns: vec![AnnotatedIntron {
                start: 10,
                end: 20,
                rank: 0,
                cds_offset: Some(off),
            }],
            frame_mismatches: 0,
            source: "test".into(),
        };
        let proj = join_annotations(
            &[mk("A", 18), mk("B", 3), mk("C", 18)],
            &ps,
            &JoinOptions {
                length_slack: u64::MAX,
                ..Default::default()
            },
        )
        .unwrap();
        let a = crate::msa::Alignment::from_fasta_str(
            ">A\nMKVLGAWEQP\n>B\n-----AWEQP\n>C\nMKVLGAWEQP\n",
            "fam",
        )
        .unwrap();
        let ap = project_into_alignment(&a, &proj, &ps, &ProjectOptions::default());
        let mut m = SiteMatrix::build(&ap);
        m.alen = a.len;
        (m, ap, proj, ps)
    }

    #[test]
    fn introns_table_carries_every_coordinate_view() {
        let (_, _, proj, _) = fixture();
        let t = introns_table(&proj);
        assert_eq!(t.header.len(), 18);
        assert_eq!(t.len(), 3);
        let s = t.to_tsv();
        let row = s.lines().nth(1).unwrap();
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols[0], "A", "sorted by protein id");
        assert_eq!(cols[10], "18", "cds_offset");
        assert_eq!(cols[11], "0", "phase");
        assert_eq!(cols[12], "6", "codon");
        assert_eq!(cols[13], "6", "residue_before");
        assert_eq!(cols[14], "7", "residue_after");
        assert_eq!(cols[16], "10", "protein_len");
        assert_eq!(cols[17], "yes");
    }

    #[test]
    fn matrix_table_has_one_column_per_taxon() {
        let (m, ..) = fixture();
        let t = matrix_table(&m);
        assert_eq!(t.header.len(), 6 + 3);
        assert_eq!(&t.header[6..], &["A", "B", "C"]);
        assert_eq!(t.len(), 1, "the shared position collapses to one site");
        assert_eq!(t.rows[0][3], "3", "n_taxa");
        assert_eq!(t.rows[0][5], "111", "pattern");
        assert_eq!(&t.rows[0][6..], &["1", "1", "1"]);
    }

    #[test]
    fn sites_and_taxa_tables_align_with_the_matrix() {
        let (m, ap, _proj, ps) = fixture();
        let s = sites_table(&ap, &ps);
        assert_eq!(s.len(), 3, "one row per taxon-site");
        assert!(
            s.rows.iter().all(|r| r[3] == "6"),
            "all three land on column 6: {:?}",
            s.rows
        );
        let tt = taxa_table(&ap, &ps, &m);
        assert_eq!(tt.len(), 3);
        assert!(tt.rows.iter().all(|r| r[9] == "yes"));
    }

    #[test]
    fn pairwise_and_spectrum_tables_render() {
        let (m, ..) = fixture();
        let p = pairwise_table(&m.pairwise(0, true));
        assert_eq!(p.len(), 3);
        assert!(p.rows.iter().all(|r| r[4] == "1"), "{:?}", p.rows);
        let sm = sharing_matrix_table(&m, 0);
        assert_eq!(sm.header.len(), 2 + 3);
        assert_eq!(sm.rows[0][3], "1", "diagonal holds the taxon's own count");
        let sp = spectrum_table(&m);
        assert_eq!(
            sp.rows,
            vec![vec![
                "3".to_string(),
                "1".to_string(),
                "1.000000".to_string()
            ]]
        );
    }

    #[test]
    fn scoring_tables_render() {
        let (m, ..) = fixture();
        let tree = Tree::from_newick("((A,C),B);").unwrap();
        let (data, _) = prepare(&m, &tree, &[], 1).unwrap();
        let d = score_dollo(&data, true);
        let dt = dollo_table(&d);
        assert_eq!(dt.len(), 1);
        assert_eq!(dt.rows[0][5], "1", "one gain");
        assert_eq!(dt.rows[0][7], "yes", "monophyletic");
        assert_eq!(dt.rows[0][8], "yes", "one gain, no losses");
        let mp = score_mp(&data, true);
        let mt = mp_table(&mp);
        assert_eq!(
            mt.rows[0][5], "0",
            "all three taxa carry the only site, so MP needs no change"
        );
        let ct = counts_table(&score_counts(&m));
        assert_eq!(ct.rows[0][7], "A,B,C");
        assert_eq!(ct.rows[0][5], "0", "no missing cells in this fixture");
        let rt = reconstruction_table(
            &data,
            &d.sites
                .iter()
                .map(|s| ("dollo", s.label.as_str(), s.reconstruction.clone()))
                .collect::<Vec<_>>(),
        );
        assert_eq!(rt.len(), data.tree.len(), "one row per node");
    }

    #[test]
    fn diagnostics_and_unmatched_tables_render() {
        let (_, _, proj, _) = fixture();
        let u = unmatched_table(&proj);
        assert!(u.is_empty(), "everything joined");
        let d = diagnostics_tally_table(&proj.diagnostics);
        assert!(d.is_empty() || d.rows[0][1].parse::<u64>().is_ok());
    }

    #[test]
    fn tree_with_gain_counts_annotates_the_gaining_branch() {
        let (m, ..) = fixture();
        let tree = Tree::from_newick("(((A,C)ac,B)b,((D,E)de)f);").unwrap();
        // A/C/B/D/E are in the tree; only A, B, C have data.
        let (data, rep) = prepare(&m, &tree, &[], 1).unwrap();
        assert_eq!(
            rep.tree_tip_without_data.len(),
            2,
            "D and E have no matrix row"
        );
        let d = score_dollo(&data, false);
        let nwk = tree_with_gain_counts(&data, &d);
        assert!(nwk.contains("[gains=1]"), "{nwk}");
        assert!(nwk.ends_with(';'));
        Tree::from_newick(&nwk.replace([']', '['], ""))
            .expect("still parses after stripping comments");
    }

    #[test]
    fn nexus_output_matches_the_matrix() {
        let (m, ..) = fixture();
        let s = crate::matrix::characters_nexus(&m);
        assert!(s.contains("NTAX=3 NCHAR=1"), "{s}");
        assert!(s.contains("A\t1"), "{s}");
    }

    #[test]
    fn table_writes_a_file() {
        let dir = std::env::temp_dir().join(format!("ip_tables_{}", std::process::id()));
        let (m, ..) = fixture();
        let n = matrix_table(&m).write(dir.join("sub/matrix.tsv")).unwrap();
        assert_eq!(n, 1);
        let text = std::fs::read_to_string(dir.join("sub/matrix.tsv")).unwrap();
        assert_eq!(text.lines().count(), 3, "comment + header + row");
        assert!(text.starts_with("# format=2 alignment=fam"), "{text}");
    }

    #[test]
    fn protein_site_display_and_keys_agree_with_tables() {
        let s = ProteinSite::from_cds_offset(18);
        assert_eq!(s.site_key(), "18");
        assert_eq!(s.to_string(), "6:0");
    }

    #[test]
    fn streaming_matrix_writer_is_byte_identical_and_leaves_no_temp_file() {
        let (m, ..) = fixture();
        let dir = std::env::temp_dir().join(format!("ip_stream_{}", std::process::id()));
        let path = dir.join("matrix.tsv");
        // Writing over an existing file must replace it atomically and leave nothing behind.
        crate::util::prepare_parent(&path).unwrap();
        std::fs::write(&path, "stale").unwrap();
        let n = write_matrix_tsv(&m, &path).unwrap();
        assert_eq!(n, m.sites.len());
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            matrix_table(&m).to_tsv()
        );
        let leftovers: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter(|n| n.contains(".tmp"))
            .collect();
        assert!(leftovers.is_empty(), "{leftovers:?}");
    }

    #[test]
    fn streaming_introns_writer_is_byte_identical_to_the_table() {
        let (_, _, proj, _) = fixture();
        assert!(!proj.rows.is_empty());
        let dir = std::env::temp_dir().join(format!("ip_introns_stream_{}", std::process::id()));
        let path = dir.join("introns.tsv");
        let n = write_introns_tsv(&proj, &path).unwrap();
        assert_eq!(n, proj.rows.len());
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            introns_table(&proj).to_tsv()
        );
    }
}
