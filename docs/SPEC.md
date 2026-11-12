# intronfate — design specification

Version 0.1 · Nov 2026 · status: projection + Dollo/MP/count scoring implemented and tested;
database sources and model-based scoring planned.

> `Malin/` in the working tree is third-party reference data (MALIN distribution + test datasets),
> not project source: it is git-ignored, and it is the basis of the interop/validation plan in §6b.

---

## 1. Problem

Genome annotation tells us where the introns are **in genome coordinates**. Comparative questions
about introns are asked in **protein coordinates** — "do these four orthologs have an intron at the
*same place*?" — and genome coordinates are useless for that: orthologous introns sit at completely
different chromosome positions in different species, and even orthologous *genes* shuffle positions
through indels and rearrangement.

This tool converts each annotated intron into a position on the protein its transcript encodes,
projects that position onto a multiple sequence alignment of the orthologs, and reports which
positions are shared. Shared-vs-unique is then scored under explicit evolutionary models.

**Why it matters.** The interpretation of a shared intron position depends on the model. Under
Dollo's law an intron at an exact position-and-phase is a complex character that arises once and is
lost often, so sharing implies common ancestry and *non*-sharing implies loss (or convergence). That is an
*assumption*, not a finding: parallel gain is estimated at ~8 % of shared positions overall and up to
~20 % between distant lineages such as plants and animals [R5, R6], and Dollo parsimony is known to bias
ancestral intron density downward relative to likelihood methods [R9]. Under
unrestricted parsimony, absence can equally be ancestral. The tool therefore does not pick a
philosophy — it produces one matrix and scores it several ways so the models can be compared.

### Primary use case

Given a gene family (orthologous proteins) aligned, plus a GFF3 for each species:

1. list every intron position in protein space for each member;
2. build the `site × taxon` matrix over alignment columns;
3. score sharing (Dollo gains/losses, MP steps, frequencies, pairwise conservation);
4. emit tables that are useful without re-running anything.

---

## 2. The coordinate model (the crux of the design)

An intron interrupts a coding sequence at a specific nucleotide. Walk the CDS in **transcript order**
(5'→3' along the transcribed strand, which means descending genomic coordinates on the minus strand)
and count coding bases before the intron. That count, `cds_offset`, is the canonical position: it is
strand-independent, UTR-independent, and lossless.

```
CDS (transcript order):   [ exon1: 12 nt ][ intron ][ exon2: 15 nt ][ intron ][ exon3: 12 nt ]
coding bases before:         0········12              12··········27            27········39

intron A: cds_offset = 12  → codon index0 = 12/3 = 4, phase = 12%3 = 0  (boundary: after codon 4)
intron B: cds_offset = 27  → codon index0 = 9,      phase = 0           (boundary: after codon 9)

If instead exon1 were 13 nt: cds_offset = 13 → codon index0 = 4, phase = 1
                              → the intron splits codon 5 after its first base.
```

Derived views, all emitted as separate columns so no downstream consumer has to re-derive them:

| field             | definition                                                     |
|-------------------|----------------------------------------------------------------|
| `cds_offset`      | coding nucleotides 5' of the intron (0-based)                   |
| `phase`           | `cds_offset % 3` ∈ {0,1,2}                                      |
| `codon_index0`    | `cds_offset / 3`, 0-based codon the site anchors to             |
| `codon`           | 1-based codon number used in the `codon:phase` label            |
| `residue_before`  | 1-based residue 5' of the intron (the interrupted one if phase>0)|
| `residue_after`   | 1-based residue 3' of the intron                                |
| `site_key`        | `cds_offset` as text — the stable machine key                   |

`phase 0` means the intron sits exactly at a codon boundary; `phase 1`/`2` mean it splits the codon
after 1 or 2 bases. Anchoring uses `cds_offset / 3` uniformly for every phase, which is what makes
the alignment lookup a single array index (§4).

**Convention and its edge case.** A phase-0 intron lies *between* residues but is anchored to the
*following* residue's column; phase 1/2 introns anchor to the interrupted residue. If a taxon has a gap
or indel exactly at a phase-0 boundary, its following residue can sit in a different column from its
orthologs', so the same intron lands on a different column. `--tolerance` absorbs this in pairwise
statistics only; the matrix and the Dollo/MP scoring see exact columns. MALIN uses the same convention: its alignment position `3k` lies
*between codon k and k+1*, with 1–2 inside the first codon (header of the `#MAP` block in `Chromo23.table`;
MALIN guide §2.3.3), which matches anchoring to the following residue. Whether to anchor to the preceding
residue instead remains a possible refinement for the gap-boundary case above.

### Site identity — three decisions that change the numbers

1. **The alignment column, not the residue number, is the homology unit.** Species A's residue 42 and
   species B's residue 39 can be the same position; only the alignment says so. Two introns are the
   same site iff they land on the same column.
2. **Phase is part of identity by default.** An intron at `codon 5, phase 0` and one at `codon 5,
   phase 1` are one nucleotide apart but are *different sites* — conflating them is how a
   "conserved position" signal gets manufactured. `--ignore-phase` collapses them when the
   looser definition is wanted.
3. **Drift tolerance is a statistic, not a matrix property.** Apparent intron position shifts near
   alignment gaps are common, but most are annotation artefacts rather than true sliding [R10; see also R11], so
   annotations drift by a column at indel boundaries. The matrix keeps exact columns (so it stays a faithful record); the *pairwise*
   comparison takes `--tolerance N` and matches greedily within ±N columns, so one site can never
   support two shared calls. Greedy matching is order-dependent, so counts near the tolerance can
   differ by processing order; treat tolerance results as approximate.

### Missing data is not absence

A taxon with no annotation for this gene contributes **no information**. Scoring treats its cells as `MISSING`, never `0`.
Scoring that treats "not annotated" as "intron lost" fabricates losses and is the single most
dangerous silent error in this analysis; the matrix therefore stores a per-taxon `annotated` flag and
every statistic that means anything divides by `n_annotated`, not `n_taxa`.

**Per-cell observability (implemented).** Missing data is per `(site, taxon)` cell, not per taxon. Each
site carries an `observable` bitset next to its `presence` bitset, and a cell is `Present` (`1`),
`Absent` (`0`) or `Missing` (`?`). A cell is *observable* when the taxon is annotated **and** aligned
around the site, following MALIN's rule (guide §2.3.8):

| site phase | taxon must have a residue in columns |
|---|---|
| 1 or 2 (inside the residue at column `c`) | `c-n .. c+n` |
| 0 (between residues) | `c-m .. c-1` **and** `c .. c+m-1`, with `m = max(n, 1)` |

`n` is `--flank` (default 1, an *unvalidated default* — the MALIN guide gives the rule but not a value; tune
it against `Malin/test`). Because the phase-0 columns must be adjacent, a boundary intron sitting beside an
indel is missing rather than anchored to a different column in each taxon.

Consequences worth knowing:

- A gap at the shared intron's column gives `?`, so Dollo reports **no loss** on that branch; the same taxon
  aligned through the site gives `0` and one loss. This pair is pinned by
  `crates/intronfate-core/tests/per_cell_missing.rs`.
- An annotated intron whose taxon is *not* aligned around it is **not** kept as presence. It becomes `?` and
  is counted in `masked_present` (and reported as `present_not_observable`). Keeping presence while
  discarding absence in poorly aligned regions would bias the matrix towards gains.
- An alignment row whose residues differ from the protein it was joined to is not an observation: it is
  scored missing and reported as `row_sequence_mismatch`, rather than projecting through a trimmed row.
- Matrices read from older files (no `?`, no `format=2`) still load: an unannotated taxon's `0` cells are
  read as missing. A `1` in an unannotated taxon is an error.

`matrix.tsv` is parsed by header name (`pattern` marks the end of the fixed columns). Under
`--ignore-phase` a merged site pools introns of every phase, so an *absence* needs the union of the two
windows (`c-m … c+n`) while each *carrier* is judged by the window of its own phase — judging absences by the
narrower boundary window would turn `?` into `0` and invent losses. Under `--clamp-tail` a clamped carrier
(placed on the last residue's column although its codon lies beyond the row) is exempt from the window, because
its right-hand flank is past the end of the row by construction; clamped sites are approximate and are
counted separately (`clamped_sites` in `taxa.tsv`).

---

## 3. Inputs, and the honesty problem

### Annotation → transcripts

GFF3/GTF is read as a stream and grouped by transcript (`Parent`, else `transcript_id`, else…;
attribute order is configurable with `--transcript-attrs`, because MAKER/NCBI/Ensembl/Augustus/protein2genome all disagree).
Introns are **derived**, not read: CDS features are ordered in transcript direction, merged, and each
gap becomes an intron whose `cds_offset` is the running coding length. Explicit `intron` features are
used for cross-checking rather than as the source of truth, because phase math needs the CDS.

Three independent checks run on every transcript, and disagreements are *reported*, not resolved
silently:

- **Reading-frame check.** By spec a CDS feature's `phase` is the number of bases to skip to reach the next codon boundary, `(3 - coding_bases_before % 3) % 3` (so 1 leftover base → phase 2, 2 leftover → phase 1; *not* `coding_bases_before % 3`). The
  derived offsets must reproduce the annotated frames; a mismatch means the CDS set is incomplete or
  inconsistent (`frame_mismatch`). **5'-partial models:** if the first CDS feature (in transcript order)
  has phase 1 or 2, its leading `phase` bases belong to a codon cut off at the 5' end. Protein space starts
  at the first *complete* codon, so every `cds_offset` is counted from there (offset = coding bases − skip),
  `cds_len` excludes the skipped bases, and later features are checked against `expected_frame(cum − skip)`.
  An intron inside the skipped partial codon is not projectable (`intron_before_first_codon`).
- **Transcript keys.** `Parent` is tried before `transcript_id` (GFF3), so NCBI files — whose exon rows carry
  `transcript_id=NM_x` but whose CDS rows carry only `Parent=rna-NM_x` — group under one transcript; a
  multi-valued `Parent=t1,t2` feeds every parent instead of creating a phantom transcript named `t1,t2`.
- **Exon cross-check.** Every coding gap should also be an exon gap. Exon gaps that are *not* coding
  gaps are UTR introns: kept, ranked, and explicitly non-projectable (`cds_offset = NULL`).
- **Length cross-check.** After joining to a protein, `len(protein)` must be compatible with
  `cds_len / 3` (± `--length-slack`, default 3, because annotation differs on whether the stop codon
  is part of the CDS). This is the cheapest detector of a join that landed on the **wrong paralog** —
  which otherwise produces entirely plausible-looking tables.

### Joining annotation to sequence

The fragile edge. Resolution order: `--id-map` file → `protein_id`/`protein` attribute → a
`protein`/`polypeptide` feature whose `Parent` is the CDS or transcript → transcript id (flagged as a
fallback). Lookups try exact, then version-insensitive (`ENSP000001.3` ↔ `ENSP000001`); version
stripping only splits on `.`, so `NP_001234` is never truncated to `NP`, and the transform is
idempotent because it is applied to both sides. Both sides of every failure are written to
`unmatched.tsv`: annotations with no sequence, sequences with no annotation.

### Alignments

Aligned FASTA and Stockholm. Each row precomputes `res2col` (column of its n-th residue), which
reduces projection to one indexed lookup.

Stockholm `#=GS <id> SF/SE` is honoured: those fields say the row carries **source residues SF..SE**,
so protein residue `SF-1` maps to the first aligned residue and residues outside the range map to
nothing. Treating that range as "trim these insert columns" instead of "renumber from here" shifts
every projected position in the row; the semantics are pinned by a unit test.

---

## 4. Architecture

```
GFF3/GTF ──► gff::GffReader ──► annotation::GffSource ─┐
                                                       ├─► annotation::AnnotationSource (trait)
   (future: SQL/Chado/Ensembl reader) ─────────────────┘          │ transcript-level structure
                                                                  ▼
protein FASTA ─► seqio::ProteinSet ──────────► project::Joiner / JoiningSink (streamed)
                                                  │  protein_id resolution + length checks
                                                  ▼
                                              Projection (per-protein site sets)
                                                                  │
MSA ─► msa::Alignment (residue↔column maps) ──► project::project_into_alignment
                                                  │
                                                  ▼
                                          matrix::SiteMatrix  (column×phase × taxon, bitsets)
                                             │           │            │
                              pairwise/──────┘  score::prepare ───────└──► tables::*  (TSV/NEXUS)
                              Jaccard                │
                                              score::{score_dollo, score_mp, score_counts}
                                                     ▲
                                     tree::Tree ◄────┘  (Newick, --outgroup re-rooting)
```

### The seam

```rust
pub trait AnnotationSource {
    fn describe(&self) -> String;
    fn stream_into(&mut self, sink: &mut dyn AnnotationSink) -> Result<()>;
}
pub trait AnnotationSink {
    fn push_transcript(&mut self, tx: TranscriptAnnotation) -> Result<()>;
    fn diagnostic(&mut self, level: Level, kind: &'static str, subject: String, message: String);
}
```

Visitor-shaped rather than `Iterator` because both real sources are push-based: GFF can only
finalise a transcript at EOF, and a SQL cursor hands rows over as they stream. A future
`SqlSource` implements the same trait and touches nothing downstream. `AnnotationSet` is the
in-memory sink; `MemorySource` is the cached source, which also makes fixtures cheap.

### Modules

| module          | responsibility                                                       |
|-----------------|----------------------------------------------------------------------|
| `gff`           | byte-level GFF3/GTF line parsing, zero allocation in the hot loop     |
| `annotation`    | grouping by transcript, CDS ordering/merging, intron + offset derivation, the checks in §3 |
| `seqio`         | protein FASTA, dual (exact/version-stripped) id index                  |
| `msa`           | alignment rows with `res2col` maps and SF/SE windows                   |
| `project`       | annotation↔sequence join, protein site sets, projection to columns     |
| `matrix`        | `SiteMatrix`, bitsets, conservation spectrum, pairwise sharing, NEXUS  |
| `tree`          | Newick parse, MRCA, monophyly, edge-splitting re-rooting               |
| `score`         | Dollo, weighted-Dollo, Fitch MP, counts, missing-data handling         |
| `tables`        | TSV/NEXUS serialisation of everything                                   |
| `report`        | diagnostics with counts, capped retained rows                          |

Workspace: `crates/intronfate-core` (the library, publishable) and `crates/intronfate-cli` (the
`intronfate` binary). Dependencies are deliberately few: `rayon`, `flate2`, `rustc-hash`,
`thiserror`, plus `clap`/`anyhow` in the CLI. Parsers are hand-written because the formats are
line-oriented and control over tolerance/performance matters more than reuse.

---

## 5. Scoring models

Each site is a binary character over taxa. What differs between models is the cost scheme.

### Dollo parsimony (implemented — classic form)

Assumption: a position-with-phase arises **once** and is never regained after loss.

The gain is pinned at the **MRCA of the taxa carrying the site** — any deeper placement can only add
losses, any shallower one cannot cover all carriers. Below the MRCA a two-state irreversible DP
minimises losses (`present → lost` allowed, `lost → present` forbidden):

```
costP(v) = Σ_children min( costP(c), costL(c) + 1 )     // +1 = a loss on the edge to c
costL(v) = Σ_children costL(c)                          // lost lineages stay lost
tips:  present → (0, ∞)   absent → (∞, 0)   missing → (0, 0)
losses = costP(mrca)
```

Reported per site: `gains` (always 1 under classic Dollo, so `dollo_gains` equals the site count and carries no information — the loss count is the result), `losses`, `monophyletic` (carriers are exactly one clade),
`compatible` (`losses == 0`), `gained_at`, `lost_in`, and optionally the full per-node
reconstruction. `tree_gains.nwk` puts `[gains=N]` on the branches where origins land.

**A degeneracy worth knowing.** The tempting formulation — "3-state DP, gain cost 1, loss cost 0,
minimise gains" — is *useless*: one gain at the root plus free losses explains any pattern, so the
minimum is 1 for every site. This was implemented first, and caught only because a test asserted
that three independent islands need more than one event. Under classic Dollo the informative
statistic is the **loss count**, and the only way to get a real gain/loss trade-off is weights (§5.2).

### 5.2 Weighted Dollo (implemented, for sensitivity analysis)

Three states — `U` not yet gained, `P` present, `L` lost-and-terminal — with a gain costing
`gain_weight` and a loss `loss_weight`; regaining a lost position is structurally impossible.

```
U(v) = Σ min( U(c), gw + P(c) )      P(v) = Σ min( P(c), lw + L(c) )      L(v) = Σ L(c)
```

With `gw=lw=1`, two independent origins (cost 2) beat one origin plus two losses (cost 3); with
`gw=3, lw=1` the deep single origin wins. Sweeping this ratio is how you show a conclusion is not an
artefact of weighting.

### 5.3 Maximum parsimony (implemented)

Fitch down-pass/up-pass, symmetric unit cost — the unrestricted comparison for the same matrix. Nodes
with more than two children (an unrooted Newick has a trifurcating root) use Hartigan's counting rule: a
node takes the state(s) held by the most informed children and every other informed child costs one step.
Intersecting children pairwise is only right for binary nodes — children `{0},{0},{1},{1},{1}` need 2 steps,
not 1. Missing tips are free.
Reports per-site steps, parsimony-informativeness, ambiguity, the total length and consistency
index, plus optional ancestral states.

### 5.4 Counts (implemented)

No tree: per-site `n_present`, `n_absent`, frequency, carrier list; per-pair `shared`, Jaccard,
`shared/min(n_a,n_b)` (commonly reported as a fraction of conserved positions; no single source for the exact
denominator was found, so it is stated here as this tool's definition); the conservation
spectrum (how many sites are shared by exactly *k* taxa).

### 5.5 Planned

- **Continuous-time gain/loss likelihood** on branch lengths (ML, and a likelihood-ratio test of
  "ancient intron" vs "convergent insertion"), which needs `Tree` branch lengths — already parsed.
- **Rate heterogeneity across sites** (Γ over positions); intron gain rates vary strongly
  between lineages [R5, R12], so a homogeneous model over-reads losses.
- **Bias correction**: positions in poorly aligned or gappy columns are under-called. Plan:
  per-column observability from the alignment and an explicit "scorable site" set.
- **Probabilistic projection**: propagate annotation uncertainty into the sharing matrix rather than
  treating every boundary as exact.

---

## 6. Performance

Stated target: thousands of proteins and genome-scale GFF, comfortably multi-species.

- Streaming I/O throughout; gzip is transparent (`.gz`/`.bgz` and magic-byte detection), so
  multi-GB annotation never lands on disk uncompressed.
- GFF parsing allocates nothing per line: attributes are byte spans into a reused buffer.
- Interned indices (`u32`) instead of strings in inner loops; presence sets are bitsets, so
  intersections are popcounts.
- Per-transcript finalisation, per-site scoring and per-character Fitch are embarrassingly parallel
  and run under `rayon` (`into_par_iter`, `par_iter`).
- Projection into a row is one array index per site; matrix construction is one pass over projected
  sites.

**Complexities.** Parse O(features); finalise per transcript O(s log s) for the CDS sort; project
O(S log R) per taxon (S sites, R residues); matrix build O(S log S); Dollo/MP O(S · N) for N tree
nodes; pairwise O(t² · S̄) for t taxa.

### Measured

Scoring the intron table distributed with MALIN as `rogozin_introns_table` (derived from the 684-KOG,
8-genome analysis of Rogozin et al. 2003 [R4]; the Ecdysozoa tree comes from Roy & Gilbert 2005 [R7]) (7 taxa × 7 236 alignment columns → 7 221 scorable
sites) from a `matrix.tsv`: convert 0.04 s, score (counts + Dollo + MP + all writers) 0.05 s. The 23-taxon
`Chromo23.table` over 405 alignments: import 0.74 s, scoring all 405 alignments one process each 8.3 s, export
0.06 s. **Measured on the HPCC** (job 32011208, `epyc` partition, release build, 2026-11-11; full report in
`bench/V1_malin_validation/REPORT.2026-11-11.txt`); these are wall times of the whole process including
start-up. Parse vs project vs score is not separated within a run (only whole-process wall time exists).

**Projection on real annotation (V2, HPCC `highclock`, 2026-11-11, `bench/V2_gff_projection/REPORT.2026-11-11.txt`).**
Nine fungal genomes (funannotate GFF3, 5.9–74 k CDS rows each): GFF parse + protein-space projection takes
0.05–0.32 s and 18–125 MB peak RSS per genome. 739 BUSCO ortholog families (6 363 protein rows) through
`analyze` in one process: 8.4 s and 507 MB (a binary built before the per-alignment memory release, parallel
pairwise, O(k) `mrca` and streaming-writer fixes took 34 s and 3.7 GB on the same data — two changes, not
isolated). **Genome scale (V3, same node type, `bench/V3_gff_scale`).** The largest GFF3 among 23 600 funannotate annotations is
14 MB gzipped / 165 MB plain (*Phakopsora pachyrhizi*, ~108 k transcripts, ~520 k introns): `project` takes ~2.0 s
and ~400 MB (parse 1.2 s ≈ 140 MB/s plain; gzip decoding is ~11 % of parse). Wall time does not improve with threads
(1→16: 1.00→1.02×; parsing is single-threaded). 10 / 100 / 1 000 genomes in one process: 1.8 s, 0.3 GB / 19.7 s,
2.9 GB / 226 s, 32.5 GB (10.9 M transcripts, 23.4 M introns) — linear in time, memory-bound; before the streaming
`introns.tsv` writer and early transcript release the 1 000-genome run needed 55.8 GB. After files were read in
parallel and the annotation streamed into the join, it takes **52 s and 7.0 GB** (0.64 KB per transcript; with 100
genomes 3.5 s from 12.3 s as threads go 1→16). Still open: a mammalian-sized
GFF (the BFD has none), `matrix`/`analyze` scaling in taxa and sites, and a synthetic generator.

---

## 6b. Relationship to MALIN (and the interop plan)

MALIN (Csűrös 2008 [R1]; BSD-licensed) is a closely related prior tool and ships test data in this working tree
(`Malin/test/`). Comparing the two clarifies both what is new here and what must be compatible.

**Different input side.** MALIN consumes intron positions that someone else already derived, carried
in FASTA headers:

```
>gi|124423860|emb|CAK88653.1| … [Paramecium tetraurelia] /organism=Ptet … {i 384,1535,1584 i} /intron-length=28,23,37
```

`{i … i}` lists the intron positions on the source sequence and `/organism=` is the taxon key. This
tool's contribution is on the other side of that arrow: deriving those positions from a **GFF3/GTF**
by CDS arithmetic, with the frame/exon/length cross-checks of §3. MALIN itself does not read genome
annotation. Other tools do work on that side — **GenePainter** [R13] maps intron positions and phase onto
protein alignments (its documentation lists GFF input), and **CIWOG** [R14] and **GECA** [R15] build intron
sets across orthologs — so this tool's novelty claim should be made against those as well as MALIN, and a
head-to-head comparison with GenePainter is worth adding to §6b's validation plan.

**Matrix encodings.** MALIN's `Chromo23.table` has one row per taxon and one character per intron
*site*, where a site is a nucleotide-resolution alignment position: a `#MAP site kog alipos` block maps
each site to a KOG alignment and a position (`3k` is between codon k and k+1; 1–2 are inside the first
codon). The alphabet is `0` (no intron), `1|2|3` (intron present; `3` = phase 0, i.e. at a codon boundary,
`1` and `2` = phase 1 and 2, the phase being fixed by `alipos mod 3`), and `?` (missing). The guide states
that 1, 2 and 3 are currently all treated as "present". This tool writes one row per `(column, phase)` and
one `0/1` column per taxon. A MALIN site is therefore already a (column, phase) pair, and the two
encodings correspond directly: `alipos = 3·codon_index0 + phase` (checked against the first three `#MAP`
rows of `Chromo23.table`; a full cross-check is part of the validation plan). A taxon has its own cell in
each site row, so two introns of one taxon on one column with different phases (very short exons, or
`TailPolicy::Clamp`) are two separate MALIN sites, not a conflict. What does not round-trip yet is the `?`
per-cell state (next section) and the KOG/`#MAP` bookkeeping.

*Not re-verified:* alphabet frequencies for `Chromo23.table` (an earlier version of this spec quoted
counts and described phase 2 as common; both were wrong or unchecked). Phase 0 is the commonest class in
the literature, and phase 2 the rarest [R16]. Recount on the HPCC before quoting any figure.

### Missing data in MALIN tables (now supported on the matrix side)

MALIN's tables mark `?` **per (taxon, site)** (the MALIN guide, §2.4.1, ties it to a missing ortholog or
alignment uncertainty). `SiteMatrix` now stores that per cell (§2) and `matrix.tsv` round-trips `?`.
`scripts/malin_to_matrix.py` converts a MALIN table with `#MAP` lines into one `matrix.tsv` per alignment,
using `(column, phase) = (alipos // 3, alipos % 3)` and the phase codes `3`→0, `1`→1, `2`→2; it rejects a table
whose phase code contradicts its `alipos`. MALIN also requires a minimum number of unambiguous sequences per
site before calling presence/absence (guide §2.3.8); the corresponding `--min-observed` filter is not
implemented yet.

A related correction made on the way: "annotated" used to be derived from "has introns", so a gene that lost *every*
intron was scored as missing data — discarding exactly the negative evidence intron-loss analysis is made of.
`Projection` now records annotated proteins independently of intron count (`n_transcripts_of` /
`is_annotated`); an annotated intronless gene is an observation of absence.

### Validation plan

1. `rogozin_introns_table` (7 taxa × 7 236 columns, alphabet `0/1` only, so the conversion is exact)
   converted to `matrix.tsv` and scored on `Rogozin-ecdysozoa.tre`: 7 221 sites, 7 221 gains,
   **2 770 independent losses**, 1 619 positions (22 %) not explainable without at least one loss,
   9 007 MP steps — **re-derived on the HPCC and confirmed** (`bench/V1_malin_validation`, 2026-11-11). The source analysis [R4] concluded that many ancestral introns were lost in fungi, nematodes and
   arthropods (not in Ecdysozoa as a clade); the Ecdysozoa-vs-Coelomata reading of this table is
   contested [R7, R8]. **Cross-check against MALIN itself:** the guide (p. 27) shows MALIN's own Dollo run on
   this table. All nine per-branch loss counts visible there — Sp 796, Ecdysozoa 568, Diptera 436, Ce 411,
   Ag 194, Dm 156, Hs 112, At 46, Pf 0 — are reproduced **exactly** by `dollo_branches.tsv`. The visible rows
   sum to 2 719; the remaining 51 losses are on two branches the guide's page does not show (Bilatera 27,
   Opisthokont 24), which gives our total of 2 770. This agrees with MALIN's Dollo on branch counts, not just
   on the total (it is evidence that classic Dollo with the one-gain-at-MRCA rule matches MALIN's; it says
   nothing about the likelihood models, which are not implemented). Site-by-site agreement with MALIN's own
   output files has not been compared.
2. Native `--from-malin` / `--format malin` (the Python converter covers import today). Per-cell missing
   data has landed and the import is native: `Chromo23.table` (23 taxa, 9 256 sites over **405** alignments in
   its `#MAP` block — an earlier draft of this spec said 394, which was wrong) imports, scores (8 686 of the
   9 256 sites have a carrier and are scored: 8 686 gains, 6 475 losses, 11 760 MP steps) and exports back with
   every `(alignment, alipos)` row identical. Its alphabet counts (`0` 160 418, `?` 33 102, `3` 10 493,
   `2` 4 568, `1` 4 307) were also re-counted and match.
3. Compare against `*.rates` (MALIN's maximum-likelihood gain/loss rate estimates, produced from
   `table + newick`) once §5.5 exists — that is the direct quantitative cross-check on the model side.
4. Cross-check the projection: run this tool's GFF path over a species whose MALIN `.faa` headers
   carry known positions, and diff derived `cds_offset` against theirs.

## 7. Outputs

All TSV: `# key=value` metadata first line where the format has nowhere else for it, one header row,
absolute coordinates, no wrapped lines, stable ordering (diffable between runs).

| file                       | contents                                                            |
|----------------------------|---------------------------------------------------------------------|
| `introns.tsv`              | every coding intron: genomic coords + `cds_offset`/`phase`/`codon`/`residue_before`/`residue_after`/`site_key` |
| `unmatched.tsv`            | both sides of every failed join, with a reason code                  |
| `diagnostics.tsv` / `_detail.tsv` | counts by kind, then the individual records                  |
| `<aln>/sites.tsv`          | long form: taxon × projected site (keeps the protein-space audit)    |
| `<aln>/taxa.tsv`           | per taxon: annotated?, aligned residues, sites placed, tail dropped  |
| `<aln>/matrix.tsv`         | the binary `site × taxon` matrix + `pattern` string; re-readable by `score` |
| `<aln>/presence.tsv`       | presence-only long form                                              |
| `<aln>/pairwise.tsv`, `sharing_matrix.tsv` | sharing counts, Jaccard, fraction-of-smaller, `n_both_observable` |
| `<aln>/spectrum.tsv`       | sites shared by exactly *k* taxa                                     |
| `<aln>/counts_sites.tsv`   | tree-free frequency table                                            |
| `<aln>/dollo_sites.tsv`    | gains, losses, monophyletic, compatible, gained_at, lost_in          |
| `<aln>/mp_sites.tsv`       | steps, informative, ambiguous                                        |
| `<aln>/reconstructions.tsv`| per-node ancestral states (with `--reconstructions`)                 |
| `<aln>/tree_gains.nwk`     | Newick with `[gains=N]` annotations                                  |
| `<aln>/characters.nex`     | NEXUS discrete-character block for external phylogenetics            |
| `combined/*`               | with `--combine`: all families stacked into one matrix (`family`, `family_column`) and scored once; see below |

**Combined matrix (`--combine`, `SiteMatrix::combine`).** A run over many gene families needs the totals per lineage; summing per-family tables
afterwards is error-prone, and the right unit for it is one matrix. Design decisions: (1) Dollo, Fitch and counts are per-site, so scoring the stacked
sites once equals summing the families; this is asserted by a test on every branch. (2) A site's identity is *(family, column, phase)*, never
*(column, phase)*: the same column in two families is two sites. To keep `matrix.tsv` readable by the existing reader (which sorts and merges by
`(column, phase)`), every family's columns are shifted by its alignment length plus a pad of 10,000, so the stored `column` is unique and ascending, the
true column is a separate `family_column`, and pairwise matching at a tolerance below the pad can never pair sites of different families (a larger
tolerance is refused). (3) Taxa are the union over families, and a taxon a family lacks is `?` at that family's sites: missing, never absent. (4) Families
must share phase mode and flank and have distinct names; anything else would make sites incomparable or ambiguous, so it is an error. The combined
matrix is held in memory; the per-family folders are still streamed and written.

`matrix.tsv` is the pipeline's interchange format: `score` reads it back (including the metadata line
that says which taxa were unannotated), so re-scoring under different rootings, tolerances or weight
schemes never re-parses a GFF.

---

## 8. CLI

```
intronfate project --gff sp.gff.gz -p sp.faa -o out                 # protein space only
intronfate matrix  --gff sp.gff.gz -p sp.faa --msa fam.faa -o out   # + alignment matrices
intronfate score   --matrix out/fam/matrix.tsv -t fam.nwk --outgroup SP_X
intronfate analyze --gff ... -p ... --msa ... -t ... --reconstructions
intronfate convert  --from-malin t.table -o dir                      # MALIN table <-> matrix.tsv
intronfate demo     -o intron_demo                                   # built-in worked example
```

Multiple `--gff` / `--proteins` / `--msa` values are accepted for multi-species runs; ids must be
globally unique; a transcript id repeated across annotation files is kept as separate records and reported (`transcript_id_collision`) — files are independent and read in parallel, which also bounds parse memory to one file per worker thread. `--ignore-phase`, `--clamp-tail`, `--tolerance`,
`--min-taxa`, `--length-slack none`, `--strict`, and the attribute-configurable join keys expose
every modelling decision made above. `demo` is a self-check: it prints three sites whose expected
Dollo/MP values are asserted by the test suite.

---

## 9. Correctness strategy

245 tests pass (180 core unit, 23 + 4 core integration in `crates/intronfate-core/tests/per_cell_missing.rs` and `combined.rs`, 10 CLI unit, 25 CLI end-to-end in `crates/intronfate-cli/tests/cli.rs`, 3 doc) plus 3 Python converter tests (`scripts/test_malin_to_matrix.py`); last verified 2026-11-12. The layers are described in [`TESTING.md`](TESTING.md). Layers:

1. **Property/consistency tests on the coordinate model.** `cds_offset == codon_index0*3 + phase` for
   every offset; every derived view agrees with the emitted table columns.
2. **Brute-force cross-check of the hard DP.** For five topologies × all 2^N presence patterns, the
   Dollo DP's loss count is compared against exhaustive enumeration over every legal
   origin-placement and every irreversible present/lost labelling. This validates both the DP and the
   "place the gain at the MRCA" claim, and it is what would catch a future edit to the recurrences.
3. **Reconstruction self-consistency.** The traceback's gains/losses are recounted from the emitted
   node states, and violations (loss before gain, regained lost state, tip/state mismatch) panic.
   Reconstructed MP states must realise exactly the reported step count.
4. **Round-trip tests.** `matrix.tsv` → `SiteMatrix` → identical statistics; Newick → text → tree.
   Re-rooting preserves every tip, all edge lengths, and the MRCA relations it should not touch.
5. **Adversarial input fixtures.** GTF quoting and `=` inside quoted values, percent-escapes, CRLF,
   missing columns, malformed lines (counted vs fatal under `--strict`), 7-column rejection,
   overlapping/contiguous CDS, missing strand, transcripts with no CDS, colliding ids across files,
   ragged and duplicate-id alignments, all-gap rows, multi-block Stockholm.
6. **A worked biological fixture** (`demo`) with hand-computed expected gains, losses and MP steps,
   asserted by tests — the numbers in §10 are test assertions, not prose.

### Bugs the tests caught (each now has a regression test)

- **Degenerate Dollo**: minimising gains with free losses always yields 1 → replaced with the
  classic one-gain/minimise-losses formulation (§5.1).
- **Minus-strand CDS merging** compared segments in genomic order, so every minus-strand intron
  looked "contiguous" and vanished → merging now happens in transcript-space coordinates.
- **Re-rooting at an outgroup tip** made that tip an internal node, silently dropping it from the
  analysis → re-rooting splits the *edge* and inserts a new root, keeping all tips.
- **`strip_version`** split on `_`/`-`, truncating `NP_001234` to `NP` → `.` only, idempotent.
- **Dropped transcripts lost their diagnostics**, because finalisation returned only
  `Option<TranscriptAnnotation>` → `(Option<…>, Vec<Diagnostic>)`.
- **Tip order** was reversed (a postorder already lists tips left-to-right), which would have
  permuted characters against tree tips.

Found by the multi-agent review on 2026-11-11 (each has a regression test):

- **Fitch undercounted on polytomies** (pairwise intersection of children); now Hartigan counting, checked
  against brute force over every `{0,1,missing}` tip pattern on eight tree shapes including polytomies. The
  Dollo DP is checked the same way.
- **Observability was never wired**: `project_into_alignment` left the per-column residue sets empty, so
  gaps would still have scored as losses. Now populated, with `None` meaning "hand-built projection".
- **5'-partial CDS** shifted every `cds_offset` by the first feature's phase (see §3).
- **NCBI exon/CDS split** and **`Parent=t1,t2` phantom transcripts** (see §3).
- **`pair_stat` missed shared sites** when one taxon had two phases at a column (matching is now per phase).
- **`strip_version` was not idempotent** although applied to both sides of a join; it now removes every
  trailing `.<digits>` segment. The version-insensitive fallback refuses a key shared by several sequences
  instead of silently joining the first.
- **`TailPolicy::Clamp`** snapped non-tail sites (before a Stockholm window, or on an all-gap row) onto the last
  column.
- **Annotated intronless proteins** were still listed as "no annotation" in `unmatched.tsv`.
- **Panics / aborts** on a `%` escape beside a multibyte character, on a Newick error snippet cut inside a
  multibyte label, and on one Latin-1 byte anywhere in a GFF/FASTA/alignment (now read lossily and counted).
- **Newick**: `[...]` comments, quoted UTF-8 labels, numeric internal labels (support values) used as branch
  names, multi-taxon outgroups that are the complement of a clade, and a unary node left behind by re-rooting.
- **`dollo_weighted`** lost a gain at very large weights (an "impossible" sentinel smaller than a real cost).
- **`matrix.tsv` order** was assumed, not checked; rows are now sorted and merged on read.
- **Annotation diagnostics were never reported** (found while restructuring the loader for memory): the CLI kept
  only the transcripts and dropped the reader's diagnostics, so `frame_mismatch`, `no_cds`, `strand_unknown`,
  `protein_id_fallback_transcript`, `intron_before_first_codon`, `transcript_id_collision` and `invalid_utf8` never
  reached `diagnostics.tsv` (the usage guide documented several of them). Fixed; the annotation is now streamed into
  the join (`JoiningSink`), which carries its diagnostics.
- **Second review round (2026-11-11, three Opus reviewers + a mutation run):**
  - a gain on an internal node labelled `root` was drawn on the real root (names are now unique, gains are
    indexed by node id);
  - `--ignore-phase` judged absences by the narrower boundary window, turning `?` into `0` and inventing losses
    — absences now use the union window and carriers their own phase's;
  - `--clamp-tail` could never keep a phase-1/2 site at the default flank (clamped carriers are now exempt);
  - `mp_sites.ambiguous` was true whenever any Fitch set was `{0,1}`; it now means several optimal labellings,
    counted exactly;
  - `convert --from-malin` could overwrite or write outside `-o`; taxon names with commas/spaces broke the
    `unannotated=` list; isoforms with nested intron sets were not flagged; MALIN export of an empty matrix.
- **MALIN converter** mapped phase codes one off (`3`→phase 2) and treated table positions as alignment columns.

---

## 10. Worked example (the built-in demo, asserted by tests)

Five orthologs, 15-column alignment; intron positions at CDS offsets
species1 `{12,27}`, species2/3 `{12,18}`, species4 `{27}`, outgroup `{12,27}`, all phase 0. Projected:

| site  | column | phase | carriers          | reading                                    |
|-------|--------|-------|-------------------|--------------------------------------------|
| `4:0` | 4      | 0     | s1,s2,s3,outgroup | ancestral, lost in species4                 |
| `6:0` | 6      | 0     | s2,s3             | derived insertion on the s2/s3 stem         |
| `9:0` | 9      | 0     | s1,s4,outgroup    | ancestral, lost once in the s2/s3 ancestor  |

Rooted on `outgroup`: **3 Dollo gains, 2 losses**; 2 of 3 sites need at least one loss;
**3 MP steps** — cheaper per site here only because the losing taxa happen to be sisters. This is the
comparison the tool exists to make mechanically repeatable.

---

## 11. Assumptions and limitations

- **CDS features define protein space.** UTR introns are reported but never projected. With no CDS
  at all, `--derive-cds-from-exons` assumes translation starts at the first exon base and flags it.
- **One protein per sequence id.** Isoforms sharing a protein id are unioned and flagged
  (`isoform_intron_conflict`); they are not modelled as separate entities.
- **Exact-name taxon matching between matrix and tree.** Fuzzy matching here would place introns on
  the wrong species, so mismatches are errors/warnings instead.
- **Alternative assemblies/annotation versions** are out of scope; the annotation and the proteins
  must correspond, and the length check is the guard rail.
- Projection into an alignment inherits the alignment's errors. No alignment is inferred here.
- GFF3/GTF and aligned FASTA/Stockholm in; no BAM/spliced-alignments, no GenBank `join()` CDS
  (`gff`-only for now), no variation of the genetic code.
- Branch lengths are parsed but used only for re-rooting; no model uses them yet.

## 12. Roadmap

1. ~~Per-cell missing data~~ — done (§2).
2. **MALIN interop** — table import/export is done (`convert`, `malin.rs`, §6b). Still open: reading
   `/organism=` + `{i … i}` FASTA headers as an `AnnotationSource`, which would make the two tools directly
   comparable on MALIN's own inputs, and running `Chromo23.table` as a regression fixture on the HPCC.
3. **Benchmark + synthetic generator** — scale taxa/sites independently; time the GFF path on a real
   genome (§6).
4. **`AnnotationSource` for SQL/Chado/Ensembl** — the seam exists; needs a schema decision.
5. **Likelihood-based scoring** (gain/loss CTMC on branch lengths, LRT, rate heterogeneity),
   validated against MALIN's `*.rates` output (§6b).
6. ~~Column observability~~ — done (§2: per-cell observability with a flank rule), as are `--min-observed` and
   the `n_both_observable` column in `pairwise.tsv`. Remaining: tuning the `--flank` default.
7. **Stdin/stdout support** (`--gff -`, `-o -`) for pipeline composition.
8. **Python bindings** returning numpy arrays, if exploration in pandas becomes the bottleneck.

## 13. References

Entries marked ✓ were checked against PubMed (title, authors, journal, year, DOI) on 2026-11-11.

- **R1** ✓ Csűrös M. Malin: maximum likelihood analysis of intron evolution in eukaryotes. *Bioinformatics* 24(13):1538–9 (2008). doi:10.1093/bioinformatics/btn226
- **R2** MALIN User's Guide, M. Csűrös, 10 Apr 2008 (`Malin/malin_guide.pdf`). Section numbers cited above are the guide's own.
- **R3** Csűrös M. Count: evolutionary analysis of phylogenetic profiles with parsimony and likelihood. *Bioinformatics* 26(15):1910–2 (2010). doi:10.1093/bioinformatics/btq315 (PMID 20551134; not intron-specific; not re-checked).
- **R4** ✓ Rogozin IB, Wolf YI, Sorokin AV, Mirkin BG, Koonin EV. Remarkable interkingdom conservation of intron positions and massive, lineage-specific intron loss and gain in eukaryotic evolution. *Curr Biol* 13(17):1512–7 (2003). doi:10.1016/s0960-9822(03)00558-x
- **R5** ✓ Carmel L, Rogozin IB, Wolf YI, Koonin EV. Patterns of intron gain and conservation in eukaryotic genes. *BMC Evol Biol* 7:192 (2007). doi:10.1186/1471-2148-7-192
- **R6** Sverdlov AV et al. *Nucleic Acids Res* (2005), PMID 15788746, doi:10.1093/nar/gki316 (cited by the review agent; not re-checked).
- **R7** ✓ Roy SW, Gilbert W. Resolution of a deep animal divergence by the pattern of intron conservation. *Proc Natl Acad Sci USA* 102(12):4403–8 (2005). doi:10.1073/pnas.0409891102
- **R8** ✓ Zheng J, Rogozin IB, Koonin EV, Przytycka TM. Support for the Coelomata clade of animals from a rigorous analysis of the pattern of intron conservation. *Mol Biol Evol* 24(11):2583–92 (2007). doi:10.1093/molbev/msm207
- **R9** ✓ Csűrös M, Rogozin IB, Koonin EV. A detailed history of intron-rich eukaryotic ancestors inferred from a global survey of 100 complete genomes. *PLoS Comput Biol* 7(9):e1002150 (2011). doi:10.1371/journal.pcbi.1002150
- **R10** ✓ Sêton Bocco S, Csűrös M. Splice sites seldom slide: intron evolution in oomycetes. *Genome Biol Evol* 8(8):2340–50 (2016). doi:10.1093/gbe/evw157
- **R11** Stoltzfus A et al. *Proc Natl Acad Sci USA* (1997), PMID 9380704, doi:10.1073/pnas.94.20.10739 (intron sliding; not re-checked).
- **R12** ✓ Rogozin IB, Carmel L, Csűrös M, Koonin EV. Origin and evolution of spliceosomal introns. *Biol Direct* 7:11 (2012). doi:10.1186/1745-6150-7-11 (review; also Roy & Gilbert 2005 *PNAS* doi:10.1073/pnas.0500383102, not re-checked).
- **R13** GenePainter: PMID 23496949 (doi:10.1186/1471-2105-14-77); v2 PMID 25434742 (doi:10.1093/bioinformatics/btu798). Not re-checked; GFF-input claim is from its help page.
- **R14** Wilkerson MD et al. CIWOG. *Brief Bioinform* (2009), PMID 19933210, doi:10.1093/bib/bbp051 (not re-checked).
- **R15** Fawal N et al. GECA. *Bioinformatics* (2012), PMID 22467908, doi:10.1093/bioinformatics/bts153 (not re-checked).
- **R16** Fedorov A et al. *Nucleic Acids Res* (1992), PMID 1598214, doi:10.1093/nar/20.10.2553 (phase distribution; not re-checked).
- **R17** ✓ Roy SW, Gilbert W. The pattern of intron loss. *Proc Natl Acad Sci USA* 102(3):713–8 (2005). doi:10.1073/pnas.0408274102 (phase-zero introns more likely lost; phase 0 = does not interrupt a codon).
- **R18** ✓ Csűrös M, Rogozin IB, Koonin EV. Extremely intron-rich genes in the alveolate ancestors inferred with a flexible maximum-likelihood approach. *Mol Biol Evol* 25(5):903–11 (2008). doi:10.1093/molbev/msn039 (the Chromo23 dataset; attribution from the review agent).
- **R19** GFF3 specification, The Sequence Ontology (phase = bases to remove to reach the next codon start). https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

**Not verified at all:** ROPE, Iwamoto/Hoshiyama-style intron databases, IntronDB and IGLAS were named in the
review brief but not found; none is cited here. A 2022 Communications Biology paper on maximum-likelihood
intron mapping was seen only as a search snippet and is not cited.

**Related methods worth reading when the likelihood models (§5.5) are built:** Roy & Gilbert 2005 EM
(PMID 15687506); Nguyen et al. 2005 (PMID 16389300); Csűrös, Holey & Rogozin 2007 (PMID 17646350);
Stajich, Dietrich & Roy 2007 (PMID 17949488); Poverennaya et al. 2020 (PMID 33308147). PMIDs only, from the
review agent, unchecked.
