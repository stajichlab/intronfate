#!/usr/bin/env python3
"""S1: simulate a matrix with known truth, annotation errors and junction reads (standard library only).

Truth. The real strict-family matrix (J1) supplies the sites and each site's missing-data mask. Intron presence at every tip is simulated on the species tree
with the two-state model the likelihood fits (root presence probability, per-branch gain and loss probabilities; constant rates, taken from the J1 constant-rate
hard fit's `likelihood_branches.tsv`), so each site's number of origins (root presence + 0 -> 1 changes) is known. A site is kept only when some observable tip is
*annotated* present (the model's own ascertainment): its truth is resampled until it is.
Annotation error. Each tip's true intron is missed with probability `miss`, and a true absence is annotated present with the tip's false rate. The false rate is
set per taxon so that false introns are the scenario's share of that taxon's annotated-present cells (`false_rates`, over a universe of all positions a false intron
could land on, seen or not): one rate for all taxa would swamp taxa with few introns. Scenarios (miss, false share): none (0, 0), low (0.05, 1%),
mid (0.10, 3%), high (0.15, 8%).
Reads. For the six genomes with RNA-seq in J1, each annotated-present cell gets a read count. `empirical` draws it from the distribution of reads among the cells that
agree with RefSeq (true) or disagree (false) in J2b (`J2b_bins.tsv`, near12 level); this makes false introns look better supported than they probably are, because some
disagreeing cells are real introns RefSeq lacks. `clean` gives a false intron one read with probability 0.05 and none otherwise, true introns as `empirical`. The soft
table is `p(n) = 1 - (1 - p0) k / (n + k)` for (p0, k) = (0.8, 3), (0.8, 0.2) and (0.5, 0.2) (the last separates unsupported from supported cells most); the other three taxa (no RNA-seq) stay hard.
Output folder: matrix.tsv, truth.tsv (site, n_origins, n_true_tips), errors.tsv (the true error rates, for `--error-rates`), soft_0.8_3.tsv, soft_0.8_0.2.tsv, soft_0.5_0.2.tsv.
Run: simulate.py MATRIX BRANCHES READ_BINS OUTDIR SCENARIO READS SEED UNIVERSE | --selftest
"""
import collections, csv, os, random, sys

# (miss probability of a true intron, share of annotated-present cells that are false); the false rate per absent cell is set per taxon from that share
SCENARIOS = {"none": (0.0, 0.0), "low": (0.05, 0.01), "mid": (0.10, 0.03), "high": (0.15, 0.08)}
BIN_RANGE = {"0": (0, 0), "1": (1, 1), "2": (2, 2), "3-4": (3, 4), "5-9": (5, 9), "10+": (10, 40)}
SOFT = (("0.8", "3"), ("0.8", "0.2"), ("0.5", "0.2"))
CLEAN_FALSE_READ = 0.05


def parse_newick(s):
    """-> nested lists of tip names; a tip is a string."""
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [node()]
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            assert s[pos] == ")"
            pos += 1
            while pos < len(s) and s[pos] not in ",)":   # skip an internal label or branch length
                pos += 1
            return kids
        start = pos
        while pos < len(s) and s[pos] not in ",):":
            pos += 1
        name = s[start:pos]
        if pos < len(s) and s[pos] == ":":
            while pos < len(s) and s[pos] not in ",)":
                pos += 1
        return name

    return node()


def tips_of(t):
    return [t] if isinstance(t, str) else [x for k in t for x in tips_of(k)]


def load_branches(path):
    """-> (root presence probability, {frozenset(tips): (gain, loss)}) from likelihood_branches.tsv"""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    root, par = None, {}
    for r in rows:
        tips = frozenset(r["tips"].split(","))
        if r["branch"] == "root":
            root = float(r["gain_prob"])
        else:
            par[tips] = (float(r["gain_prob"]), float(r["loss_prob"]))
    return root, par


def draw_truth(tree, root_p, par, rng):
    """one site -> ({tip: 0/1}, n_origins)"""
    state, org = {}, 0
    top = 1 if rng.random() < root_p else 0
    org += top

    def walk(t, s):
        nonlocal org
        if isinstance(t, str):
            state[t] = s
            return
        for k in t:
            g, l = par.get(frozenset(tips_of(k)), (0.0, 0.0))
            if s == 1:
                c = 0 if rng.random() < l else 1
            else:
                c = 1 if rng.random() < g else 0
                org += c
            walk(k, c)

    walk(tree, top)
    return state, org


def read_bins(path):
    """-> {label: ([true bins], [true weights], [false bins], [false weights])} from J2b_bins.tsv (near12 level)"""
    by = collections.defaultdict(lambda: ([], [], [], []))
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["level"] != "near12":
            continue
        n, a = int(r["n_cell_introns"]), int(r["n_agree"])
        b = by[r["label"]]
        b[0].append(r["reads"]); b[1].append(a); b[2].append(r["reads"]); b[3].append(n - a)
    return dict(by)


def draw_reads(bins, is_true, mode, rng):
    tb, tw, fb, fw = bins
    if is_true or (mode == "empirical" and sum(fw) > 0):
        b = rng.choices(tb if is_true else fb, weights=tw if is_true else fw)[0]
        lo, hi = BIN_RANGE[b]
        return rng.randint(lo, hi)
    return 1 if rng.random() < CLEAN_FALSE_READ else 0


def false_rates(present, share, universe):
    """{taxon: annotated-present cells in the real matrix} -> {taxon: false rate per truly absent position}. `universe` is the number of positions a
    false intron could land on, seen or not (matrix sites + the fit's unobserved positions): a false intron at a position nobody else annotates is
    itself a matrix site, so the rate has to be per position of the whole universe, not per observed absent cell. The expected false introns are then
    `share` of the taxon's annotated-present cells: rate = share/(1-share) * present / (universe - present)."""
    return {t: (share / (1 - share)) * p / (universe - p) for t, p in present.items()}


def p_of_n(n, p0, k):
    return 1 - (1 - p0) * k / (n + k)


def simulate(matrix, branches, bins_path, outdir, scenario, mode, seed, nwk, universe):
    rng = random.Random(seed)
    miss, share = SCENARIOS[scenario]
    root_p, par = load_branches(branches)
    tree = parse_newick(open(nwk).read())
    bins = read_bins(bins_path)
    os.makedirs(outdir, exist_ok=True)
    lines = open(matrix).read().split("\n")
    header = [l for l in lines if l.startswith("#")]
    rows = [l.split("\t") for l in lines if l and not l.startswith("#")]
    cols, body = rows[0], rows[1:]
    taxa = cols[8:]
    cnt = {t: 0 for t in taxa}
    for r in body:
        for t, c in zip(taxa, r[8:]):
            cnt[t] += c == "1"
    fp = false_rates(cnt, share, universe)
    fm = open(os.path.join(outdir, "matrix.tsv"), "w")
    ft = open(os.path.join(outdir, "truth.tsv"), "w")
    fs = [open(os.path.join(outdir, f"soft_{p}_{k}.tsv"), "w") for p, k in SOFT]
    for h in header:
        fm.write(h + "\n")
    fm.write("\t".join(cols) + "\n")
    ft.write("site\tn_origins\tn_true_tips\n")
    for f, (p, k) in zip(fs, SOFT):
        f.write(f"# S1 simulated soft evidence, p0={p}, k={k}\nsite\ttaxon\tp_present\treads\tintrons\n")
    for r in body:
        mask = [c == "?" for c in r[8:]]
        for _attempt in range(2_000_000):
            truth, org = draw_truth(tree, root_p, par, rng)
            ann = {}
            for t, m in zip(taxa, mask):
                if m:
                    ann[t] = "?"
                    continue
                s = truth[t]
                ann[t] = "1" if (rng.random() >= miss if s else rng.random() < fp[t]) else "0"
            if "1" in ann.values():
                break
        else:
            raise SystemExit(f"site {r[0]}: no annotated intron in 2,000,000 draws (scenario {scenario})")
        cells = [ann[t] for t in taxa]
        r[5] = str(cells.count("1"))
        r[6] = str(sum(c != "?" for c in cells))
        r[7] = "".join(cells)
        r[8:] = cells
        fm.write("\t".join(r) + "\n")
        ft.write(f"{r[0]}\t{org}\t{sum(truth.values())}\n")
        for t in taxa:
            if ann[t] == "1" and t in bins:
                n = draw_reads(bins[t], truth[t] == 1, mode, rng)
                for f, (p, k) in zip(fs, SOFT):
                    f.write(f"{r[0]}\t{t}\t{p_of_n(n, float(p), float(k)):.6f}\t{n}\t1\n")
    with open(os.path.join(outdir, "errors.tsv"), "w") as fe:
        fe.write("taxon\tsensitivity\tfalse_rate\n")
        for t in taxa:
            fe.write(f"{t}\t{1 - miss}\t{fp[t]:.6f}\n")
    for f in (fm, ft, *fs):
        f.close()


def selftest():
    t = parse_newick("(A,(B,(C,D)));")
    assert t == ["A", ["B", ["C", "D"]]] and tips_of(t) == ["A", "B", "C", "D"]
    par = {frozenset(["A"]): (0.0, 0.0), frozenset(["B", "C", "D"]): (0.0, 0.0), frozenset(["B"]): (0.0, 0.0),
           frozenset(["C", "D"]): (0.0, 0.0), frozenset(["C"]): (0.0, 0.0), frozenset(["D"]): (0.0, 0.0)}
    rng = random.Random(1)
    st, org = draw_truth(t, 1.0, par, rng)
    assert st == {"A": 1, "B": 1, "C": 1, "D": 1} and org == 1                      # root present, no loss, no gain
    par[frozenset(["C"])] = (0.0, 1.0)
    st, org = draw_truth(t, 1.0, par, rng)
    assert st["C"] == 0 and st["D"] == 1 and org == 1                               # a loss changes no origin count
    par = {k: (1.0, 0.0) for k in par}
    st, org = draw_truth(t, 0.0, par, rng)
    assert all(v == 1 for v in st.values()) and org == 2, (st, org)                 # gains at the two children of the root: A and (B,(C,D))
    assert abs(p_of_n(0, 0.8, 3) - 0.8) < 1e-12 and abs(p_of_n(1, 0.8, 3) - 0.85) < 1e-12 and abs(p_of_n(1, 0.8, 0.2) - 0.8 - 0.2 * 1 / 1.2) < 1e-9
    fr = false_rates({"A": 100, "B": 0}, 0.1, 1100)
    assert abs(fr["A"] - (0.1 / 0.9) * 100 / 1000) < 1e-12 and fr["B"] == 0.0
    assert abs(fr["A"] * (1100 - 100) / (fr["A"] * (1100 - 100) + 100) - 0.1) < 1e-9       # false introns are 10% of A's annotated cells
    b = (["0", "1"], [0, 5], ["0", "1"], [3, 0])
    assert draw_reads(b, True, "empirical", rng) == 1 and draw_reads(b, False, "empirical", rng) == 0
    assert draw_reads(b, False, "clean", random.Random(2)) in (0, 1)
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1:2] == ["--selftest"]:
        selftest()
    else:
        a = sys.argv[1:]
        nwk = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "V6_orthology_nine_species", "species.nwk")
        simulate(a[0], a[1], a[2], a[3], a[4], a[5], int(a[6]), nwk, int(a[7]))
