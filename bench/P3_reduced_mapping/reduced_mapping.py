#!/usr/bin/env python3
"""Map a binary trait (reduced intron density) onto the genus tree and count independent origins.

Trait per tip (genus representative genome, tip = LOCUSTAG): reduced = genome-wide introns per kb CDS < THR.
Tips below the core-set floors (CDS bp < 50% of the class median, or fewer than 100 of the 453 marker Pfams found) are
missing (state unknown). Fitch parsimony on the tree gives the minimum number of state changes (independent of the
root); with missing tips all states are allowed. Maximal reduced clades are listed under a rooting on the early-diverging
flagellated phyla when they form a branch (else midpoint is not used and the clades are read on the unrooted tree).
Run: reduced_mapping.py TREE DENSITY_TSV MARKER_TSV | --selftest. Read-only."""
import sys, re, csv, collections, statistics as st

def parse_tree(nwk):
    toks = re.findall(r"\(|\)|,|;|[^(),;:]+(?::[0-9eE.+-]+)?|:[0-9eE.+-]+", nwk.strip())
    tips, children, stack, after_close = [], [], [[]], False
    node_of_tip = {}
    for t in toks:
        if t == "(":
            stack.append([]); after_close = False
        elif t == ",":
            after_close = False
        elif t == ")":
            kids = stack.pop(); children.append(kids); stack[-1].append(("n", len(children) - 1)); after_close = True
        elif t == ";":
            break
        else:
            name = t.split(":")[0]
            if name == "" or after_close:
                after_close = False; continue
            tips.append(name); stack[-1].append(("t", len(tips) - 1))
    root_kids = stack[0]
    if len(root_kids) == 1 and root_kids[0][0] == "n":
        root = root_kids[0][1]
    else:
        children.append(root_kids); root = len(children) - 1
    return tips, children, root

def masks(tips, children):
    m = [0] * len(children)
    for i, kids in enumerate(children):          # children are created before their parent: index order is post-order
        x = 0
        for k, j in kids: x |= (1 << j) if k == "t" else m[j]
        m[i] = x
    return m

def fitch(tips, children, root, state):
    """state: dict tip-> 0/1/None. Returns (changes, node sets)."""
    S = [None] * len(children); changes = 0
    def leafset(j):
        s = state.get(tips[j]); return {0, 1} if s is None else {s}
    for i, kids in enumerate(children):
        sets = [leafset(j) if k == "t" else S[j] for k, j in kids]
        cur = sets[0]
        for s in sets[1:]:
            inter = cur & s
            if inter: cur = inter
            else: cur = cur | s; changes += 1
        S[i] = cur
    return changes, S

def selftest():
    t, ch, r = parse_tree("((a,b),(c,(d,e)));")
    state = {"a": 1, "b": 1, "c": 0, "d": 1, "e": 0}
    assert fitch(t, ch, r, state)[0] == 2, fitch(t, ch, r, state)[0]
    assert fitch(t, ch, r, {"a": 1, "b": 1, "c": 0, "d": 0, "e": 0})[0] == 1
    assert fitch(t, ch, r, {"a": 1, "b": None, "c": 0, "d": 1, "e": 0})[0] == 2
    t2, ch2, r2 = parse_tree("((a:1,b:1)0.9:1,(c:1,d:1)0.8:1);")
    assert t2 == list("abcd")
    print("selftest ok")

def main(tree, dens_tsv, mark_tsv):
    tips, children, root = parse_tree(open(tree).read())
    full = (1 << len(tips)) - 1
    dens = {r["LOCUSTAG"]: r for r in csv.DictReader(open(dens_tsv), delimiter="\t")}
    mark = {r["LOCUSTAG"]: int(r["n_marker_proteins"]) for r in csv.DictReader(open(mark_tsv), delimiter="\t")}
    bycls = collections.defaultdict(list)
    for r in dens.values(): bycls[r["CLASS"]].append(float(r["cds_bp"]))
    med = {c: st.median(v) for c, v in bycls.items()}
    ok = {}
    for t in tips:
        r = dens.get(t)
        ok[t] = bool(r) and float(r["cds_bp"]) >= 0.5 * med[r["CLASS"]] and mark.get(t, 0) >= 100
    print(f"tips {len(tips)}; with density {sum(t in dens for t in tips)}; pass floors (CDS>=50% class median, >=100 markers) {sum(ok.values())}")
    # outgroup: early-diverging flagellated phyla (+Microsporidia) if they form one branch
    tax = {r["LOCUSTAG"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv"))}
    OUT = {"Chytridiomycota", "Neocallimastigomycota", "Monoblepharomycota", "Blastocladiomycota", "Microsporidia", "Rozellomycota", "Aphelidiomycota"}
    O = 0
    for i, t in enumerate(tips):
        if t in tax and tax[t]["PHYLUM"] in OUT: O |= 1 << i
    M = masks(tips, children)
    allm = set(M) | {full ^ x for x in M}
    cand = {"flagellated+Microsporidia": O,
            "Chytridiomycota+Neocallimastigomycota+Monoblepharomycota+Blastocladiomycota": sum(1 << i for i, t in enumerate(tips) if t in tax and tax[t]["PHYLUM"] in OUT - {"Microsporidia", "Rozellomycota", "Aphelidiomycota"})}
    for k, v in cand.items(): print(f"  outgroup candidate '{k}': {bin(v).count('1')} tips; is a branch: {v in allm}")
    for thr in (0.3, 0.5, 1.0):
        for label, usefloor in (("all tips", False), ("core tips only", True)):
            state = {}
            for t in tips:
                r = dens.get(t)
                if r is None or (usefloor and not ok[t]): state[t] = None
                else: state[t] = 1 if float(r["introns_per_kb_cds"]) < thr else 0
            known = [s for s in state.values() if s is not None]
            ch, S = fitch(tips, children, root, state)
            print(f"\nthreshold < {thr}/kb, {label}: tips with state {len(known)}, reduced {sum(known)}; minimum number of state changes (Fitch) = {ch}")
            if thr == 0.5 and usefloor:
                # maximal reduced clades on the tree as rooted in the file (arbitrary root) -> report clades whose known tips are all reduced
                red = sum(1 << i for i, t in enumerate(tips) if state[t] == 1)
                nonred = sum(1 << i for i, t in enumerate(tips) if state[t] == 0)
                clades = [(x, bin(x).count("1")) for x in allm if x & red and not x & nonred and x != full]
                clades.sort(key=lambda c: -c[1])
                maximal = []
                for x, n in clades:
                    if not any((x & y) == x and x != y for y, _ in maximal): maximal.append((x, n))
                print(f"   maximal reduced-only branches (either side of a split, no non-reduced known tip): {len(maximal)}")
                for x, n in maximal[:40]:
                    cls = collections.Counter(tax[tips[i]]["CLASS"] or "(none)" for i in range(len(tips)) if (x >> i) & 1)
                    print(f"     n={n:4d} reduced tips; classes: {dict(cls.most_common(4))}")
if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2], sys.argv[3])
