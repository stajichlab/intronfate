"""V13 core: read intronfate outputs and trees, map *Antarctolithica* (Alith) intron positions onto gains, losses and placements.

No plotting here (see v13_figures.py), so every function can be fixture-tested (tests/test_v13.py).

Conventions
- A node is identified by the frozenset of its tip names. intronfate's `dollo_branches.tsv` gives every node's full tip list, so its
  (possibly truncated) branch names in `dollo_sites.tsv` (`gained_at`, `lost_in`) and in `reconstructions.tsv` (`node<A-B-C-D…(n tips)>`)
  are translated to tip sets through that table.
- A lineage made of several strains (Alith) is one unit: present if any strain carries the position, absent if every observed strain
  lacks it, missing if no strain is observed.
"""
import csv
import random


# ---------------------------------------------------------------- trees
def parse_newick(text):
    """Newick -> nested dicts {'kids': [...], 'name': str|None}; branch lengths and internal labels are dropped."""
    text = text.strip().rstrip(";")
    i, n = 0, len(text)

    def skip_len():
        nonlocal i
        if i < n and text[i] == ":":
            i += 1
            while i < n and text[i] not in ",()":
                i += 1

    def node():
        nonlocal i
        if text[i] == "(":
            i += 1
            kids = [node()]
            while text[i] == ",":
                i += 1
                kids.append(node())
            i += 1
            while i < n and text[i] not in ",():;":
                i += 1
            skip_len()
            return {"kids": kids, "name": None}
        j = i
        while j < n and text[j] not in ",():":
            j += 1
        name = text[i:j].strip()
        i = j
        skip_len()
        return {"kids": [], "name": name}

    return node()


def leaves(x):
    return [x["name"]] if not x["kids"] else [l for k in x["kids"] for l in leaves(k)]


def nodes(x):
    yield x
    for k in x["kids"]:
        yield from nodes(k)


def tipset(x):
    return frozenset(leaves(x))


def path_to(root, target_tips):
    """Nodes from the root down to the node whose tip set is target_tips (inclusive), or None."""
    if tipset(root) == target_tips:
        return [root]
    for k in root["kids"]:
        if target_tips <= tipset(k):
            p = path_to(k, target_tips)
            return None if p is None else [root] + p
    return None


def depth_of(root, tips):
    """Number of edges from the root to the node with this tip set."""
    p = path_to(root, tips)
    if p is None:
        raise ValueError("node not in tree")
    return len(p) - 1


def edge_distance(root, a, b):
    """Number of edges on the path between the nodes with tip sets a and b."""
    pa, pb = path_to(root, a), path_to(root, b)
    if pa is None or pb is None:
        raise ValueError("node not in tree")
    k = 0
    while k < min(len(pa), len(pb)) and pa[k] is pb[k]:
        k += 1
    return (len(pa) - k) + (len(pb) - k)


# ---------------------------------------------------------------- intronfate tables
def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def read_branches(path):
    """dollo_branches.tsv -> list of dicts with 'branch', 'node' (int), 'tips' (frozenset), 'gains', 'losses' (int)."""
    out = []
    for r in read_tsv(path):
        out.append({"branch": r["branch"], "node": int(r["node"]), "tips": frozenset(r["tips"].split(",")),
                    "tip_list": r["tips"].split(","), "gains": int(r["gains"]), "losses": int(r["losses"])})
    return out


def branch_index(branches):
    """branch name -> tip set (names as written in dollo_sites.tsv gained_at / lost_in)."""
    return {b["branch"]: b["tips"] for b in branches}


def recon_name(tip_list):
    """The node name intronfate writes in reconstructions.tsv for an internal node (tables::reconstruction_table)."""
    names = ",".join(tip_list[:4]) + (f"…({len(tip_list)} tips)" if len(tip_list) > 4 else "")
    return "node" + names.replace(",", "-")


def read_dollo_sites(path, bidx):
    """dollo_sites.tsv -> {site: {'gains', 'losses', 'gained_at' (tipset|None), 'lost_in' [tipsets]}}."""
    out = {}
    for r in read_tsv(path):
        g = r["gained_at"]
        out[r["site"]] = {
            "gains": int(r["gains"]), "losses": int(r["losses"]),
            "gained_at": bidx[g] if g not in ("", "-") else None,
            "lost_in": [bidx[x] for x in r["lost_in"].split(";") if x and x != "-"],
        }
    return out


def read_site_losses(path):
    """dollo_sites.tsv -> {site: losses} (all that the placement comparison needs)."""
    return {r["site"]: int(r["losses"]) for r in read_tsv(path)}


def read_matrix(path):
    """matrix.tsv (format=2) -> (taxa, {site: {taxon: '1'|'0'|'?'}}). Taxa are the columns after 'pattern'."""
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
    head = lines[0].split("\t")
    k = head.index("pattern") + 1
    taxa = head[k:]
    rows = {}
    for l in lines[1:]:
        c = l.split("\t")
        rows[c[0]] = dict(zip(taxa, c[k:]))
    return taxa, rows


def unit_state(cells, members):
    """Present/absent/missing of a lineage made of several tips: '1' if any carries, '0' if all observed lack it, '?' otherwise."""
    v = [cells[m] for m in members]
    if "1" in v:
        return "1"
    if "0" in v:
        return "0"
    return "?"


def ancestral_counts(recon_path, branches):
    """reconstructions.tsv -> {method: {tipset: n_present}} for internal nodes and tips."""
    by_name = {recon_name(b["tip_list"]): b["tips"] for b in branches if len(b["tips"]) > 1}
    counts = {}
    with open(recon_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["state"] != "present":
                continue
            ts = frozenset([r["node"]]) if r["is_tip"] == "true" else by_name[r["node"]]
            m = counts.setdefault(r["method"], {})
            m[ts] = m.get(ts, 0) + 1
    return counts


# ---------------------------------------------------------------- analyses
def origin_profile(sites, rows, members, root):
    """Where the positions a lineage carries were gained, and where the ones it lost were gained.

    For every ancestor a of the lineage (from the root down to the lineage's stem), counts positions gained at a that the lineage
    carries, lacks (observed absent) or cannot be seen in; plus positions gained inside the lineage. Returns a list of dicts in
    root-to-tip order.

    `by_construction` marks the lineage's parent node and its own stem: classic Dollo places a gain at the common ancestor of the
    carriers, so a position gained at the parent (carriers on both sides of the split, one side being the lineage) or on the stem is
    carried by the lineage by definition. Only origins above the parent say anything about retention.
    """
    mset = frozenset(members)
    path = path_to(root, mset)
    if path is None:
        raise ValueError("the lineage is not a clade of the tree")
    anc = [tipset(x) for x in path]                     # root ... lineage stem
    prof = {a: {"node": a, "present": 0, "absent": 0, "missing": 0, "by_construction": i >= len(anc) - 2} for i, a in enumerate(anc)}
    inside = {"node": "inside", "present": 0, "absent": 0, "missing": 0, "by_construction": True}
    for s, d in sites.items():
        g = d["gained_at"]
        if g is None:
            continue
        st = unit_state(rows[s], members)
        key = {"1": "present", "0": "absent", "?": "missing"}[st]
        if g in prof:
            prof[g][key] += 1
        elif g < mset:
            inside[key] += 1
    return [prof[a] for a in anc] + [inside]


def retention_by_origin(sites, rows, unit_members, root):
    """For a lineage, retention (present / (present + absent)) of positions gained at each of its ancestors."""
    out = []
    for p in origin_profile(sites, rows, unit_members, root):
        obs = p["present"] + p["absent"]
        out.append(dict(p, retention=(p["present"] / obs) if obs else float("nan")))
    return out


def placement_support(losses_by_tree, carrier_sites):
    """Per placement: total losses over all sites and over the lineage's carrier sites, and how many sites each placement is (one of)
    the cheapest for. losses_by_tree: {tree_id: {site: losses}} (same sites in every tree). Returns (per_tree rows, per_site rows)."""
    trees = sorted(losses_by_tree)
    sites = sorted(losses_by_tree[trees[0]])
    per_tree = {t: {"tree": t, "losses_all": 0, "losses_carrier": 0, "best_all": 0, "best_carrier": 0, "best_absent": 0,
                    "unique_best_all": 0} for t in trees}
    per_site = []
    for s in sites:
        v = {t: losses_by_tree[t][s] for t in trees}
        lo = min(v.values())
        best = [t for t in trees if v[t] == lo]
        informative = len(best) < len(trees)
        carrier = s in carrier_sites
        for t in trees:
            per_tree[t]["losses_all"] += v[t]
            if carrier:
                per_tree[t]["losses_carrier"] += v[t]
        if informative:
            for t in best:
                per_tree[t]["best_all"] += 1
                per_tree[t]["best_carrier" if carrier else "best_absent"] += 1
            if len(best) == 1:
                per_tree[best[0]]["unique_best_all"] += 1
        per_site.append({"site": s, "carrier": carrier, "informative": informative, "min_losses": lo, "n_best": len(best),
                         "best": ",".join(best) if informative else "all"})
    rows = sorted(per_tree.values(), key=lambda r: (r["losses_all"], r["tree"]))
    for key in ("losses_all", "losses_carrier"):
        lo = min(r[key] for r in rows)
        for r in rows:
            r["delta_" + key.split("_")[1]] = r[key] - lo
    return rows, per_site


# ---------------------------------------------------------------- calibration: thinning a known lineage to Alith's intron count
def thin_matrix(src, dst, taxon, keep_prob, seed):
    """Copy matrix.tsv, turning each '1' of `taxon` into '0' with probability 1 - keep_prob(site).

    keep_prob: a number, or a function site -> probability (age-matched thinning). The pattern/n_taxa columns are recomputed for the
    changed rows; the header comment is kept. A thinned intron becomes an observed absence, as a real loss would.
    """
    rng = random.Random(seed)
    kp = keep_prob if callable(keep_prob) else (lambda s, p=keep_prob: p)
    with open(src) as f, open(dst, "w") as out:
        head = None
        for line in f:
            if line.startswith("#"):
                out.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = c
                ti = head.index(taxon)
                pi = head.index("pattern")
                ni = head.index("n_taxa")
                out.write(line)
                continue
            if c[ti] == "1" and rng.random() >= kp(c[0]):
                c[ti] = "0"
                cells = c[pi + 1:]
                c[pi] = "".join(cells)
                c[ni] = str(cells.count("1"))
            out.write("\t".join(c) + "\n")


def carrier_only(src, dst, taxon):
    """Keep only the rows where `taxon` carries the position ('1')."""
    with open(src) as f, open(dst, "w") as out:
        head = None
        for line in f:
            if line.startswith("#"):
                out.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = c
                ti = head.index(taxon)
                out.write(line)
                continue
            if c[ti] == "1":
                out.write(line)
