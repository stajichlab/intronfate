import csv, statistics as st, collections
rows = list(csv.DictReader(open("version_audit.tsv"), delimiter="\t"))
print("rows", len(rows), "errors", sum(1 for r in rows if r.get("error")))
def f(r, k):
    try: return float(r[k])
    except (KeyError, ValueError): return None
def q(v):
    v = sorted(x for x in v if x is not None)
    if not v: return "n=0"
    return f"n={len(v):3d} p10/p50/p90={v[int(.1*(len(v)-1))]:.2f}/{v[len(v)//2]:.2f}/{v[int(.9*(len(v)-1))]:.2f}"
ok = [r for r in rows if f(r, "cur_rec") is not None and f(r, "cur_ref_introns") and f(r, "cur_ref_introns") >= 200]
print("genomes with >=200 reference introns:", len(ok))
for m in ("cur_prec", "cur_rec", "cur_rec_ge200", "cur_none_frac", "cur_mrna_ratio"):
    print(f"\n{m} by version")
    for v in sorted({r["version"] for r in ok}):
        sel = [r for r in ok if r["version"] == v]
        x = [f(r, m) for r in sel if (m != "cur_rec_ge200" or (f(r, "cur_n_long") or 0) >= 20)]
        print(f"  {v:14s} {q(x)}")
print("\nby PASA training output present (has_pasa_gff3)")
for flag in ("1", "0"):
    sel = [r for r in ok if r.get("has_pasa_gff3") == flag]
    print(f" pasa={flag}: rec {q([f(r,'cur_rec') for r in sel])} | prec {q([f(r,'cur_prec') for r in sel])}")
print("\nversion x has_pasa: median recall")
for v in sorted({r["version"] for r in ok}):
    for flag in ("1", "0"):
        sel = [f(r, "cur_rec") for r in ok if r["version"] == v and r.get("has_pasa_gff3") == flag]
        if sel: print(f"  {v:14s} pasa={flag} n={len(sel):3d} median rec={st.median(sel):.2f}")
print("\nby class (n>=15): median rec / prec / rec>=200 / none_frac")
byc = collections.defaultdict(list)
for r in ok: byc[r["CLASS"]].append(r)
for c, v in sorted(byc.items(), key=lambda kv: -len(kv[1])):
    if len(v) >= 15:
        g = lambda k: st.median([f(r, k) for r in v if f(r, k) is not None])
        print(f"  {c[:22]:22s} n={len(v):3d} rec={g('cur_rec'):.2f} prec={g('cur_prec'):.2f} rec200={g('cur_rec_ge200'):.2f} none={g('cur_none_frac'):.2f}")
print("\nworst 25 genomes by recall (current): asm, version, class, rec, prec, rec>=200, none_frac, pasa")
for r in sorted(ok, key=lambda r: f(r, "cur_rec"))[:25]:
    print(f"  {r['ASMID'][:28]:28s} {r['version']:14s} {r['CLASS'][:16]:16s} rec={f(r,'cur_rec'):.2f} prec={f(r,'cur_prec'):.2f} r200={f(r,'cur_rec_ge200'):.2f} none={f(r,'cur_none_frac'):.2f} pasa={r.get('has_pasa_gff3')}")
pair = [r for r in ok if f(r, "pre_rec") is not None]
print("\nPAIRED: genomes with a pre_rc6 prediction:", len(pair))
if pair:
    print("pre_version counts:", collections.Counter(r.get("pre_version", "") for r in pair).most_common(6))
    for m in ("rec", "prec", "rec_ge200", "none_frac"):
        d = [f(r, "cur_" + m) - f(r, "pre_" + m) for r in pair if f(r, "cur_" + m) is not None and f(r, "pre_" + m) is not None]
        d.sort()
        print(f"  {m:10s} cur-pre: median {st.median(d):+.3f}  p10 {d[int(.1*(len(d)-1))]:+.3f}  p90 {d[int(.9*(len(d)-1))]:+.3f}  better {sum(x>0.02 for x in d)} worse {sum(x<-0.02 for x in d)} same {sum(abs(x)<=0.02 for x in d)}")
    print("  by current version, median recall change:")
    for v in sorted({r["version"] for r in pair}):
        d = [f(r, "cur_rec") - f(r, "pre_rec") for r in pair if r["version"] == v]
        print(f"    {v:14s} n={len(d):3d} median {st.median(d):+.3f}")
    print("\n  genomes where current is clearly worse than pre_rc6 (rec drop > 0.05):")
    for r in sorted(pair, key=lambda r: f(r, "cur_rec") - f(r, "pre_rec"))[:15]:
        print(f"    {r['ASMID'][:28]:28s} {r['version']:14s} pre={r.get('pre_version',''):14s} rec pre {f(r,'pre_rec'):.2f}->{f(r,'cur_rec'):.2f} prec {f(r,'pre_prec'):.2f}->{f(r,'cur_prec'):.2f}")
ref7 = ("ASM1142v1", "CNA3", "NC12", "R64", "ASM294v3", "MYCGR", "ASM14353v4")
print("\nthe reference genomes (current vs pre_rc6): rec / prec / rec>=200")
for r in rows:
    if any(x in r["ASMID"] for x in ref7):
        print(f"  {r['ASMID'][:26]:26s} cur {r['version']:14s} {f(r,'cur_rec') or 0:.2f}/{f(r,'cur_prec') or 0:.2f}/{f(r,'cur_rec_ge200') or 0:.2f}   pre {r.get('pre_version',''):14s} {f(r,'pre_rec') or 0:.2f}/{f(r,'pre_prec') or 0:.2f}/{f(r,'pre_rec_ge200') or 0:.2f}")
