#!/usr/bin/env python3
"""Tabulate a V3 run (results/*) into the performance questions: throughput, memory, stages, scaling."""
import glob, os

def kv(path):
    d = {}
    if os.path.exists(path):
        for l in open(path):
            if "\t" in l: k, v = l.rstrip("\n").split("\t", 1); d[k] = v
    return d

def mem(path):
    """stage -> (rss_mb, peak_mb) from the tool's [mem] lines (Linux only)."""
    out = {}
    for l in open(path):
        if l.startswith("[mem] "):
            k, rss, hwm = l[6:].rstrip("\n").split("\t"); out[k] = (int(rss) / 1024, int(hwm) / 1024)
    return out

def stages(path):
    out = {}
    for l in open(path):
        if l.startswith("[time] "):
            k, v = l[7:].rstrip("\n").split("\t"); out[k] = float(v)
    return out

runs = {}
for d in sorted(glob.glob("results/*/")):
    tag = os.path.basename(d.rstrip("/"))
    t = kv(d + "time.tsv"); m = kv(d + "meta.tsv"); s = kv(d + "summary_project.tsv")
    if not t: continue
    runs[tag] = dict(t=t, m=m, s=s, st=stages(d + "stderr.txt"), mem=mem(d + "stderr.txt"))

def row(tag):
    r = runs[tag]; w = float(r["t"]["wall_s"]); rss = int(r["t"]["max_rss_kb"]) / 1024
    gz = int(r["m"].get("gff_bytes_gz", 0)) / 1e6
    return w, rss, gz, r

print("== largest GFF3s, one genome per process (default threads)")
print(f"  {'genome':44s} {'gz MB':>7s} {'wall s':>8s} {'MB/s(gz)':>9s} {'RSS MB':>8s} {'parse s':>8s} {'join s':>7s} {'write s':>8s} {'transcripts':>12s} {'introns':>9s}")
for tag in runs:
    if not tag.startswith("top_"): continue
    w, rss, gz, r = row(tag)
    st = r["st"]
    print(f"  {tag[4:]:44s} {gz:7.1f} {w:8.2f} {gz/w:9.1f} {rss:8.0f} {st.get('gff_parse_and_join',0):8.2f} {st.get('join_finish',0):7.2f} {st.get('write_introns_and_unmatched',0):8.2f} {r['s'].get('transcripts_read','?'):>12s} {r['s'].get('intron_records','?'):>9s}")

print("\n== thread scaling on the largest GFF3 (wall s; speedup vs 1 thread)")
base = None
for tag in sorted((t for t in runs if t.startswith("threads_")), key=lambda t: int(t.split("_")[1])):
    w, rss, gz, r = row(tag)
    base = base or w
    print(f"  threads {tag.split('_')[1]:>3s}: wall {w:7.2f}  speedup {base/w:5.2f}  RSS {rss:7.0f} MB  parse {r['st'].get('gff_parse_and_join',0):6.2f}  user {float(r['t']['user_s']):7.2f}")

print("\n== many genomes in one process (project)")
print(f"  {'genomes':>8s} {'gz MB':>8s} {'wall s':>8s} {'RSS MB':>8s} {'parse s':>8s} {'proteins s':>10s} {'join s':>7s} {'write s':>8s} {'transcripts':>12s} {'introns':>10s} {'out MB':>8s}")
for tag in sorted((t for t in runs if t.startswith("genomes_")), key=lambda t: int(t.split("_")[1])):
    w, rss, gz, r = row(tag); st = r["st"]
    print(f"  {r['m']['genomes']:>8s} {gz:8.0f} {w:8.1f} {rss:8.0f} {st.get('gff_parse_and_join',0):8.1f} {st.get('proteins_read',0):10.1f} {st.get('join_finish',0):7.1f} {st.get('write_introns_and_unmatched',0):8.1f} {r['s'].get('transcripts_read','?'):>12s} {r['s'].get('intron_records','?'):>10s} {int(r['m']['output_bytes'])/1e6:8.0f}")

print("\n== thread scaling with 100 genomes in one process (wall s; speedup vs 1 thread)")
base = None
for tag in sorted((t for t in runs if t.startswith("gthreads_")), key=lambda t: int(t.split("_")[1])):
    w, rss, gz, r = row(tag)
    base = base or w
    print(f"  threads {tag.split('_')[1]:>3s}: wall {w:7.2f}  speedup {base/w:5.2f}  RSS {rss:7.0f} MB  parse {r['st'].get('gff_parse_and_join',0):6.2f}  proteins {r['st'].get('proteins_read',0):6.2f}")

for tag in ("gzip_largest", "plain_largest"):
    if tag in runs:
        w, rss, gz, r = row(tag)
        print(f"\n{tag}: wall {w:.2f} s  parse {r['st'].get('gff_parse_and_join', 0):.2f} s  RSS {rss:.0f} MB")
if "gzip_largest" in runs and "plain_largest" in runs:
    g = runs["gzip_largest"]["st"].get("gff_parse_and_join", 0); p = runs["plain_largest"]["st"].get("gff_parse_and_join", 0)
    print(f"  gzip decoding adds {g-p:.2f} s to the parse stage ({(g-p)/g*100:.0f}% of it)")

print("\n== resident memory after each stage (MB: now / peak so far), per run")
order = ["proteins_read", "gff_parse_and_join", "join_finish", "write_introns_and_unmatched"]
for tag in sorted(runs, key=lambda t: (t.split("_")[0], int(t.split("_")[1]) if t.split("_")[1].isdigit() else 0, t)):
    mm = runs[tag]["mem"]
    if not mm: continue
    cells = "  ".join(f"{k.split('_')[0]}={mm[k][0]:.0f}/{mm[k][1]:.0f}" for k in order if k in mm)
    print(f"  {tag[:46]:46s} {cells}")
