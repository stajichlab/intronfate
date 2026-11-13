import csv, collections, statistics as st
NAMES={"GCF_000011425.1_ASM1142v1":"A.nidulans","GCF_000219625.1_MYCGR_v2.0":"Z.tritici","GCF_000146045.2_R64":"S.cerevisiae","GCF_001661235.1_Picme2":"Picme2(Pichio)","GCF_000151645.1_ASM15164v1":"GCF_0001516(Pezizo)","GCF_000143535.2_ASM14353v4":"B.cinerea","GCF_000182925.2_NC12":"N.crassa","GCF_000149245.1_CNA3":"C.neoformans","GCF_000002945.2_ASM294v3":"S.pombe"}
G=list(csv.DictReader(open("genes.tsv"),delimiter="\t")); E=list(csv.DictReader(open("events.tsv"),delimiter="\t"))
CL=["identical","ends_only","shifted_boundary","exon_difference","retained_intron","end_difference","different_structure","split","merged","opposite_strand_only","missing_gene","other"]
print("GENE-LEVEL CLASS (% of reference genes); rows genome/annotation")
print("%-26s"%"genome ann n"+"".join("%8s"%c[:7] for c in CL))
by=collections.defaultdict(list)
for r in G: by[(r["asmid"],r["annotation"])].append(r)
for (a,ann),v in sorted(by.items(), key=lambda kv:(NAMES[kv[0][0]],kv[0][1])):
    c=collections.Counter(r["class"] for r in v); n=len(v)
    print("%-26s"%f"{NAMES[a][:15]:15s} {ann} {n:5d}"+"".join("%8.1f"%(100*c[k]/n) for k in CL))
print("\nAMINO ACIDS for genes with an overlapping BFD gene and not identical: median missing aa / extra aa; % of ref aa missing (mean over all matched genes); genes losing >=10% of aa")
for (a,ann),v in sorted(by.items(), key=lambda kv:(NAMES[kv[0][0]],kv[0][1])):
    m=[r for r in v if r["class"] not in ("missing_gene","opposite_strand_only")]
    ni=[r for r in m if r["class"]!="identical"]
    if not m: continue
    fr=[int(r["missing_aa"])/max(int(r["ref_aa"]),1) for r in m]
    print("  %-26s matched=%5d nonident=%5d  med missing aa=%5.1f extra aa=%5.1f | mean %% ref aa missing=%5.2f | >=10%% lost: %4.1f%% | >=50%% lost: %4.1f%%"%(f"{NAMES[a][:15]} {ann}",len(m),len(ni),st.median([int(r['missing_aa']) for r in ni]) if ni else 0,st.median([int(r['extra_aa']) for r in ni]) if ni else 0,100*st.fmean(fr),100*sum(x>=.1 for x in fr)/len(fr),100*sum(x>=.5 for x in fr)/len(fr)))
print("\nBY CLASS: median missing aa / extra aa (current only, pooled over genomes)")
cur=[r for r in G if r["annotation"]=="cur"]
for k in CL:
    v=[r for r in cur if r["class"]==k]
    if v: print("  %-22s n=%6d  missing aa median=%5.1f  extra aa median=%5.1f  mean ref aa=%5.0f"%(k,len(v),st.median([int(r['missing_aa']) for r in v]),st.median([int(r['extra_aa']) for r in v]),st.fmean([int(r['ref_aa']) for r in v])))
print("\nINTRON EVENTS in genes NOT identical (current): count and %")
ev=[r for r in E if r["annotation"]=="cur"]
c=collections.Counter(r["event"] for r in ev); n=sum(c.values())
for k,x in c.most_common(): print("  %-26s %7d %5.1f%%"%(k,x,100*x/n))
print("\nREFERENCE introns >=200 bp absent from BFD: what happened (current), by genome")
for a in NAMES:
    v=[r for r in ev if r["asmid"]==a and r["event"] in ("retained","shifted","intron_differs","outside_span") and int(r["intron_len"])>=200]
    s=[r for r in ev if r["asmid"]==a and r["event"] in ("retained","shifted","intron_differs","outside_span") and int(r["intron_len"])<200]
    c=collections.Counter(r["event"] for r in v); c2=collections.Counter(r["event"] for r in s)
    f=lambda c: " ".join("%s=%d"%(k[:9],c[k]) for k in ("retained","shifted","intron_differs","outside_span"))
    print("  %-20s >=200bp: %-50s | <200bp: %s"%(NAMES[a][:20],f(c),f(c2)))
print("\nclass of genes with a missed >=200 bp ref intron (current, all genomes pooled):",collections.Counter(r["class"] for r in ev if r["event"] in ("retained","shifted","intron_differs","outside_span") and int(r["intron_len"])>=200).most_common(6))
print("\nPRE vs CUR event mix (all genomes pooled): reference introns missing")
for ann in ("pre","cur"):
    c=collections.Counter(r["event"] for r in E if r["annotation"]==ann and r["event"] in ("retained","shifted","intron_differs","outside_span"))
    n=sum(c.values()); print("  %s: "%ann+"  ".join("%s %.0f%% (%d)"%(k,100*x/n,x) for k,x in c.most_common()))
