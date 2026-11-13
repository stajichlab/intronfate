import csv, json, os, gzip
R="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
outs={r["ASMID"]:r["out"] for r in csv.DictReader(open("../P1d_version_audit/version_audit.tsv"),delimiter="\t")}
for a in ["GCF_000011425.1_ASM1142v1","GCF_000219625.1_MYCGR_v2.0","GCF_001661235.1_Picme2","GCF_000151645.1_ASM15164v1","GCF_000143535.2_ASM14353v4","GCF_000182925.2_NC12","GCF_000002945.2_ASM294v3"]:
    o=outs[a]; d=f"{R}genome_annotation/{o}/predict_results/"
    p=[d+f for f in (o+".parameters.json",o+".parameters.json.gz") if os.path.exists(d+f)]
    print("==",a,o)
    if p:
        j=json.load(open(p[0])) if p[0].endswith(".json") else json.load(gzip.open(p[0],"rt"))
        def flat(x,pre=""):
            for k,v in (x.items() if isinstance(x,dict) else []):
                if isinstance(v,dict): yield from flat(v,pre+k+".")
                else: yield pre+k, v
        for k,v in flat(j):
            s=str(v)
            if any(t in k.lower() for t in ("intron","train","pasa","rna","augustus","genemark","snap","glimmer","weight","busco","min_","max_","source","stringtie","bam")): print("  ",k,"=",s[:110])
    st=d+o+".stats.json"
    if os.path.exists(st): 
        j=json.load(open(st)); print("   stats keys:", list(j)[:12])
