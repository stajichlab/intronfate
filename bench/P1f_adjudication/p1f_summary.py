import csv, collections
N={"GCF_000011425.1_ASM1142v1":"A.nidulans","GCF_000219625.1_MYCGR_v2.0":"Z.tritici","GCF_000146045.2_R64":"S.cerevisiae","GCF_001661235.1_Picme2":"Picme2","GCF_000143535.2_ASM14353v4":"B.cinerea","GCF_000182925.2_NC12":"N.crassa","GCF_000149245.1_CNA3":"C.neoformans","GCF_000002945.2_ASM294v3":"S.pombe"}
I=list(csv.DictReader(open("P1f_intron_support.tsv"),delimiter="\t")); V=list(csv.DictReader(open("P1f_gene_verdicts.tsv"),delimiter="\t"))
print("INTRON LEVEL (introns lying inside a Trinity transcript span). % with exact transcript junction / conflicting junction / none")
print("%-13s %-9s %-6s %6s %7s %9s %6s"%("genome","set","len","n","exact%","conflict%","none%"))
for a,nm in N.items():
    for s in ("shared","bfd_only","ref_only"):
        for lc in ("<200",">=200"):
            r=[x for x in I if x["asmid"]==a and x["set"]==s and x["length_class"]==lc]
            if r and int(r[0]["n"])>=10:
                n=int(r[0]["n"]); print("%-13s %-9s %-6s %6d %6.0f%% %8.0f%% %5.0f%%"%(nm,s,lc,n,100*int(r[0]["ta_exact"])/n,100*int(r[0]["ta_conflict"])/n,100*int(r[0]["ta_none"])/n))
print("\nGENE LEVEL (matched non-identical genes, classes other than ends_only/missing): verdict counts per genome")
print("%-13s %9s %9s %7s %12s | among ref_right: PASA has ref (EVM overrode) / PASA lacks ref | among bfd_right: PASA has BFD / lacks"%("genome","bfd_right","ref_right","mixed","no_junction"))
tot=collections.Counter()
for a,nm in N.items():
    c=collections.Counter()
    for r in V:
        if r["asmid"]==a: c[r["verdict"]]+=int(r["n"]); c[r["pasa_status"] or "-"]+=int(r["n"]) if r["pasa_status"] else 0
    tot.update(c)
    print("%-13s %9d %9d %7d %12d | %d / %d | %d / %d"%(nm,c["bfd_right"],c["ref_right"],c["mixed"],c["no_junction_support"],c["PASA_has_ref"],c["PASA_lacks_ref"],c["PASA_has_bfd"],c["PASA_lacks_bfd"]))
print("%-13s %9d %9d %7d %12d | %d / %d | %d / %d"%("ALL",tot["bfd_right"],tot["ref_right"],tot["mixed"],tot["no_junction_support"],tot["PASA_has_ref"],tot["PASA_lacks_ref"],tot["PASA_has_bfd"],tot["PASA_lacks_bfd"]))
print("\nref_right and bfd_right by class (all genomes pooled):")
cl=collections.defaultdict(collections.Counter)
for r in V: cl[r["class"]][r["verdict"]]+=int(r["n"])
for k,c in sorted(cl.items(), key=lambda kv:-sum(kv[1].values())):
    print("  %-20s bfd_right=%5d ref_right=%4d mixed=%3d no_junction=%6d"%(k,c["bfd_right"],c["ref_right"],c["mixed"],c["no_junction_support"]))
print("\nref_right by class and PASA status (pooled): where the correct (RNA-supported) reference structure was lost")
for k,c in sorted(cl.items()):
    pass
d=collections.defaultdict(collections.Counter)
for r in V:
    if r["verdict"]=="ref_right": d[r["class"]][r["pasa_status"]]+=int(r["n"])
for k,c in d.items(): print("  %-20s PASA_has_ref=%3d PASA_lacks_ref=%3d"%(k,c["PASA_has_ref"],c["PASA_lacks_ref"]))
