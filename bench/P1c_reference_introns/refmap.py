import pyarrow.parquet as pq, os, re
T="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
S="/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"
dirs=os.listdir(S)
want={"GCF_000011425":"A.nidulans","GCF_000182965":"C.albicans","GCF_000149245":"C.neoformans H99","GCF_000182925":"N.crassa OR74A","GCF_000146045":"S.cerevisiae S288C","GCF_000002945":"S.pombe 972h-","GCF_000219625":"Z.tritici IPO323"}
sp=pq.read_table(T+"species.parquet").to_pylist()
print("species cols:", list(sp[0].keys()))
for acc,name in want.items():
    num=acc.split("_")[1]
    src=[d for d in dirs if num in d]
    bfd=[(r["LOCUSTAG"],r["ASMID"],r["SPECIES"]) for r in sp if num in str(r["ASMID"]) or num in str(r.get("ACCESSION",""))]
    print(name, "| source dirs:", src, "| BFD:", bfd[:3])
