import pyarrow.parquet as pq, pandas as pd, glob, collections
T="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
f=pq.ParquetFile(T+"pfam.parquet"); print("pfam row groups",f.num_row_groups)
print(f.read_row_group(0,columns=["protein_id","pfam_acc","domain_i_evalue"]).slice(0,3).to_pandas())
s=pq.read_table(T+"samples.parquet").to_pandas()
print(s[s.SPECIES.str.contains("Antarctolithica",na=False)][["ASMID","STRAIN","LOCUSTAG","PHYLUM","SUBPHYLUM","CLASS"]])
print("locustag unique:", s.LOCUSTAG.nunique(), "of", len(s))
print(s.groupby(["PHYLUM"]).size().sort_values(ascending=False).head(12))
b=pq.ParquetFile(T+"busco_genome.parquet"); print(b.schema_arrow)
ids=[l[1:].split()[0] for l in open(glob.glob("../V9_antarctolithica_expanded/proteomes/*Alith*")[0]) if l[0]==">"][:3] if glob.glob("../V9_antarctolithica_expanded/proteomes/*Alith*") else None
print("V9 proteome files", sorted(glob.glob("../V9_antarctolithica_expanded/proteomes/*"))[:5], ids)
