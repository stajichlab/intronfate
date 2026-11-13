import pyarrow.parquet as pq, sys
T="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
for n in ("pfam","species","samples","gene_proteins"):
    f=pq.ParquetFile(T+n+".parquet"); print("==",n,f.metadata.num_rows,"rows"); print(f.schema_arrow)
