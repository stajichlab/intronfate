import pyarrow.parquet as pq
T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
for n in ["gene_proteins","gene_transcripts","gene_info"]:
    f = pq.ParquetFile(T+n+".parquet"); print(n, f.metadata.num_rows, [(x.name,str(x.type)) for x in f.schema_arrow])
    b = next(f.iter_batches(batch_size=2)); print(b.to_pylist())
