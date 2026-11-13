-- usage: duckdb < pfam_counts.sql   (run from this folder; reads selected_genomes.tsv; read-only on BFD tables)
CREATE TEMP TABLE sel AS SELECT locustag FROM read_csv('selected_genomes.tsv', delim='\t', header=true);
COPY (
  SELECT split_part(protein_id,'_',1) AS locustag, pfam_acc, count(*) AS n_domains, count(DISTINCT protein_id) AS n_proteins
  FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/pfam.parquet')
  WHERE domain_i_evalue < 1e-5 AND split_part(protein_id,'_',1) IN (SELECT locustag FROM sel)
  GROUP BY 1,2
) TO 'pfam_counts.tsv' (DELIMITER '\t', HEADER);
