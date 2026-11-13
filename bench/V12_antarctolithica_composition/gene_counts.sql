-- usage: duckdb < gene_counts.sql  (run in this folder; reads ../V10_antarctolithica_domains/selected_genomes.tsv; read-only on BFD tables)
CREATE TEMP TABLE sel AS SELECT locustag FROM read_csv('../V10_antarctolithica_domains/selected_genomes.tsv', delim='\t', header=true);
COPY (SELECT LOCUSTAG AS locustag, count(*) AS n_genes
      FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/gene_info.parquet')
      WHERE gene_type = 'protein_coding' AND LOCUSTAG IN (SELECT locustag FROM sel) GROUP BY 1) TO 'gene_counts.tsv' (DELIMITER '\t', HEADER);
