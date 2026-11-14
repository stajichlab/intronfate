#!/usr/bin/env python3
"""P1 pilot (a): intron density per genome from the BFD tables.

Reads gene_introns and gene_CDS (streamed in batches, memory-light), counts introns
and CDS bases per genome (genome = first 8 characters of the ID, the LOCUSTAG),
joins taxonomy and assembly stats, and writes one TSV. Read-only on the BFD tables.

Density = introns per kb of CDS, counted over all transcripts (isoforms, if any, count
in both numerator and denominator). This is a screening metric to find reduced-intron
lineages, not a final estimate.
"""
import sys
import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow as pa

T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
OUT = sys.argv[1] if len(sys.argv) > 1 else "intron_density_per_genome.tsv"
BATCH = 4_000_000


def count_by_genome(path, id_col, value_expr=None, cols=None):
    """Stream a parquet file; return {LOCUSTAG: sum} (rows counted if value_expr is None)."""
    pf = pq.ParquetFile(path)
    totals = {}
    seen = 0
    for batch in pf.iter_batches(batch_size=BATCH, columns=cols):
        tbl = pa.Table.from_batches([batch])
        tag = pc.utf8_slice_codeunits(tbl[id_col], 0, 8)
        if value_expr is None:
            val = pa.array([1] * tbl.num_rows, type=pa.int64())
        else:
            val = value_expr(tbl)
        g = pa.table({"tag": tag, "v": val}).group_by("tag").aggregate([("v", "sum")])
        for t, s in zip(g["tag"].to_pylist(), g["v_sum"].to_pylist()):
            totals[t] = totals.get(t, 0) + s
        seen += tbl.num_rows
        print(f"{path.split('/')[-1]}: {seen:,} rows", file=sys.stderr, flush=True)
    return totals


introns = count_by_genome(T + "gene_introns.parquet", "intron_id", cols=["intron_id"])
cds_bp = count_by_genome(
    T + "gene_CDS.parquet",
    "transcript_id",
    value_expr=lambda t: pc.add(pc.subtract(t["end"], t["start"]), 1),
    cols=["transcript_id", "start", "end"],
)

sp = pq.read_table(T + "species.parquet").to_pylist()
asm = {r["ASMID"]: r for r in pq.read_table(T + "asm_stats.parquet").to_pylist()}

cols = ["LOCUSTAG", "ASMID", "PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS",
        "SPECIES", "n_introns", "cds_bp", "introns_per_kb_cds", "N50_bp", "contig_count",
        "annotation_status"]
n_bad = 0
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in sp:
        tag = r["LOCUSTAG"]
        ni = introns.get(tag, 0)
        bp = cds_bp.get(tag, 0)
        if bp == 0:
            continue
        dens = ni / (bp / 1000.0)
        if not (0 <= dens <= 20):  # sanity bound; flag, do not drop
            n_bad += 1
        a = asm.get(r["ASMID"], {})
        row = [tag, r["ASMID"], r["PHYLUM"], r["SUBPHYLUM"], r["CLASS"], r["ORDER"],
               r["FAMILY"], r["GENUS"], r["SPECIES"], ni, bp, f"{dens:.4f}",
               a.get("N50_bp", ""), a.get("contig_count", ""), a.get("annotation_status", "")]
        fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")
print(f"done; genomes with CDS: {sum(1 for t in cds_bp if cds_bp[t])}; "
      f"density outside 0-20/kb: {n_bad}", file=sys.stderr)
