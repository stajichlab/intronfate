# 04_public_genomes: inputs for the public-genome how-to

The genome lists and species trees for [`docs/HOWTO.md`](../../docs/HOWTO.md). Unlike examples 01-03, this one downloads complete
genomes and needs OrthoFinder and MAFFT, so it is not part of `cargo test`. It is tested by running the how-to page as written on a
cluster ([`bench/D1_howto_public_genomes`](../../bench/D1_howto_public_genomes/README.md)).

| File | What |
|---|---|
| `fungi5.tsv`, `fungi5.nwk` | Example 1: five fungi, outgroup *Cryptococcus neoformans* |
| `fungi9.tsv`, `fungi9.nwk` | Example 2: nine fungi, outgroup *Spizellomyces punctatus* (the species set of the README's worked result) |

The lists have three tab-separated columns: the label (tree tip and output name), the NCBI RefSeq assembly as `<accession>_<assembly name>`
(the folder name under `https://ftp.ncbi.nlm.nih.gov/genomes/all/`), and the organism. All are NCBI "reference genome" assemblies with
RefSeq annotation, checked on 2026-11-13.

The trees are the accepted taxonomy with no branch lengths. They are not inferred from these data, and the Dollo gain placements depend on them.
