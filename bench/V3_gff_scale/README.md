# V3 — genome-scale GFF performance

Question: how does the projection path (GFF3 → protein space) perform on real, large annotation — in time, memory
and scaling? V2 used nine small genomes (5.9–74 k CDS rows each); this finds the largest GFF3s in the BFD and
also scales the number of genomes.

**Inputs** (read-only, from `BFD/Fungi_BFD_runs/genome_annotation/<genome>/predict_results/`): every funannotate
`*.gff3.gz` with its `*.proteins.fa.gz`. `run.sbatch` lists their sizes once (slow on the shared filesystem; cached in
`data/gff_sizes.tsv`) and `select.py` picks, deterministically:

- the `TOP_K` (8) largest GFF3s, each run alone → throughput (MB/s of compressed input), peak RSS, and the stage
  times the tool prints with `-v` (`proteins_read`, `gff_parse_and_join`, `join_finish`,
  `write_introns_and_unmatched`);
- thread scaling on the largest (`THREADS` = 1 2 4 8 16, via `RAYON_NUM_THREADS`);
- many genomes in one process: 10, 100, 1000 genomes drawn at random with a fixed seed (`G_STEPS`) → does time or
  memory grow worse than linearly with the number of genomes?
- gzip vs plain text of the largest → what decoding costs.

```sh
mkdir -p logs && sbatch -p highclock run.sbatch     # ~ minutes to an hour; the listing pass is the slow part
python3 summarize_scale.py                           # re-tabulate any time (results/ holds timings only)
```

**What the numbers can and cannot say.** Funannotate GFF3s are modest (tens of MB compressed at most), smaller than
a mammalian Ensembl file, so the largest here is "large for fungi", not "large". The G axis is where scale shows:
1 000 genomes is ~10⁷ CDS rows. For a truly huge single file, point `--gff` at one (a human or plant GFF3) and
read the same `[time]`/`[size]` lines. Only `project` is measured; `matrix`/`analyze` scale with alignments and
taxa (see `docs/HANDOFF.md`, performance backlog). Compare runs on the same node type; run each at least twice
before trusting a difference smaller than ~10 %. Peak RSS is `/usr/bin/time` max resident set.

The tool prints stage timers only with `-v`; they go to stderr as `[time] <stage>\t<seconds>` and
`[size] <name>\t<value>`, one per line (tested in `crates/intron-cli/tests/cli.rs`).

## Baseline result (2026-11-11, job 32025162, `highclock` node hz01, `REPORT.2026-11-11.before.txt`)

Binary built before the streaming-`introns.tsv` writer and the early release of per-transcript structures.

| | |
|---|---|
| Largest GFF3 in the BFD (23 600 genomes) | 14.2 MB gz = **165 MB** uncompressed (*Phakopsora pachyrhizi*, ~108 k transcripts, ~520 k introns) |
| One genome, `project` | 2.5–2.9 s wall, ~980 MB peak RSS; parse 1.17 s (≈ 140 MB/s uncompressed), proteins 0.35 s, join 0.22 s, write 0.6 s |
| Threads 1 → 16 | **no speed-up** (1.00 → 1.02×): reading and parsing are single-threaded; only the (cheap) finalise step is parallel |
| gzip vs plain text | decoding adds 0.15 s, 13 % of the parse stage |
| 10 / 100 / 1 000 genomes in one process | 2.5 s, 0.5 GB / 27 s, 5.0 GB / **294 s, 55.8 GB** (10.9 M transcripts, 23.4 M introns, 883 MB gz) |

Reading it: time is linear in the number of genomes (≈ 0.27 s each) but **memory is the limit** — ~5 KB per transcript,
~2.4 KB per intron, so 1 000 genomes nearly fills a 64 GB request. The largest contributors, from the code:
`introns_table` made one `String` per cell (18 per intron, ~9·10⁶ allocations for the largest genome) and then copied
the whole file into one more `String`; `SiteRow` clones four strings per intron; every `TranscriptAnnotation` stayed
alive after the join. The first and last are fixed (streaming writer; transcripts dropped after the join); the
"after" run below measures them. Not addressed: `SiteRow` string clones, parallel GFF parsing, `zlib-rs`.

## After the two memory fixes (job 32027709, same inputs and node type, `REPORT.2026-11-11.after.txt`)

One replicate each (the baseline and this run, same partition), so differences under ~10 % are noise; these are not.
Output sizes and intron counts are identical before and after.

| | before | after |
|---|---|---|
| Largest genome, peak RSS | 980 MB | **402 MB** (−59 %) |
| Largest genome, wall / write stage | 2.6 s / 0.6 s | 2.0 s / 0.17 s |
| 100 genomes, RSS / wall | 5.0 GB / 27.3 s | **2.9 GB / 19.7 s** |
| 1 000 genomes, RSS / wall | 55.8 GB / 294 s | **32.5 GB / 226 s** (−42 % memory) |

Still ~3 KB per transcript at 1 000 genomes after that round. Remaining contributors (from the code, not yet measured separately):
`SiteRow` clones four strings per intron; `ProteinSet` holds every proteome (10.9 M proteins) although `project`
only needs ids and lengths after the join; the parse stage (single-threaded, 87 s of 226 s) and protein read (88 s)
now dominate time. Next measurement: parse the protein FASTA once per genome in parallel with the GFF.

## Rounds 2–4: parallel files, streaming join, leaner join state (`REPORT.2026-11-11.after2/3/4.txt`)

Changes, each checked output-identical (the golden examples, the V2 report, and equivalence tests): GFF files are
read and finalised **in parallel**, one independent annotation per file, in thread-width batches; `-p` files are read
in parallel and, for `project`, as ids and lengths only; the annotation is **streamed into the join**
(`JoiningSink`) so no transcript list is ever held; per-protein hash maps and `BTreeSet`s became dense arrays and a
sorted row index; the version-stripped id index is built only on the first inexact lookup; the isoform check and
final sort stopped hashing per comparison; unused header descriptions are no longer stored.

| 1 000 genomes in one process | wall | peak RSS |
|---|---|---|
| baseline | 294 s | 55.8 GB |
| + streaming `introns.tsv`, early transcript release | 226 s | 32.5 GB |
| + parallel files, lengths-only proteins, leaner join rows | 77 s | 22.3 GB |
| + streamed join, dense per-protein state, lazy stripped index | **52 s** | **7.0 GB** |

Other rows of the final run: largest genome 1.6 s / 222 MB; 100 genomes 3.5 s / 0.9 GB; with 100 genomes wall time
falls 12.3 s → 3.5 s from 1 to 16 threads (3.5×; parse 4.0×, proteins 6.4×); a single genome still does not scale
with threads (one file). Memory is now ~640 B per transcript and ~300 B per protein id.

What dominates now (final run, 1 000 genomes): proteome ids (3.3 GB resident after `proteins_read`, ~300 B each:
the id is stored in the sequence record and in the exact-match map) and the GFF parse (30 s of 54 s, parallel over
files only). Not done: sharing one allocation per id (`Arc<str>`), parallelism inside one large file, `zlib-rs`.
