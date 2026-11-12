# Example 3 — MALIN tables

`tiny.table` is a MALIN intron table: four taxa, five intron sites over two alignments (`K1`, `K2`), with `#MAP`
lines placing each site and `?` for cells that could not be observed. `tiny.nwk` is `((A,B),(C,D));`.

```
A  3103?          codes: 0 no intron · ? not observed · 3 / 1 / 2 intron of phase 0 / 1 / 2
B  31201          #MAP site alignment alipos     alipos = 3·column + phase
C  302?0          0 K1 3   -> column 1, phase 0
D  0000?          1 K1 7   -> column 2, phase 1     2 K1 14 -> column 4, phase 2
                  3 K2 9   -> column 3, phase 0     4 K2 10 -> column 3, phase 1
```

`run.sh` imports the table (`convert --from-malin`), scores each alignment, exports the matrices back
(`convert --to-malin`) and checks that the data rows are identical: **round trip: identical**.

## Expected, by hand

`K1`, sites `1:0` (A,B,C), `2:1` (A,B), `4:2` (B,C):

- `1:0`: carriers' ancestor is the root, D is observed absent → 1 loss.
- `2:1`: A and B are sisters → compatible, 0 losses.
- `4:2`: B and C are not sisters; A and D are absent → 2 losses.

Total **3 gains, 3 losses**. In `K2` every cell that would be a loss is `?` (C and D are unobserved where the
only carriers are) → **2 gains, 0 losses**.
