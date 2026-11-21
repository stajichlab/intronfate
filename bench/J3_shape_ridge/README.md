# J3: is the 4-class gamma fit on a flat ridge in the shape?

[J1](../J1_junction_support/README.md) saw `-m likelihood --ml-rate-categories 4` creep for tens of thousands of iterations (shape 2.1 → 4.3, ln L
rising by ~0.006 per iteration) and called it a flat ridge in the gamma shape. This job tests that directly: the shape is **fixed** at each value of a
grid (`--ml-gamma-shape`), every other parameter is fitted, and the ln L profile is read off. `run.sbatch`; table `results/profile.tsv`
(J1's combined matrix, 15,697 sites, nine taxa; hard cells, and soft cells with `--unsupported-p 0.8 --half-reads 3`; at most 20,000 iterations per fit).

## Result: the shape is well determined; what creeps is something else

| shape | hard ln L | converged | soft ln L | converged |
|---|---|---|---|---|
| 0.5 | −47,968.4 | no | −44,742.7 | no |
| 1 | −46,700.5 | no | −43,680.8 | no |
| 2 | −46,592.3 | no | −43,638.0 | no |
| 3 | −46,579.1 | no | −43,628.4 | no |
| 4 | −46,580.3 | yes (13,248) | **−43,625.4** | no |
| 5 | −46,586.7 | yes | −43,625.6 | no |
| 6 | −46,594.7 | yes | −43,626.3 | no |
| 8 | −46,610.2 | yes | −43,628.7 | no |
| 12 | −46,632.9 | yes | −43,633.8 | no |
| 20 | −46,657.5 | yes | −43,640.6 | no |
| 50 | −46,686.2 | yes | −43,649.4 | no |
| free | **−46,578.5** (shape 3.34, converged at 19,988) | yes | −43,625.8 (shape 4.33) | no |
| constant rates | −46,710.0 | yes | −43,656.6 | yes |

- **The profile is peaked, not flat.** Hard: ln L falls by 1.8 from shape 3.3 to 4, by 13.8 to 2, by 31 to 8, by 108 to 50. Soft: it falls by 3.3 from
  shape 4 to 8 and by 12.7 to 2. Taking a drop of 1.92 (half of χ²₁ at 95%) as the interval, the shape is roughly 2.7–4.1 (hard) and 3–6 (soft). It is not
  running off towards infinity (constant rates), and nothing in the profile says the shape is unidentified.
- **Rate variation is overwhelmingly supported.** Free shape against constant rates: 131 ln L units (hard) and 31 (soft) for one extra parameter.
- **The slow part is the EM, not the shape.** With the shape fixed at 3 or below, the hard fits still do not converge in 20,000 iterations, and at 4 and above they
  do (13,248 iterations at 4, 5,750 at 5, ~2,500 at 12). Strong rate variation makes the E/M iteration slow; weak variation does not. The soft fits do not
  converge at any fixed shape from 0.5 to 50 within 20,000 iterations (constant rates take 12,387).
- **How much is unfinished.** The free soft fit at 20,000 iterations (−43,625.8) is 1.6 ln L units below the 50,000-iteration fit J1 kept (−43,624.1), and the shape
  moved from about 4.28 to 4.33 meanwhile. So the unconverged part is about 2 units out of 43,600 and does not change the shape estimate by more than ~0.05.
  *Which* parameters creep was not measured here.

## What this changes

- J1's wording "flat ridge in the gamma shape" was wrong; it is slow EM convergence with a well-determined shape. J1 now says so.
- Use the free-shape fit for inference (the iteration cap costs about 2 ln L units), or fix the shape near the profile maximum and finish in far fewer iterations
  for large scans. Rankings from J1's 4-class fits (Spearman 0.03–0.80 against hard) came from unconverged fits and should be read with that in mind.
- Not done: which parameters move in the last iterations (the `dx` column of `-vv`), and whether a better acceleration or a looser `--ml-tol` would be a
  sound default. Both are questions about the optimiser, not the model.
