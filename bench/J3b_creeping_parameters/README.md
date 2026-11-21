# J3b: which parameters creep in the 4-class fits, and which stopping rule is safe?

[J3](../J3_shape_ridge/README.md) showed the gamma shape is well determined and that the slowness is in the EM, but did not say which parameters are still moving at the iteration cap.
This job stops the same 4-class fit (J1 matrix, 15,697 sites; hard, and soft p0 0.8 / k 3) at 500, 1,000, 2,000, 5,000, 10,000, 20,000 and 50,000 iterations and with looser tolerances
(`--ml-tol`/`--ml-param-tol`: 1e-5/1e-6, 1e-4/1e-5, 1e-3/1e-4), and `analyse.py` (selftest in `bench/test_support_vs_reference.py`) compares each with the 50,000-iteration fit.
`run.sbatch` (15 minutes, one core per fit); `results/creep.tsv`.

## Results

Hard fit (the 50,000 cap is the converged fit: it stopped at 25,292 iterations; the default tolerances converge at ~20,000):

| stopped at | ln L below the final | gamma shape off by | gains / losses off by | largest branch change | Spearman of ln LR / top 100 kept |
|---|---|---|---|---|---|
| 500 | 0.11 | 0.033 | −0.4% / −2.5% | 0.17 | 0.985 / 100 |
| 1,000 | 0.037 | 0.018 | −0.3% / −1.5% | 0.11 | 0.992 / 100 |
| 2,000 | 0.010 | 0.009 | −0.2% / −0.8% | 0.057 | 0.998 / 100 |
| 5,000 | 0.0011 | 0.003 | −0.06% / −0.3% | 0.018 | 0.9994 / 100 |
| 10,000 | 0.0001 | 0.001 | −0.02% / −0.08% | 0.006 | 0.99994 / 100 |
| tol 1e-4 / 1e-5 (2,618 its) | 0.005 | 0.006 | −0.14% / −0.6% | 0.041 | 0.9988 / 100 |
| tol 1e-5 / 1e-6 (10,895 its) | 0.0001 | 0.001 | −0.02% / −0.07% | 0.005 | 0.99995 / 100 |
| tol 1e-3 / 1e-4 (480 its) | 0.12 | 0.033 | −0.4% / −2.6% | 0.17 | 0.985 / 100 |

- **One parameter creeps in the hard fit: the loss probability of the branch above *S. cerevisiae* + *N. crassa* + *A. nidulans*.** At every stopping point it is the largest change: 0.25 after
  500 iterations, 0.18 after 1,000, 0.094 after 5,000, 0.0819 after 10,000, against 0.0760 at the end. It trades off against gains on its daughter branches, and the data barely constrain it
  (moving it by 0.17 costs 0.11 ln L). Everything else is within 0.001 of its final value by 5,000 iterations.
- **The headline outputs are stable long before the fit converges.** From 2,000 iterations on, gains and losses are within 0.8%, the shape within 0.01, and the 100 highest ln LR sites are the
  same 100; Spearman of ln LR is ≥ 0.998 (≥ 0.9994 at 5,000).
- **Soft fit (p0 0.8, k 3): the same kind of creep, on a different branch.** The lagging parameter is the loss probability above *N. crassa* + *A. nidulans* (0.572 final; 0.70 at 2,000, 0.577 at
  10,000, 0.574 at 20,000). The soft fit is much slower and is not stable below 5,000 iterations (Spearman 0.68 at 500, 0.78 at 1,000, 0.98 at 2,000; shape off by 2.4 at 500). From 5,000 on
  it is: shape within 0.007, losses within 0.04%, Spearman ≥ 0.9999, top 100 identical. ln L is still 5.1 below the 50,000-iteration value at 5,000 iterations and 1.6 at 20,000, and the 50,000
  fit is itself unfinished (the gap about halves per doubling), but that ln L has no visible effect on the outputs.
- **The looser tolerances.** Hard: 1e-5/1e-6 stops at 10,895 iterations with outputs indistinguishable from converged; 1e-4/1e-5 at 2,618 with ln L 0.005 below, losses −0.6%, top 100 intact;
  1e-3/1e-4 stops at 480 iterations, which is too early (the creeping branch is off by 0.17). Soft: only 1e-3/1e-4 stops before 50,000 (24k iterations), and its outputs match the cap.

## What this means

- Defaults stay safe; they are slow, not wrong. `--ml-tol 1e-5 --ml-param-tol 1e-6` halves the hard 4-class fit time (25,292 → 10,895 iterations) at no visible cost. I do not recommend the
  looser settings as the default (the stop depends on the slowest parameter, which here is poorly determined). **Not changed in the code; it is the owner's decision.**
- A reported loss probability on the branch above *S. cerevisiae* + *N. crassa* + *A. nidulans* (hard) or *N. crassa* + *A. nidulans* (soft) should not be read to more than about two decimals from a capped
  fit; gains and losses totals and the origin ranking are stable from ~5,000 iterations.
- Soft fits need at least ~5,000 iterations; below that the ranking and the shape are unreliable. A convergence warning on a soft 4-class fit does not by itself mean the outputs are unsettled
  after ~20,000 iterations.
- Limit: one data set (the J1 matrix) and one soft setting; the 'final' soft fit is not converged, so the soft gaps are relative to the 50,000-iteration state.
