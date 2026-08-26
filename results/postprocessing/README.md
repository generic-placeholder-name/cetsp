# Mennell postprocessing results

## Protocol

- Release build, `CETSP_SET_BACKEND=AUTO` (resolved to Abseil).
- Fixed solver seed `123456789`; automatic repeat parallelism on 16 logical CPUs.
- 1,000 CETSP repeats per instance, except 10,000 for `bonus1000.txt`.
- One LKH ordering pass followed by one grouped-circle Gurobi/SOCP position pass.
- LKH 3 and Gurobi 12.0.2 were invoked as external executables.
- Lei–Hao values and historical CETSP values came from the repository's
  `data/mennell/reference_tours/` and `results/standalone/comparison.csv`,
  respectively.

The repository notes that its publicly available Mennell files are not
necessarily identical to the files used in every earlier paper. These numbers
therefore retain the same comparability limitations as the existing table.

## Aggregate result

| Measure | Result |
| --- | ---: |
| Instances | 62 |
| Valid final tours | 62 |
| Improved by LKH | 2 |
| Improved by SOCP after LKH | 3 |
| Improved by either postprocessor | 5 |
| Mean postprocessing improvement | 0.009103% |
| Median postprocessing improvement | 0% |
| Historical CETSP mean gap to Lei–Hao | 1.110684% |
| Current pre-postprocessing mean gap to Lei–Hao | 1.070744% |
| Final mean gap to Lei–Hao | 1.061437% |
| Current final solutions beating Lei–Hao | 1 of 62 |

Across instances, the final run was on average 0.046428% shorter than the
stored historical CETSP result, but it was shorter on 30 instances and longer
on 32. This is a stochastic comparison; the AUTO/Abseil backend does not promise
cross-process replay.

## Instances changed by postprocessing

| Instance | Before | After LKH | After SOCP | Gain | Final gap to Lei–Hao |
| --- | ---: | ---: | ---: | ---: | ---: |
| rotatingDiamonds5 | 1544.178436 | 1538.252460 | 1538.252460 | 0.383762% | 1.820451% |
| team1_100rdmRad | 397.767591 | 397.138929 | 397.138929 | 0.158047% | 2.213928% |
| dsj1000_or2 | 924.458778 | 924.458778 | 924.296547 | 0.017549% | 1.626665% |
| bubbles3 | 531.148800 | 531.148800 | 531.132716 | 0.003028% | 0.222229% |
| pcb442_or2 | 334.333597 | 334.333597 | 334.326934 | 0.001993% | 4.527471% |

`bonus1000.txt` finished at 379.133385 versus the Lei–Hao value 384.365, a
-1.361106% gap. Postprocessing did not improve that candidate. The stored
historical CETSP result is 378.6217.

## Runtime

| Stage | Total over 62 instances |
| --- | ---: |
| Parallel CETSP search | 123.128 s |
| LKH | 3.476 s |
| Gurobi/SOCP | 24.737 s |

Postprocessing added 28.213 seconds, or about 22.9% relative to the parallel
search time in this run. The detailed per-instance data is in
`mennell_postprocess.csv`.
