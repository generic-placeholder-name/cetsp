# Mennell postprocessing comparison

`run_mennell.ps1` runs the current parallel CETSP heuristic, then LKH ordering
and grouped-circle Gurobi/SOCP point optimization for every locally comparable
Mennell instance. It records the current pre/post distances alongside the
historical CETSP result in `results/standalone/comparison.csv` and the Lei–Hao
value shipped in `data/mennell/reference_tours/`.

The script uses 1,000 randomized repeats per instance and 10,000 for
`bonus1000.txt` by default, matching the historical comparison protocol.

The detailed measurements are in
[`results/postprocessing`](../../results/postprocessing).
