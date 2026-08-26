# Results

- [`mennell_hybrid`](./mennell_hybrid) contains the 20-run seeded MA-CETSP
  comparison, statistical analysis, figures, and 62 validated best tours.
- [`postprocessing`](./postprocessing) contains the LKH/SOCP-only comparison.
- [`standalone`](./standalone) contains the construction heuristic comparison.
- [`runtime`](./runtime) contains the raw logs used for scaling measurements.

Scripts that regenerate these artifacts are under [`experiments`](../experiments)
and [`tools`](../tools). Large per-run work directories remain under ignored
`build-*` or `output/` paths.
