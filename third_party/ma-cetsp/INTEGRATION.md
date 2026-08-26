# CETSP integration

This directory contains MA-CETSP from
<https://github.com/leizy1008/MA-CETSP> at commit
`d391a83871d4323e96536c92612e32d11e2faaa0`. Its MIT license is in
`LICENSE`.

The executable accepts these additional options:

- `--instance_file`: a node count followed by `x y radius` rows, with the
  depot at index zero;
- `--seed_file`: `MA_CETSP_SEEDS 1`, node and tour counts, then `TOUR`
  records containing one `id x y` row per node;
- `--result_file`: the output tour path;
- `--lkh_executable` and `--lkh_temp_root`: runtime LKH configuration;
- `--patience`: maximum generations without an incumbent improvement. Its
  default is one tenth of `--iteration`;
- `--initial_patience`: generations allowed for the first improvement. Its
  default is `--patience`;
- `--solver_threads`: Gurobi threads used by each continuous optimization.
  Zero retains Gurobi's default;
- `--snapshot_dir`: writes population seed files, member statistics, and a
  manifest for later analysis;
- `--snapshot_interval`: records the population every given number of
  generations. Improvements, initialization, and termination are always
  recorded when `--snapshot_dir` is set.

The CMake project discovers current Gurobi installations through
`GUROBI_ROOT`. QCP output points are projected into their corresponding disks
before their tour length is accepted.

Initialization output reports how many supplied seed tours remain after local
optimization and diversity filtering, how many generated tours fill the
population, and the initialization duration. The summary reports time to the
best solution and total runtime.
