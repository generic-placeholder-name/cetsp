# MA-CETSP hybrid

The hybrid initializes Lei and Hao's memetic CETSP solver with tours produced
by this project. Each construction tour is expanded to one point per circle.
Expanded points remain on the original polygonal path and inside their assigned
circles, so expansion does not increase the tour length.

The exporter omits circles satisfied by the depot and larger circles that
contain another circle. Result checking and the final SOCP pass still use the
original instance.

Imported MA-CETSP points with at most a `1e-5 * max(1, radius)` boundary error
are projected into their assigned circle before SOCP. Larger misses are
rejected, and the final tour must pass the regular strict coverage check.

Build the exporter and result checker:

```powershell
cmake -S . -B build-hybrid -G Ninja -DCMAKE_BUILD_TYPE=Release `
  -DCETSP_BUILD_POSTPROCESSING=ON `
  -DCETSP_BUILD_MA_CETSP_INTEGRATION=ON
cmake --build build-hybrid -j 8
```

Build the vendored MA-CETSP executable separately with MSVC because the Windows
Gurobi C++ library uses the MSVC ABI:

```powershell
cmake -S third_party/ma-cetsp -B build-ma-cetsp `
  -G "Visual Studio 17 2022" -A x64 `
  -DGUROBI_ROOT=C:/gurobi1202/win64
cmake --build build-ma-cetsp --config Release -j 8
```

Run the Mennell instances with `scripts/run_mennell.ps1`. The script writes one row
after every completed instance and supports `-Resume`.

```powershell
experiments/ma_cetsp/scripts/run_mennell.ps1 `
  -SeedExporter build-hybrid/CETSP_ma_cetsp_seed_export.exe `
  -MaCetspExecutable build-ma-cetsp/Release/MA-CETSP.exe `
  -LkhExecutable C:/path/to/LKH.exe `
  -GurobiExecutable C:/gurobi/bin/gurobi_cl.exe `
  -ResultChecker build-hybrid/CETSP_ma_cetsp_result_check.exe
```

`-Instances rotatingDiamonds5,bonus1000` restricts a run to named instances.

Use `scripts/run_mennell_protocol.ps1` for the repeated benchmark protocol. It runs 20
independent populations per instance, permits 1,500 generations for the first
seeded improvement, restores the 500-generation stagnation limit afterward,
and retains the best validated tour. Runs execute concurrently and are
checkpointed individually.

```powershell
experiments/ma_cetsp/scripts/run_mennell_protocol.ps1 `
  -SeedExporter build-hybrid/CETSP_ma_cetsp_seed_export.exe `
  -MaCetspExecutable build-ma-cetsp/Release/MA-CETSP.exe `
  -LkhExecutable C:/path/to/LKH.exe `
  -GurobiExecutable C:/gurobi/bin/gurobi_cl.exe `
  -ResultChecker build-hybrid/CETSP_ma_cetsp_result_check.exe
```

The defaults match Lei and Hao's population size, generation limit, fitness
weight, edit-distance threshold, neighbor count, and best-of-20 reporting.
`-Resume` continues from the per-run metrics stored in the work directory.
`-ParallelJobs`, `-ConstructionThreads`, and `-SolverThreads` control how the
available cores are divided among independent runs.

Plot a seed population with the repository's Python environment:

```powershell
.venv/Scripts/python.exe experiments/ma_cetsp/analysis/plot_seed_tours.py `
  path/to/instance.txt path/to/population.seeds path/to/tours.png
```

MA-CETSP can record the complete population periodically and whenever it finds
a new best tour:

```powershell
build-ma-cetsp/Release/MA-CETSP.exe `
  --snapshot_dir path/to/snapshots `
  --snapshot_interval 100
```

Each snapshot has a reusable `.seeds` population file, a per-member CSV, and a
row in `manifest.csv`. Plot one or more histories with:

```powershell
.venv/Scripts/python.exe experiments/ma_cetsp/analysis/plot_population_history.py `
  path/to/history.png `
  "Seeded=path/to/seeded/manifest.csv" `
  "Native=path/to/native/manifest.csv"
```

`--patience N` sets the no-improvement stopping window independently of the
generation limit. Without it, the window is one tenth of `--iteration`.
`--initial_patience N` can grant a longer window only until the first
improvement; it then returns to `--patience`.

Compare how populations place their assigned points inside circles with:

```powershell
.venv/Scripts/python.exe experiments/ma_cetsp/analysis/plot_point_positions.py `
  path/to/instance.txt path/to/positions.png `
  "Seeded=path/to/seeded.seeds" `
  "Native=path/to/native.seeds"
```

The checked-in benchmark tables, diagnostic figures, and best validated tours
are in [`results/mennell_hybrid`](../../results/mennell_hybrid).
