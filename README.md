# CETSP

CETSP is a randomized, expected-$O(n\log n)$ construction heuristic for the
close-enough traveling salesman problem. It is designed for fast standalone
solutions and for generating initial populations for higher-cost optimization.

## How it works

The solver uses a heuristic based on the [pair-center algorithm](https://www.sciencedirect.com/science/article/pii/S1877750324002175). It repeatedly merges two circles into a representative circle, then reconstructs a tour by expanding the hierarchy back to the original targets.

During reconstruction, selected visits are removed and their circles reinserted,
and visit positions are locally optimized against neighboring tour points.

## Requirements

You must have Boost installed. `Point` is currently a Boost.Geometry point in
the supported API, so Boost remains a public dependency. Abseil is an optional,
private implementation dependency used for the fastest set backend; when it is
unavailable, the default build falls back to `std::unordered_set`.

## Library interface

The canonical library include is:

```cpp
#include <cetsp/cetsp.hpp>
```

The supported headers live in `include/cetsp/`. They define the input geometry,
solver options, solver entry points, and tour verification helpers. In
particular, `Point` remains a Boost.Geometry type for source compatibility.

Everything under `src/` is solver implementation: merge-tree state, tour
ownership, spatial-index payloads, reconstruction, and internal geometry
helpers. Applications should include headers from `include/cetsp/`, not
headers from `src/`.

## Project layout

```text
include/cetsp/  supported library headers
src/            private implementation, including optional postprocessors
app/            command-line application
tests/          regression tests and explicit quality benchmarks
tools/          dataset utilities and MA-CETSP integration executables
data/           Mennell instances and reference tours
experiments/    reproducible benchmark and analysis scripts
experimental/   exploratory algorithms outside the supported solver
results/        checked-in measurements, figures, and best tours
paper/          manuscript source, figures, and compiled PDF
third_party/    vendored external solver sources and licenses
```

Configure `CETSP_BUILD_POSTPROCESSING=ON` to build the optional
`cetsp::postprocess` target. LKH optimizes visit order, while Gurobi/SOCP
optimizes visit positions for a fixed order. `cetsp::core` does not require
either external solver.

The optional MA-CETSP hybrid adapter and its benchmark runner are documented in
[`experiments/ma_cetsp`](./experiments/ma_cetsp). The vendored MA-CETSP solver
is built separately because its Gurobi C++ dependency must use a compatible
compiler ABI.

## Building (Windows, MinGW; UNIX/POSIX similar)

1. Open a terminal (PowerShell, Command Prompt, or UNIX shell).
2. Make a build directory (if it doesn’t exist) and navigate into it:

```bash
mkdir build   # Windows: mkdir build
cd build
```

3. Configure the project with CMake using **MinGW**:

```bash
cmake .. -G "MinGW Makefiles"
```

4. Build the project:

```bash
cmake --build .
```

The build also creates a reusable `cetsp_core` library. When testing is enabled
(the default), run the focused regression suite with:

```bash
ctest --test-dir build --output-on-failure
```

5. Go back to the base directory:

```bash
cd ..
```

**Notes:**

* On UNIX/POSIX, you can replace `MinGW Makefiles` with `Unix Makefiles` or use `Ninja` (`-G "Ninja"`).
* You can also pass options to CMake:

```bash
cmake .. -DCETSP_SET_BACKEND=AUTO -DDEBUG=OFF -G "MinGW Makefiles"
```

For an Abseil-free Release build, explicitly select the standard unordered
backend:

```bash
cmake .. -DCETSP_SET_BACKEND=UNORDERED -DDEBUG=OFF -G "MinGW Makefiles"
```

For fixed-seed replay, configure a separate ordered build:

```bash
cmake .. -DCETSP_SET_BACKEND=ORDERED -DDEBUG=OFF -G "MinGW Makefiles"
```

`CETSP_SET_BACKEND` accepts:

* `AUTO` (default) → uses Abseil `flat_hash_set` when available and otherwise falls back to `std::unordered_set`.
* `ABSEIL` → requires and uses Abseil `flat_hash_set`.
* `UNORDERED` → uses `std::unordered_set` without Abseil; iteration order is not guaranteed reproducible.
* `ORDERED` → uses `std::set` and enables the fixed-seed replay regression test.
* `CETSP_BUILD_POSTPROCESSING=ON` → builds the optional `cetsp::postprocess`
  library and `CETSP_postprocess_benchmark`; LKH and Gurobi are invoked through
  executable paths supplied at runtime.
* `CETSP_BUILD_MA_CETSP_INTEGRATION=ON` → builds the MA-CETSP seed exporter and,
  when postprocessing is enabled, its result validator.
* `DEBUG=ON` → enables debug logging (`DBG()` macros) (`OFF` by default). (Warning: the amount of debug data is massive. This is mostly used to help debug code when it's not working.)
* Copy `settings.example.txt` to `settings.txt` and add `seed <unsigned integer>` to select the solver's random stream. Exact replay between processes also requires `CETSP_SET_BACKEND=ORDERED`.
* Add `maxThreads <count>` to cap repeat workers. The default, `0`, uses `std::thread::hardware_concurrency()` (falling back to one thread), and the solver always caps the count at `numRepeats`.

## Running

* Place your local `settings.txt` file in the **base directory** (same level as the `build` folder). Generated tours are written under `output/` by default.
* Run the executable:

```bash
./build/CETSP_project.exe
```

* On Windows, you can also double-click the `.exe`.

For an isolated solution-quality check, run:

```bash
./build/CETSP_benchmark.exe ./data/mennell/instances/bonus1000.txt 10000 [seed] [max-threads]
```

The second argument is the number of randomized repetitions, followed by the
optional seed and maximum thread count. A maximum thread count of zero (the
default) selects hardware concurrency. Quality comparisons should use a Release build
because the historical `bonus1000` result used 10,000 repetitions. With
`CETSP_SET_BACKEND=ORDERED`, the optional seed replays a run using the same
executable, standard-library implementation, and target CPU, independently of
the worker count and scheduling.
It is not a portable cross-toolchain guarantee because floating-point behavior
and standard-library random distributions may differ. With Abseil or the
unordered backend, the seed still controls the solver RNG, but container
iteration can produce a different tour. When omitted, the runner generates and
prints a seed.

## Results

The construction solver is optimized for throughput. The optional MA-CETSP
hybrid trades additional runtime for better tour orderings; its 62-instance
Mennell comparison is in
[`results/mennell_hybrid`](./results/mennell_hybrid). The best validated hybrid
tours are stored in [`results/mennell_hybrid/tours`](./results/mennell_hybrid/tours).

## Limitations

The construction tour order can remain locally suboptimal, particularly on
structured instances such as the larger `bubbles` cases. The optional hybrid
uses LKH, continuous point optimization, and MA-CETSP when solution quality is
more important than construction time.
