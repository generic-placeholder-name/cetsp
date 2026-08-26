# Optional tour postprocessing

`cetsp::postprocess` applies optional improvements to a valid CETSP tour. LKH
and Gurobi are supplied as executable paths at runtime.

- LKH keeps the visit points fixed and improves their cyclic order.
- SOCP assigns every input circle to a covering visit and uses Gurobi to
  optimize all visit locations for the fixed order. A visit assigned multiple
  circles is constrained to their intersection.
- The pipeline can alternate both stages. Every stage measures the actual tour,
  verifies coverage where points can move, and rejects a worse result.

Enable the separate `cetsp::postprocess` target with:

```sh
cmake -S . -B build -DCETSP_BUILD_POSTPROCESSING=ON
cmake --build build --target CETSP_postprocess_benchmark
```

Supply executable paths to `CETSP_postprocess_benchmark` or through the C++
options in `<cetsp/postprocess.hpp>`.

The benchmark syntax is:

```text
CETSP_postprocess_benchmark <instance> [repetitions] [seed] [max-threads]
    [lkh-executable|-] [gurobi-executable|-] [rounds]
```
