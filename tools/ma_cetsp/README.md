# MA-CETSP adapters

`seed_export.cpp` creates MA-CETSP populations from construction-heuristic
tours. `result_check.cpp` imports an MA-CETSP tour, performs the final fixed-order
SOCP pass, and validates it against the original instance.

Enable `CETSP_BUILD_MA_CETSP_INTEGRATION` to build the exporter. The result
checker is also built when `CETSP_BUILD_POSTPROCESSING` is enabled. The complete
benchmark protocol is documented in [`experiments/ma_cetsp`](../../experiments/ma_cetsp).
