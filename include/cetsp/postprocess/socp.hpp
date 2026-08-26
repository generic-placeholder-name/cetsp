#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

#include <cetsp/types.hpp>

namespace cetsp::postprocess {

struct SocpOptions {
    std::filesystem::path gurobiExecutable;
    std::size_t maxThreads = 0;
    double timeLimitSeconds = 0.0;
};

// Assigns every circle to a currently covering visit, then globally optimizes
// visit locations for the fixed cyclic order. Visits assigned multiple circles
// are constrained to their intersection. The original tour is returned unless
// the resulting tour is valid and shorter.
std::vector<Point> improvePointsWithSocp(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles,
    const SocpOptions& options);

} // namespace cetsp::postprocess
