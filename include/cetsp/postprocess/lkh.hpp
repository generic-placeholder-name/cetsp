#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

#include <cetsp/types.hpp>

namespace cetsp::postprocess {

struct LkhOptions {
    std::filesystem::path executable;
    std::size_t runs = 1;
    std::optional<std::uint64_t> seed;
};

// Reorders the existing feasible visit points with LKH. The original tour is
// returned unless LKH's permutation is shorter in the solver's true geometry.
std::vector<Point> improveOrderWithLkh(
    const std::vector<Point>& tour,
    const LkhOptions& options);

} // namespace cetsp::postprocess
