#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include <cetsp/postprocess/lkh.hpp>
#include <cetsp/postprocess/socp.hpp>
#include <cetsp/types.hpp>

namespace cetsp::postprocess {

struct PostprocessOptions {
    std::optional<LkhOptions> lkh;
    std::optional<SocpOptions> socp;
    std::size_t rounds = 1;
};

// Alternates enabled order and point-location improvements. Each backend is
// non-worsening, and the pipeline stops early when a complete round makes no
// measurable improvement.
std::vector<Point> polishTour(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles,
    const PostprocessOptions& options);

} // namespace cetsp::postprocess
