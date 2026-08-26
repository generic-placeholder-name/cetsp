#include <cetsp/postprocess/pipeline.hpp>

#include <cetsp/cetsp.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cetsp::postprocess {

std::vector<Point> polishTour(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles,
    const PostprocessOptions& options) {
    if (options.rounds == 0) {
        throw std::invalid_argument(
            "postprocessing round count must be positive");
    }
    if (!verifyTour(tour, circles)) {
        throw std::invalid_argument(
            "postprocessing requires a valid initial tour");
    }

    std::vector<Point> current = tour;
    for (std::size_t round = 0; round < options.rounds; ++round) {
        const double before = totalTourDistance(current);
        if (options.lkh) {
            current = improveOrderWithLkh(current, *options.lkh);
        }
        if (options.socp) {
            current = improvePointsWithSocp(current, circles, *options.socp);
        }
        const double after = totalTourDistance(current);
        const double scale = std::max({1.0, std::abs(before), std::abs(after)});
        if (before - after <=
            16.0 * std::numeric_limits<double>::epsilon() * scale) {
            break;
        }
    }
    return current;
}

} // namespace cetsp::postprocess
