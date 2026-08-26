#include <cetsp/postprocess/lkh.hpp>

#include "lkh_parser.hpp"
#include "process_support.hpp"

#include <cetsp/cetsp.hpp>

#include <boost/geometry/core/access.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cetsp::postprocess {
namespace {

namespace bg = boost::geometry;

struct CoordinateScale {
    double minimumX = 0.0;
    double minimumY = 0.0;
    double multiplier = 1.0;
};

CoordinateScale chooseCoordinateScale(const std::vector<Point>& points) {
    double minimumX = std::numeric_limits<double>::infinity();
    double minimumY = std::numeric_limits<double>::infinity();
    double maximumX = -std::numeric_limits<double>::infinity();
    double maximumY = -std::numeric_limits<double>::infinity();
    for (const Point& point : points) {
        minimumX = std::min(minimumX, bg::get<0>(point));
        minimumY = std::min(minimumY, bg::get<1>(point));
        maximumX = std::max(maximumX, bg::get<0>(point));
        maximumY = std::max(maximumY, bg::get<1>(point));
    }

    // LKH stores EUC_2D coordinates times PRECISION (100) in signed integers.
    constexpr double targetExtent = 1'000'000.0;
    const double extent = std::max(maximumX - minimumX, maximumY - minimumY);
    return {
        minimumX,
        minimumY,
        extent > 0.0 ? targetExtent / extent : 1.0,
    };
}

void writeProblem(
    const std::filesystem::path& path,
    const std::vector<Point>& points) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create LKH problem file");
    }

    const CoordinateScale scale = chooseCoordinateScale(points);
    output << "NAME : CETSP_POSTPROCESS\n"
           << "TYPE : TSP\n"
           << "DIMENSION : " << points.size() << '\n'
           << "EDGE_WEIGHT_TYPE : EUC_2D\n"
           << "NODE_COORD_SECTION\n";
    for (std::size_t index = 0; index < points.size(); ++index) {
        const auto x = static_cast<std::int64_t>(std::llround(
            (bg::get<0>(points[index]) - scale.minimumX) * scale.multiplier));
        const auto y = static_cast<std::int64_t>(std::llround(
            (bg::get<1>(points[index]) - scale.minimumY) * scale.multiplier));
        output << index + 1 << ' ' << x << ' ' << y << '\n';
    }
    output << "EOF\n";
    if (!output) {
        throw std::runtime_error("failed while writing LKH problem file");
    }
}

void writeParameters(
    const std::filesystem::path& path,
    const std::filesystem::path& problemPath,
    const std::filesystem::path& tourPath,
    const LkhOptions& options) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create LKH parameter file");
    }
    output << "PROBLEM_FILE = " << problemPath.string() << '\n'
           << "OUTPUT_TOUR_FILE = " << tourPath.string() << '\n'
           << "RUNS = " << options.runs << '\n'
           << "TRACE_LEVEL = 0\n";
    if (options.seed) {
        output << "SEED = " << *options.seed << '\n';
    }
    if (!output) {
        throw std::runtime_error("failed while writing LKH parameter file");
    }
}

std::vector<std::size_t> readTour(
    const std::filesystem::path& path,
    std::size_t expectedNodeCount) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("LKH did not produce a tour file");
    }
    return detail::parseLkhTour(input, expectedNodeCount);
}

} // namespace

std::vector<Point> improveOrderWithLkh(
    const std::vector<Point>& tour,
    const LkhOptions& options) {
    if (tour.size() < 4) {
        return tour;
    }
    if (options.runs == 0) {
        throw std::invalid_argument("LKH run count must be positive");
    }

    detail::TemporaryDirectory workspace("cetsp-lkh");
    const std::filesystem::path problemPath = workspace.path() / "problem.tsp";
    const std::filesystem::path parameterPath = workspace.path() / "parameters.par";
    const std::filesystem::path outputTourPath = workspace.path() / "output.tour";
    writeProblem(problemPath, tour);
    writeParameters(parameterPath, problemPath, outputTourPath, options);

    const int result = detail::runProcess(
        options.executable, {parameterPath.string()});
    if (result != 0) {
        throw std::runtime_error(
            "LKH exited unsuccessfully with code " + std::to_string(result));
    }

    const std::vector<std::size_t> permutation =
        readTour(outputTourPath, tour.size());
    std::vector<Point> candidate;
    candidate.reserve(tour.size());
    for (std::size_t index : permutation) {
        candidate.push_back(tour[index]);
    }

    return totalTourDistance(candidate) < totalTourDistance(tour)
        ? candidate
        : tour;
}

} // namespace cetsp::postprocess
