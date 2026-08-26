#include <cetsp/postprocess/socp.hpp>

#include "process_support.hpp"

#include <cetsp/cetsp.hpp>

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/core/access.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cetsp::postprocess {
namespace {

namespace bg = boost::geometry;

constexpr double toleranceMultiplier = 64.0;

double circleTolerance(const Circle& circle) {
    return toleranceMultiplier * std::numeric_limits<double>::epsilon() *
        std::max({
            1.0,
            std::abs(bg::get<0>(circle.center)) + circle.r,
            std::abs(bg::get<1>(circle.center)) + circle.r,
        });
}

struct VisitGroup {
    Point point;
    std::vector<std::size_t> circleIndices;
};

std::vector<VisitGroup> assignCircles(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles) {
    if (!verifyTour(tour, circles)) {
        throw std::invalid_argument(
            "SOCP postprocessing requires a valid initial tour");
    }

    std::vector<VisitGroup> groups;
    groups.reserve(tour.size());
    for (const Point& point : tour) {
        groups.push_back({point, {}});
    }

    for (std::size_t circleIndex = 0;
         circleIndex < circles.size(); ++circleIndex) {
        const Circle& circle = circles[circleIndex];
        std::size_t bestVisit = tour.size();
        double bestDistance = std::numeric_limits<double>::infinity();
        for (std::size_t visitIndex = 0;
             visitIndex < tour.size(); ++visitIndex) {
            const double distance = bg::distance(
                tour[visitIndex], circle.center);
            if (distance < bestDistance) {
                bestDistance = distance;
                bestVisit = visitIndex;
            }
        }
        if (bestVisit == tour.size() ||
            bestDistance > circle.r + circleTolerance(circle)) {
            throw std::logic_error(
                "tour verification and SOCP circle assignment disagree");
        }
        groups[bestVisit].circleIndices.push_back(circleIndex);
    }

    std::erase_if(groups, [](const VisitGroup& group) {
        return group.circleIndices.empty();
    });
    return groups;
}

struct CoordinateTransform {
    double originX = 0.0;
    double originY = 0.0;
    double extent = 1.0;

    [[nodiscard]] double x(const Point& point) const {
        return (bg::get<0>(point) - originX) / extent;
    }
    [[nodiscard]] double y(const Point& point) const {
        return (bg::get<1>(point) - originY) / extent;
    }
    [[nodiscard]] Point point(double normalizedX, double normalizedY) const {
        return Point{
            originX + normalizedX * extent,
            originY + normalizedY * extent};
    }
};

CoordinateTransform chooseTransform(const std::vector<Circle>& circles) {
    double minimumX = std::numeric_limits<double>::infinity();
    double minimumY = std::numeric_limits<double>::infinity();
    double maximumX = -std::numeric_limits<double>::infinity();
    double maximumY = -std::numeric_limits<double>::infinity();
    for (const Circle& circle : circles) {
        const double x = bg::get<0>(circle.center);
        const double y = bg::get<1>(circle.center);
        minimumX = std::min(minimumX, x - circle.r);
        minimumY = std::min(minimumY, y - circle.r);
        maximumX = std::max(maximumX, x + circle.r);
        maximumY = std::max(maximumY, y + circle.r);
    }
    const double extent = std::max(maximumX - minimumX, maximumY - minimumY);
    return {
        minimumX,
        minimumY,
        extent > 0.0 ? extent : 1.0,
    };
}

void writeSignedLinearTerm(
    std::ostream& output,
    double coefficient,
    const std::string& variable) {
    output << (std::signbit(coefficient) ? " - " : " + ")
           << std::abs(coefficient) << ' ' << variable;
}

void writeModel(
    const std::filesystem::path& path,
    const std::vector<VisitGroup>& groups,
    const std::vector<Circle>& circles,
    const CoordinateTransform& transform) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create SOCP model file");
    }
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "Minimize\n obj:";
    for (std::size_t index = 0; index < groups.size(); ++index) {
        output << " + t" << index;
    }
    output << "\nSubject To\n";

    std::size_t constraintIndex = 0;
    for (std::size_t visitIndex = 0;
         visitIndex < groups.size(); ++visitIndex) {
        for (std::size_t circleIndex : groups[visitIndex].circleIndices) {
            const Circle& circle = circles[circleIndex];
            const double centerX = transform.x(circle.center);
            const double centerY = transform.y(circle.center);
            const double radius =
                (circle.r + circleTolerance(circle)) / transform.extent;
            output << " circle" << constraintIndex++ << ':';
            writeSignedLinearTerm(
                output, -2.0 * centerX,
                "x" + std::to_string(visitIndex));
            writeSignedLinearTerm(
                output, -2.0 * centerY,
                "y" + std::to_string(visitIndex));
            output << " + [ x" << visitIndex << " ^ 2 + y"
                   << visitIndex << " ^ 2 ] <= "
                   << radius * radius - centerX * centerX - centerY * centerY
                   << '\n';
        }
    }

    for (std::size_t index = 0; index < groups.size(); ++index) {
        const std::size_t next = (index + 1) % groups.size();
        output << " edge" << index
               << ": [ x" << index << " ^ 2 - 2 x" << index
               << " * x" << next << " + x" << next
               << " ^ 2 + y" << index << " ^ 2 - 2 y" << index
               << " * y" << next << " + y" << next
               << " ^ 2 - t" << index << " ^ 2 ] <= 0\n";
    }

    output << "Bounds\n";
    for (std::size_t index = 0; index < groups.size(); ++index) {
        output << " x" << index << " free\n"
               << " y" << index << " free\n"
               << " t" << index << " >= 0\n";
    }
    output << "End\n";
    if (!output) {
        throw std::runtime_error("failed while writing SOCP model file");
    }
}

std::vector<Point> readSolution(
    const std::filesystem::path& path,
    std::size_t visitCount,
    const CoordinateTransform& transform) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Gurobi did not produce an SOCP solution file");
    }

    std::vector<double> x(visitCount, 0.0);
    std::vector<double> y(visitCount, 0.0);
    std::string name;
    while (input >> name) {
        if (name.starts_with('#')) {
            std::string ignored;
            std::getline(input, ignored);
            continue;
        }
        double value = 0.0;
        if (!(input >> value)) {
            throw std::runtime_error("malformed Gurobi solution file");
        }
        if (name.size() < 2 || (name[0] != 'x' && name[0] != 'y')) {
            continue;
        }
        const std::size_t index = std::stoull(name.substr(1));
        if (index >= visitCount) {
            throw std::runtime_error(
                "Gurobi solution contains an invalid visit variable");
        }
        if (name[0] == 'x') {
            x[index] = value;
        } else {
            y[index] = value;
        }
    }
    if (input.bad()) {
        throw std::runtime_error("failed while reading Gurobi solution file");
    }

    std::vector<Point> solution;
    solution.reserve(visitCount);
    for (std::size_t index = 0; index < visitCount; ++index) {
        // Gurobi's plain .sol format is sparse: omitted variables are zero.
        solution.push_back(transform.point(x[index], y[index]));
    }
    return solution;
}

std::string numberArgument(double value) {
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << value;
    return output.str();
}

} // namespace

std::vector<Point> improvePointsWithSocp(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles,
    const SocpOptions& options) {
    if (tour.empty() || circles.empty()) {
        return circles.empty() ? std::vector<Point>{} : tour;
    }
    if (options.timeLimitSeconds < 0.0 ||
        !std::isfinite(options.timeLimitSeconds)) {
        throw std::invalid_argument(
            "SOCP time limit must be finite and non-negative");
    }

    const std::vector<VisitGroup> groups = assignCircles(tour, circles);
    std::vector<Point> prunedTour;
    prunedTour.reserve(groups.size());
    for (const VisitGroup& group : groups) {
        prunedTour.push_back(group.point);
    }
    if (groups.size() < 2) {
        return prunedTour;
    }

    const CoordinateTransform transform = chooseTransform(circles);
    detail::TemporaryDirectory workspace("cetsp-socp");
    const std::filesystem::path modelPath = workspace.path() / "model.lp";
    const std::filesystem::path solutionPath = workspace.path() / "solution.sol";
    const std::filesystem::path logPath = workspace.path() / "gurobi.log";
    writeModel(modelPath, groups, circles, transform);

    std::vector<std::string> arguments{
        "OutputFlag=0",
        "LogFile=" + logPath.string(),
        "ResultFile=" + solutionPath.string(),
        "FeasibilityTol=1e-9",
        "BarConvTol=1e-10",
    };
    if (options.maxThreads != 0) {
        arguments.push_back("Threads=" + std::to_string(options.maxThreads));
    }
    if (options.timeLimitSeconds != 0.0) {
        arguments.push_back(
            "TimeLimit=" + numberArgument(options.timeLimitSeconds));
    }
    arguments.push_back(modelPath.string());

    const int result = detail::runProcess(
        options.gurobiExecutable, arguments);
    if (result != 0) {
        throw std::runtime_error(
            "Gurobi exited unsuccessfully with code " + std::to_string(result));
    }

    std::vector<Point> candidate =
        readSolution(solutionPath, groups.size(), transform);
    if (!verifyTour(candidate, circles) ||
        totalTourDistance(candidate) >= totalTourDistance(prunedTour)) {
        return prunedTour;
    }
    return candidate;
}

} // namespace cetsp::postprocess
