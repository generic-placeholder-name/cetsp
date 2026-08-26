#include <cetsp/cetsp.hpp>
#include <cetsp/postprocess/socp.hpp>

#include <boost/geometry/algorithms/distance.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace bg = boost::geometry;

constexpr double maximumImportedBoundaryError = 1e-5;

struct ImportedTour {
    std::vector<std::size_t> circleIds;
    std::vector<Point> points;
};

std::vector<Circle> readCircles(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open input instance");
    }
    std::vector<Circle> circles;
    double x = 0.0;
    double y = 0.0;
    double radius = 0.0;
    while (input >> x >> y >> radius) {
        circles.push_back({Point{x, y}, radius});
    }
    if (!input.eof() || circles.empty()) {
        throw std::runtime_error("malformed or empty input instance");
    }
    return circles;
}

std::vector<Circle> readMaCircles(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open MA-CETSP instance");
    }
    std::size_t count = 0;
    if (!(input >> count) || count == 0) {
        throw std::runtime_error("malformed or empty MA-CETSP instance");
    }
    std::vector<Circle> circles;
    circles.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        double x = 0.0;
        double y = 0.0;
        double radius = 0.0;
        if (!(input >> x >> y >> radius)) {
            throw std::runtime_error("malformed MA-CETSP instance");
        }
        circles.push_back({Point{x, y}, radius});
    }
    input >> std::ws;
    if (!input.eof()) {
        throw std::runtime_error("MA-CETSP instance has trailing data");
    }
    return circles;
}

ImportedTour readResult(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open MA-CETSP result");
    }
    std::string ignored;
    std::getline(input, ignored);
    std::getline(input, ignored);
    ImportedTour tour;
    std::size_t id = 0;
    double x = 0.0;
    double y = 0.0;
    while (input >> id >> x >> y) {
        tour.circleIds.push_back(id);
        tour.points.emplace_back(x, y);
    }
    if (!input.eof() || tour.points.empty()) {
        throw std::runtime_error("malformed or empty MA-CETSP result");
    }
    return tour;
}

double repairImportedBoundaryError(
    ImportedTour& tour,
    const std::vector<Circle>& maCircles) {
    if (tour.points.size() != maCircles.size()) {
        throw std::runtime_error(
            "MA-CETSP result and instance have different sizes");
    }

    std::vector<bool> seen(maCircles.size(), false);
    double maximumRepair = 0.0;
    for (std::size_t index = 0; index < tour.points.size(); ++index) {
        const std::size_t circleId = tour.circleIds[index];
        if (circleId >= maCircles.size() || seen[circleId]) {
            throw std::runtime_error(
                "MA-CETSP result does not contain a circle permutation");
        }
        seen[circleId] = true;

        Point& point = tour.points[index];
        const Circle& circle = maCircles[circleId];
        const double distance = bg::distance(point, circle.center);
        const double excess = distance - circle.r;
        if (excess <= 0.0) {
            continue;
        }
        // Legacy MA-CETSP results were written with six significant digits.
        const double repairLimit = maximumImportedBoundaryError *
            std::max({
                1.0,
                std::abs(bg::get<0>(circle.center)),
                std::abs(bg::get<1>(circle.center)),
                circle.r,
            });
        if (!std::isfinite(distance) || excess > repairLimit) {
            throw std::runtime_error(
                "MA-CETSP result misses its assigned circle by more than the "
                "import repair limit");
        }

        if (circle.r == 0.0 || distance == 0.0) {
            point = circle.center;
        } else {
            const double targetRadius = std::nextafter(circle.r, 0.0);
            const double scale = targetRadius / distance;
            bg::set<0>(point, bg::get<0>(circle.center) +
                scale * (bg::get<0>(point) - bg::get<0>(circle.center)));
            bg::set<1>(point, bg::get<1>(circle.center) +
                scale * (bg::get<1>(point) - bg::get<1>(circle.center)));
        }
        maximumRepair = std::max(maximumRepair, excess);
    }
    return maximumRepair;
}

void writeResult(
    const std::filesystem::path& path,
    const std::vector<Point>& tour,
    double distance) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create optimized result");
    }
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t index = 0; index < tour.size(); ++index) {
        output << index << (index + 1 == tour.size() ? '\n' : ' ');
    }
    output << "value " << distance << " running_time 0\n";
    for (std::size_t index = 0; index < tour.size(); ++index) {
        output << index << ' '
               << boost::geometry::get<0>(tour[index]) << ' '
               << boost::geometry::get<1>(tour[index]) << '\n';
    }
}

long long elapsedMilliseconds(
    std::chrono::steady_clock::time_point start) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start).count();
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 7) {
            std::cerr
                << "usage: CETSP_ma_cetsp_result_check <instance.txt> "
                   "<ma-result.txt> [gurobi-executable|-] "
                   "[optimized-result] [max-threads] [ma-instance.txt]\n";
            return 2;
        }
        const std::vector<Circle> circles = readCircles(argv[1]);
        ImportedTour imported = readResult(argv[2]);
        const double maDistance = totalTourDistance(imported.points);
        const bool maValid = verifyTour(imported.points, circles);
        double maximumRepair = 0.0;
        if (!maValid && argc >= 7) {
            maximumRepair = repairImportedBoundaryError(
                imported, readMaCircles(argv[6]));
        }
        std::vector<Point> tour = std::move(imported.points);
        const double repairedDistance = totalTourDistance(tour);
        const bool repairedValid = verifyTour(tour, circles);
        long long socpMilliseconds = 0;
        if (repairedValid && argc >= 4 && std::string(argv[3]) != "-") {
            cetsp::postprocess::SocpOptions options;
            options.gurobiExecutable = argv[3];
            if (argc >= 6) {
                options.maxThreads = std::stoull(argv[5]);
            }
            const auto start = std::chrono::steady_clock::now();
            tour = cetsp::postprocess::improvePointsWithSocp(
                tour, circles, options);
            socpMilliseconds = elapsedMilliseconds(start);
        }
        const double finalDistance = totalTourDistance(tour);
        const bool finalValid = verifyTour(tour, circles);
        if (argc >= 5 && finalValid) {
            writeResult(argv[4], tour, finalDistance);
        }
        std::cout << std::fixed << std::setprecision(10)
                  << "ma_distance=" << maDistance << '\n'
                  << "ma_valid=" << std::boolalpha << maValid << '\n'
                  << "repair_distance=" << repairedDistance << '\n'
                  << "repair_maximum=" << maximumRepair << '\n'
                  << "repair_valid=" << repairedValid << '\n'
                  << "final_distance=" << finalDistance << '\n'
                  << "final_valid=" << finalValid << '\n'
                  << "socp_ms=" << socpMilliseconds << '\n';
        return repairedValid && finalValid ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "MA-CETSP result check failed: " << error.what() << '\n';
        return 2;
    }
}
