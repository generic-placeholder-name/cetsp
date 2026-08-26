#include <cetsp/cetsp.hpp>

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/core/access.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace bg = boost::geometry;

double maximumCoverageExcess(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles) {
    double maximumExcess = -std::numeric_limits<double>::infinity();
    for (const auto& circle : circles) {
        double nearestDistance = std::numeric_limits<double>::infinity();
        for (const auto& point : tour) {
            nearestDistance = std::min(nearestDistance, bg::distance(point, circle.center));
        }
        maximumExcess = std::max(maximumExcess, nearestDistance - circle.r);
    }
    return circles.empty() ? 0.0 : maximumExcess;
}

std::uint64_t makeRandomSeed() {
    std::random_device entropy;
    const auto high = static_cast<std::uint64_t>(entropy());
    const auto low = static_cast<std::uint64_t>(entropy());
    return (high << 32U) ^ low;
}

std::string exactTourBits(const std::vector<Point>& tour) {
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const auto& point : tour) {
        output << std::setw(16)
               << std::bit_cast<std::uint64_t>(bg::get<0>(point))
               << ':'
               << std::setw(16)
               << std::bit_cast<std::uint64_t>(bg::get<1>(point))
               << ',';
    }
    return output.str();
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 5) {
            std::cerr <<
                "usage: CETSP_benchmark <instance.txt> [repetitions] [seed] "
                "[max-threads]\n";
            return 2;
        }

        const int repetitions = argc >= 3 ? std::stoi(argv[2]) : 10;
        const std::uint64_t seed =
            argc >= 4 ? std::stoull(argv[3]) : makeRandomSeed();
        const std::size_t maxThreads =
            argc == 5 ? std::stoull(argv[4]) : 0;
        std::ifstream input(argv[1]);
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
        if (!input.eof()) {
            throw std::runtime_error("malformed input instance");
        }

        CetspOptions options;
        options.numRepeats = repetitions;
        options.seed = seed;
        options.maxThreads = maxThreads;

        const auto start = std::chrono::steady_clock::now();
        const auto tour = solveCetsp(circles, options);
        const auto stop = std::chrono::steady_clock::now();
        const auto elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        const bool valid = verifyTour(tour, circles);

        std::cout << std::fixed << std::setprecision(10)
                  << "seed=" << seed << '\n'
                  << "max_threads="
                  << (maxThreads == 0 ? "auto" : std::to_string(maxThreads))
                  << '\n'
                  << "circles=" << circles.size() << '\n'
                  << "tour_points=" << tour.size() << '\n'
                  << "tour_bits=" << exactTourBits(tour) << '\n'
                  << "valid=" << std::boolalpha << valid << '\n'
                  << "maximum_coverage_excess="
                  << maximumCoverageExcess(tour, circles) << '\n'
                  << "distance=" << totalTourDistance(tour) << '\n'
                  << "elapsed_ms=" << elapsed.count() << '\n';

        return valid ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "benchmark failed: " << error.what() << '\n';
        return 2;
    }
}
