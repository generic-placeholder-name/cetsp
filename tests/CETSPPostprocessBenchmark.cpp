#include <cetsp/cetsp.hpp>
#include <cetsp/postprocess/lkh.hpp>
#include <cetsp/postprocess/socp.hpp>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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
    if (!input.eof()) {
        throw std::runtime_error("malformed input instance");
    }
    return circles;
}

long long millisecondsSince(
    std::chrono::steady_clock::time_point start) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start).count();
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 8) {
            std::cerr
                << "usage: CETSP_postprocess_benchmark <instance.txt> "
                   "[repetitions] [seed] [max-threads] "
                   "[lkh-executable|-] [gurobi-executable|-] [rounds]\n";
            return 2;
        }

        const int repetitions = argc >= 3 ? std::stoi(argv[2]) : 10;
        const std::uint64_t seed =
            argc >= 4 ? std::stoull(argv[3]) : 123456789ULL;
        const std::size_t maxThreads =
            argc >= 5 ? std::stoull(argv[4]) : 0;
        const std::string lkhExecutable = argc >= 6 ? argv[5] : "-";
        const std::string gurobiExecutable = argc >= 7 ? argv[6] : "-";
        const std::size_t rounds = argc >= 8 ? std::stoull(argv[7]) : 1;
        if (rounds == 0) {
            throw std::invalid_argument(
                "postprocessing round count must be positive");
        }

        const std::vector<Circle> circles = readCircles(argv[1]);
        CetspOptions solverOptions;
        solverOptions.numRepeats = repetitions;
        solverOptions.seed = seed;
        solverOptions.maxThreads = maxThreads;

        const auto solverStart = std::chrono::steady_clock::now();
        std::vector<Point> tour = solveCetsp(circles, solverOptions);
        const long long solverMilliseconds = millisecondsSince(solverStart);
        const double baseDistance = totalTourDistance(tour);

        double afterLkhDistance = baseDistance;
        double afterSocpDistance = baseDistance;
        long long lkhMilliseconds = 0;
        long long socpMilliseconds = 0;

        for (std::size_t round = 0; round < rounds; ++round) {
            if (lkhExecutable != "-") {
                cetsp::postprocess::LkhOptions options;
                options.executable = lkhExecutable;
                options.seed = seed + round;
                const auto start = std::chrono::steady_clock::now();
                tour = cetsp::postprocess::improveOrderWithLkh(tour, options);
                lkhMilliseconds += millisecondsSince(start);
                afterLkhDistance = totalTourDistance(tour);
            }
            if (gurobiExecutable != "-") {
                cetsp::postprocess::SocpOptions options;
                options.gurobiExecutable = gurobiExecutable;
                const auto start = std::chrono::steady_clock::now();
                tour = cetsp::postprocess::improvePointsWithSocp(
                    tour, circles, options);
                socpMilliseconds += millisecondsSince(start);
                afterSocpDistance = totalTourDistance(tour);
            }
        }

        std::cout << std::fixed << std::setprecision(10)
                  << "instance=" << std::filesystem::path(argv[1]).filename().string() << '\n'
                  << "circles=" << circles.size() << '\n'
                  << "tour_points=" << tour.size() << '\n'
                  << "base_distance=" << baseDistance << '\n'
                  << "after_lkh_distance=" << afterLkhDistance << '\n'
                  << "after_socp_distance=" << afterSocpDistance << '\n'
                  << "final_distance=" << totalTourDistance(tour) << '\n'
                  << "valid=" << std::boolalpha << verifyTour(tour, circles) << '\n'
                  << "solver_ms=" << solverMilliseconds << '\n'
                  << "lkh_ms=" << lkhMilliseconds << '\n'
                  << "socp_ms=" << socpMilliseconds << '\n';
        return verifyTour(tour, circles) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "postprocess benchmark failed: " << error.what() << '\n';
        return 2;
    }
}
