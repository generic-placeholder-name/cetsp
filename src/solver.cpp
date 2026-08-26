#include <cetsp/cetsp.hpp>

#include "circle_geometry.hpp"
#include "debug.hpp"
#include "merge.hpp"
#include "reconstruct.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {

// Magic-number multiplier for accumulated floating-point geometry error.
constexpr double floatingPointToleranceMultiplier = 64.0;

double circleTolerance(const Circle& circle) {
    const double centerX = std::abs(bg::get<0>(circle.center));
    const double centerY = std::abs(bg::get<1>(circle.center));
    return floatingPointToleranceMultiplier *
           std::numeric_limits<double>::epsilon() * std::max({
        1.0,
        centerX + circle.r,
        centerY + circle.r,
    });
}

std::mt19937_64 makeRandomEngine(const std::optional<std::uint64_t>& seed) {
    if (seed) {
        return std::mt19937_64(*seed);
    }

    std::random_device entropy;
    std::seed_seq seedSequence{
        entropy(), entropy(), entropy(), entropy(),
        entropy(), entropy(), entropy(), entropy(),
    };
    return std::mt19937_64(seedSequence);
}

std::vector<std::uint64_t> makeRepeatSeeds(
    std::size_t repeatCount,
    const std::optional<std::uint64_t>& seed) {
    auto seedGenerator = makeRandomEngine(seed);
    std::vector<std::uint64_t> repeatSeeds(repeatCount);
    for (std::uint64_t& repeatSeed : repeatSeeds) {
        repeatSeed = seedGenerator();
    }
    return repeatSeeds;
}

std::size_t chooseWorkerCount(
    std::size_t repeatCount,
    std::size_t maxThreads) noexcept {
    std::size_t available = maxThreads;
    if (available == 0) {
        available = std::thread::hardware_concurrency();
        if (available == 0) {
            available = 1;
        }
    }
    return std::min(repeatCount, available);
}

struct RepeatResult {
    double distance = std::numeric_limits<double>::infinity();
    std::size_t repeatIndex = std::numeric_limits<std::size_t>::max();
    std::vector<Point> tour;
};

bool isBetterResult(
    double distance,
    std::size_t repeatIndex,
    const RepeatResult& incumbent) noexcept {
    return distance < incumbent.distance ||
           (distance == incumbent.distance &&
            repeatIndex < incumbent.repeatIndex);
}

void validateCircles(const std::vector<Circle>& circles) {
    for (const auto& circle : circles) {
        if (!std::isfinite(bg::get<0>(circle.center)) ||
            !std::isfinite(bg::get<1>(circle.center)) ||
            !std::isfinite(circle.r) || circle.r < 0.0) {
            throw std::invalid_argument(
                "CETSP circles must have finite coordinates and non-negative radii");
        }
    }
}

} // namespace

// Verify that the tour is valid: each circle should contain a point in the tour.
bool verifyTour(const std::vector<Point>& tour, const std::vector<Circle>& circles) {
    bgi::rtree<Point, bgi::rstar<16>> rtree(tour.begin(), tour.end());

    for (const auto& c : circles) {
        const double tolerance = circleTolerance(c);
        const auto predicate = bgi::within(circleBox(c, tolerance));
        bool covered = false;
        for (auto candidate = rtree.qbegin(predicate);
             candidate != rtree.qend(); ++candidate) {
            if (bg::distance(*candidate, c.center) <= c.r + tolerance) {
                covered = true;
                break;
            }
        }
        if (!covered) {
            DBG("Circle with center (" << bg::get<0>(c.center) << ", " << bg::get<1>(c.center)
                << ") and radius " << c.r << " does not contain any tour points.");
            return false;
        }
    }

    return true;
}

double totalTourDistance(const std::vector<Point>& tour) {
    if (tour.empty()) return 0.0;

    double totalDist = 0.0;
    for (size_t i = 0; i < tour.size(); ++i) {
        const Point& p1 = tour[i];
        const Point& p2 = tour[(i + 1) % tour.size()]; // Wrap around to the start
        totalDist += bg::distance(p1, p2);
    }
    return totalDist;
}

std::vector<Point> solveCetsp(
    const std::vector<Circle>& circles,
    const CetspOptions& options) {
    validateCircles(circles);
    if (options.numRepeats <= 0) {
        throw std::invalid_argument("CETSP repetition count must be positive");
    }
    if (circles.empty()) {
        return {};
    }

    auto circlesCopy = circles;
    // Step 1: Remove covering circles to simplify the problem
    removeCoveringCircles(circlesCopy);

    // If there's only one circle left, directly return its center as the tour
    if (circlesCopy.size() == 1) {
        DBG("Only one circle remains after removing coverings. Returning its center.");
        return {circlesCopy[0].center};
    }

    const std::size_t repeatCount =
        static_cast<std::size_t>(options.numRepeats);
    const std::size_t workerCount =
        chooseWorkerCount(repeatCount, options.maxThreads);
    const std::vector<std::uint64_t> repeatSeeds =
        makeRepeatSeeds(repeatCount, options.seed);

    std::vector<RepeatResult> workerResults(workerCount);
    std::atomic<std::size_t> nextRepeat = 0;
    std::atomic<bool> stopWorkers = false;
    std::mutex exceptionMutex;
    std::exception_ptr workerException;

    auto runWorker = [&](std::size_t workerIndex) {
        try {
            RepeatResult& workerBest = workerResults[workerIndex];
            while (!stopWorkers.load(std::memory_order_relaxed)) {
                const std::size_t repeatIndex =
                    nextRepeat.fetch_add(1, std::memory_order_relaxed);
                if (repeatIndex >= repeatCount) {
                    break;
                }

                std::mt19937_64 randomEngine(repeatSeeds[repeatIndex]);
                std::vector<TreeNode> mergeTree =
                    buildMergeTree(circlesCopy, randomEngine);
                std::vector<Point> tour = reconstructTour(mergeTree);

                const double distance = totalTourDistance(tour);
                if (!isBetterResult(distance, repeatIndex, workerBest) ||
                    !verifyTour(tour, circles)) {
                    continue;
                }

                workerBest.distance = distance;
                workerBest.repeatIndex = repeatIndex;
                workerBest.tour = std::move(tour);
            }
        } catch (...) {
            {
                std::lock_guard lock(exceptionMutex);
                if (!workerException) {
                    workerException = std::current_exception();
                }
            }
            stopWorkers.store(true, std::memory_order_relaxed);
        }
    };

    if (workerCount == 1) {
        runWorker(0);
    } else {
        std::vector<std::jthread> workers;
        workers.reserve(workerCount);
        for (std::size_t workerIndex = 0;
             workerIndex < workerCount; ++workerIndex) {
            workers.emplace_back(runWorker, workerIndex);
        }
    }

    if (workerException) {
        std::rethrow_exception(workerException);
    }

    RepeatResult bestResult;
    for (RepeatResult& workerResult : workerResults) {
        if (isBetterResult(
                workerResult.distance,
                workerResult.repeatIndex,
                bestResult)) {
            bestResult = std::move(workerResult);
        }
    }

    if (bestResult.tour.empty()) {
        DBG("Tour verification failed.");
    }
    else {
        DBG("Tour verification passed. Number of elements: " << circles.size());
    }

    // Return the valid tour
    return bestResult.tour;
}

std::vector<Point> CETSP(const std::vector<Circle>& circles, int numRepeats) {
    CetspOptions options;
    options.numRepeats = numRepeats;
    return solveCetsp(circles, options);
}
