#include <cetsp/cetsp.hpp>

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/core/access.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

namespace bg = boost::geometry;

constexpr double floatingPointToleranceMultiplier = 64.0;
constexpr double fanOutValidationRelativeTolerance = 1e-10;

struct SeedNode {
    std::size_t id = 0;
    Point point;
};

struct Candidate {
    double distance = std::numeric_limits<double>::infinity();
    std::size_t repeatIndex = 0;
    std::vector<SeedNode> nodes;
};

struct VisitGroup {
    Point point;
    std::vector<std::size_t> circles;
};

struct ReducedInstance {
    std::vector<Circle> circles;
    std::size_t removed = 0;
};

double coordinate(const Point& point, std::size_t dimension) {
    return dimension == 0 ? bg::get<0>(point) : bg::get<1>(point);
}

Point interpolate(const Point& from, const Point& to, double fraction) {
    return Point{
        coordinate(from, 0) + fraction *
            (coordinate(to, 0) - coordinate(from, 0)),
        coordinate(from, 1) + fraction *
            (coordinate(to, 1) - coordinate(from, 1)),
    };
}

ReducedInstance removeDominatedCircles(
    const std::vector<Circle>& circles,
    std::size_t depotIndex) {
    std::vector<bool> removed(circles.size(), false);
    const Point& depot = circles[depotIndex].center;
    for (std::size_t outerIndex = 0;
         outerIndex < circles.size(); ++outerIndex) {
        if (outerIndex == depotIndex) {
            continue;
        }
        const Circle& outer = circles[outerIndex];
        if (bg::distance(depot, outer.center) <= outer.r) {
            removed[outerIndex] = true;
            continue;
        }
        for (std::size_t innerIndex = 0;
             innerIndex < circles.size(); ++innerIndex) {
            if (innerIndex == outerIndex || innerIndex == depotIndex) {
                continue;
            }
            const Circle& inner = circles[innerIndex];
            const double centerDistance =
                bg::distance(inner.center, outer.center);
            const bool contained =
                centerDistance + inner.r <= outer.r;
            const bool sameSize = inner.r == outer.r;
            if (contained &&
                (!sameSize || innerIndex < outerIndex)) {
                removed[outerIndex] = true;
                break;
            }
        }
    }

    ReducedInstance reduced;
    reduced.circles.reserve(circles.size());
    reduced.circles.push_back(circles[depotIndex]);
    for (std::size_t index = 0; index < circles.size(); ++index) {
        if (index != depotIndex && !removed[index]) {
            reduced.circles.push_back(circles[index]);
        }
    }
    reduced.removed = circles.size() - reduced.circles.size();
    return reduced;
}

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

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

double rayReach(
    const Point& point,
    const Point& toward,
    const Circle& circle,
    double capacity) {
    const double edgeLength = bg::distance(point, toward);
    if (edgeLength == 0.0 || capacity <= 0.0) {
        return 0.0;
    }
    const double ux = (coordinate(toward, 0) - coordinate(point, 0)) /
        edgeLength;
    const double uy = (coordinate(toward, 1) - coordinate(point, 1)) /
        edgeLength;
    const double dx = coordinate(point, 0) - coordinate(circle.center, 0);
    const double dy = coordinate(point, 1) - coordinate(circle.center, 1);
    const double projection = dx * ux + dy * uy;
    const double discriminant = std::max(
        0.0,
        projection * projection -
            (dx * dx + dy * dy - circle.r * circle.r));
    const double boundary = -projection + std::sqrt(discriminant);
    return std::clamp(boundary, 0.0, std::min(capacity, edgeLength));
}

Point moveToward(const Point& point, const Point& toward, double distance) {
    const double edgeLength = bg::distance(point, toward);
    if (edgeLength == 0.0 || distance == 0.0) {
        return point;
    }
    return interpolate(point, toward, distance / edgeLength);
}

std::vector<VisitGroup> assignCircles(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles,
    std::size_t depotIndex) {
    if (!verifyTour(tour, circles)) {
        throw std::logic_error("construction produced an invalid tour");
    }

    std::size_t depotVisit = tour.size();
    double depotDistance = std::numeric_limits<double>::infinity();
    for (std::size_t visit = 0; visit < tour.size(); ++visit) {
        const double distance = bg::distance(tour[visit], circles[depotIndex].center);
        if (distance < depotDistance) {
            depotDistance = distance;
            depotVisit = visit;
        }
    }
    if (depotVisit == tour.size() ||
        depotDistance > circles[depotIndex].r + 1e-9) {
        throw std::logic_error("tour does not contain its depot");
    }

    std::vector<VisitGroup> groups;
    groups.reserve(tour.size());
    for (const Point& point : tour) {
        groups.push_back({point, {}});
    }
    groups[depotVisit].circles.push_back(depotIndex);

    for (std::size_t circleIndex = 0;
         circleIndex < circles.size(); ++circleIndex) {
        if (circleIndex == depotIndex) {
            continue;
        }
        const Circle& circle = circles[circleIndex];
        std::size_t bestVisit = tour.size();
        double bestDistance = std::numeric_limits<double>::infinity();
        for (std::size_t visit = 0; visit < tour.size(); ++visit) {
            const double distance = bg::distance(tour[visit], circle.center);
            const double tolerance = floatingPointToleranceMultiplier *
                std::numeric_limits<double>::epsilon() *
                std::max({
                    1.0,
                    std::abs(coordinate(circle.center, 0)) + circle.r,
                    std::abs(coordinate(circle.center, 1)) + circle.r,
                });
            if (distance <= circle.r + tolerance && distance < bestDistance) {
                bestVisit = visit;
                bestDistance = distance;
            }
        }
        if (bestVisit == tour.size()) {
            throw std::logic_error("could not assign a circle to a covering visit");
        }
        groups[bestVisit].circles.push_back(circleIndex);
    }

    std::erase_if(groups, [](const VisitGroup& group) {
        return group.circles.empty();
    });
    return groups;
}

std::vector<std::size_t> makeMaIds(
    std::size_t circleCount,
    std::size_t depotIndex) {
    std::vector<std::size_t> ids(circleCount);
    ids[depotIndex] = 0;
    std::size_t next = 1;
    for (std::size_t index = 0; index < circleCount; ++index) {
        if (index != depotIndex) {
            ids[index] = next++;
        }
    }
    return ids;
}

std::vector<SeedNode> fanOut(
    const std::vector<Point>& tour,
    const std::vector<Circle>& circles,
    std::size_t depotIndex,
    std::uint64_t seed) {
    std::vector<VisitGroup> groups = assignCircles(tour, circles, depotIndex);
    const std::vector<std::size_t> maIds = makeMaIds(circles.size(), depotIndex);
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> splitDistribution(0.35, 0.65);
    std::uniform_real_distribution<double> unitDistribution(0.0, 1.0);
    std::uniform_real_distribution<double> extentDistribution(0.8, 1.0);

    std::vector<double> edgeSplits(groups.size());
    for (double& split : edgeSplits) {
        split = splitDistribution(random);
    }

    std::vector<std::vector<SeedNode>> expanded(groups.size());
    for (std::size_t groupIndex = 0;
         groupIndex < groups.size(); ++groupIndex) {
        const std::size_t previous =
            (groupIndex + groups.size() - 1) % groups.size();
        const std::size_t next = (groupIndex + 1) % groups.size();
        const VisitGroup& group = groups[groupIndex];
        const double incomingLength = bg::distance(
            groups[previous].point, group.point);
        const double outgoingLength = bg::distance(
            group.point, groups[next].point);
        const double incomingCapacity =
            (1.0 - edgeSplits[previous]) * incomingLength;
        const double outgoingCapacity =
            edgeSplits[groupIndex] * outgoingLength;

        std::size_t anchor = group.circles.front();
        if (std::find(group.circles.begin(), group.circles.end(), depotIndex) !=
            group.circles.end()) {
            anchor = depotIndex;
        } else {
            double leastFreedom = std::numeric_limits<double>::infinity();
            for (std::size_t circleIndex : group.circles) {
                const double incoming = rayReach(
                    group.point, groups[previous].point,
                    circles[circleIndex], incomingCapacity);
                const double outgoing = rayReach(
                    group.point, groups[next].point,
                    circles[circleIndex], outgoingCapacity);
                const double freedom = std::max(incoming, outgoing);
                if (freedom < leastFreedom) {
                    leastFreedom = freedom;
                    anchor = circleIndex;
                }
            }
        }

        std::vector<std::pair<double, SeedNode>> incoming;
        std::vector<std::pair<double, SeedNode>> outgoing;
        for (std::size_t circleIndex : group.circles) {
            if (circleIndex == anchor) {
                continue;
            }
            const double backward = rayReach(
                group.point, groups[previous].point,
                circles[circleIndex], incomingCapacity);
            const double forward = rayReach(
                group.point, groups[next].point,
                circles[circleIndex], outgoingCapacity);
            const double totalReach = backward + forward;
            const bool useIncoming = totalReach > 0.0 &&
                unitDistribution(random) < backward / totalReach;
            const double reach = useIncoming ? backward : forward;
            const double displacement = reach * extentDistribution(random);
            const Point point = moveToward(
                group.point,
                useIncoming ? groups[previous].point : groups[next].point,
                displacement);
            auto& side = useIncoming ? incoming : outgoing;
            side.push_back({displacement, {maIds[circleIndex], point}});
        }

        std::sort(incoming.begin(), incoming.end(),
            [](const auto& left, const auto& right) {
                return left.first > right.first;
            });
        std::sort(outgoing.begin(), outgoing.end(),
            [](const auto& left, const auto& right) {
                return left.first < right.first;
            });

        auto& nodes = expanded[groupIndex];
        nodes.reserve(group.circles.size());
        for (const auto& entry : incoming) {
            nodes.push_back(entry.second);
        }
        nodes.push_back({maIds[anchor], group.point});
        for (const auto& entry : outgoing) {
            nodes.push_back(entry.second);
        }
    }

    std::vector<SeedNode> nodes;
    nodes.reserve(circles.size());
    for (auto& group : expanded) {
        nodes.insert(
            nodes.end(),
            std::make_move_iterator(group.begin()),
            std::make_move_iterator(group.end()));
    }
    const auto depot = std::find_if(nodes.begin(), nodes.end(),
        [](const SeedNode& node) { return node.id == 0; });
    if (depot == nodes.end()) {
        throw std::logic_error("fan-out lost the depot");
    }
    std::rotate(nodes.begin(), depot, nodes.end());

    if (nodes.size() != circles.size()) {
        throw std::logic_error("fan-out did not emit one node per circle");
    }
    std::vector<bool> seen(nodes.size(), false);
    std::vector<std::size_t> circleForMaId(nodes.size());
    for (std::size_t circleIndex = 0;
         circleIndex < maIds.size(); ++circleIndex) {
        circleForMaId[maIds[circleIndex]] = circleIndex;
    }
    for (const SeedNode& node : nodes) {
        if (node.id >= nodes.size() || seen[node.id]) {
            throw std::logic_error("fan-out emitted duplicate circle IDs");
        }
        seen[node.id] = true;
        const Circle& circle = circles[circleForMaId[node.id]];
        const double scale = std::max({
            1.0,
            std::abs(coordinate(circle.center, 0)) + circle.r,
            std::abs(coordinate(circle.center, 1)) + circle.r,
        });
        if (bg::distance(node.point, circle.center) >
            circle.r + fanOutValidationRelativeTolerance * scale) {
            throw std::logic_error("fan-out placed a node outside its circle");
        }
    }
    return nodes;
}

double seedDistance(const std::vector<SeedNode>& nodes) {
    double distance = 0.0;
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        distance += bg::distance(
            nodes[index].point,
            nodes[(index + 1) % nodes.size()].point);
    }
    return distance;
}

double tourDistance(const std::vector<Point>& tour) {
    double distance = 0.0;
    for (std::size_t index = 0; index < tour.size(); ++index) {
        distance += bg::distance(tour[index], tour[(index + 1) % tour.size()]);
    }
    return distance;
}

double edgeDistance(
    const std::vector<SeedNode>& left,
    const std::vector<SeedNode>& right) {
    std::vector<std::size_t> rightPosition(right.size());
    for (std::size_t index = 0; index < right.size(); ++index) {
        rightPosition[right[index].id] = index;
    }
    std::size_t common = 0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        const std::size_t a = left[index].id;
        const std::size_t b = left[(index + 1) % left.size()].id;
        const std::size_t rightIndex = rightPosition[a];
        const std::size_t rightNext =
            right[(rightIndex + 1) % right.size()].id;
        const std::size_t rightPrevious =
            right[(rightIndex + right.size() - 1) % right.size()].id;
        if (rightNext == b || rightPrevious == b) {
            ++common;
        }
    }
    return 100.0 * static_cast<double>(left.size() - common) /
        static_cast<double>(left.size());
}

std::vector<Candidate> selectPool(
    std::vector<Candidate> candidates,
    std::size_t poolSize,
    double minimumDistance) {
    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& left, const Candidate& right) {
            return left.distance < right.distance ||
                (left.distance == right.distance &&
                 left.repeatIndex < right.repeatIndex);
        });
    std::vector<Candidate> selected;
    selected.reserve(std::min(poolSize, candidates.size()));
    std::vector<bool> used(candidates.size(), false);
    for (std::size_t index = 0;
         index < candidates.size() && selected.size() < poolSize; ++index) {
        bool diverse = true;
        for (const Candidate& existing : selected) {
            if (edgeDistance(candidates[index].nodes, existing.nodes) <=
                minimumDistance) {
                diverse = false;
                break;
            }
        }
        if (diverse) {
            used[index] = true;
            selected.push_back(std::move(candidates[index]));
        }
    }
    for (std::size_t index = 0;
         index < candidates.size() && selected.size() < poolSize; ++index) {
        if (!used[index]) {
            selected.push_back(std::move(candidates[index]));
        }
    }
    return selected;
}

void writeInstance(
    const std::filesystem::path& path,
    const std::vector<Circle>& circles,
    std::size_t depotIndex) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create MA-CETSP instance file");
    }
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << circles.size() << '\n';
    const auto writeCircle = [&](const Circle& circle) {
        output << coordinate(circle.center, 0) << ' '
               << coordinate(circle.center, 1) << ' '
               << circle.r << '\n';
    };
    writeCircle(circles[depotIndex]);
    for (std::size_t index = 0; index < circles.size(); ++index) {
        if (index != depotIndex) {
            writeCircle(circles[index]);
        }
    }
}

void writeSeeds(
    const std::filesystem::path& path,
    const std::vector<Candidate>& candidates,
    std::size_t nodeCount) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create MA-CETSP seed file");
    }
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << "MA_CETSP_SEEDS 1 " << nodeCount << ' '
           << candidates.size() << '\n';
    for (const Candidate& candidate : candidates) {
        output << "TOUR\n";
        for (const SeedNode& node : candidate.nodes) {
            output << node.id << ' '
                   << coordinate(node.point, 0) << ' '
                   << coordinate(node.point, 1) << '\n';
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4 || argc > 10) {
            std::cerr
                << "usage: CETSP_ma_cetsp_seed_export <instance.txt> "
                   "<ma-instance.txt> <seeds.txt> [repetitions] [pool-size] "
                   "[seed] [max-threads] [depot-index] [min-distance]\n";
            return 2;
        }
        const std::vector<Circle> circles = readCircles(argv[1]);
        const std::size_t repetitions =
            argc >= 5 ? std::stoull(argv[4]) : 1000;
        const std::size_t poolSize =
            argc >= 6 ? std::stoull(argv[5]) : 20;
        const std::uint64_t seed =
            argc >= 7 ? std::stoull(argv[6]) : 123456789ULL;
        std::size_t workerCount =
            argc >= 8 ? std::stoull(argv[7]) : 0;
        const std::size_t depotIndex =
            argc >= 9 ? std::stoull(argv[8]) : circles.size() - 1;
        const double minimumDistance =
            argc >= 10 ? std::stod(argv[9]) : 5.0;
        if (repetitions == 0 || poolSize == 0 || depotIndex >= circles.size()) {
            throw std::invalid_argument("invalid export options");
        }
        const ReducedInstance reduced =
            removeDominatedCircles(circles, depotIndex);
        if (workerCount == 0) {
            workerCount = std::max(1U, std::thread::hardware_concurrency());
        }
        workerCount = std::min(workerCount, repetitions);

        std::vector<std::optional<Candidate>> results(repetitions);
        std::atomic<std::size_t> nextRepeat = 0;
        std::atomic<bool> stop = false;
        std::mutex exceptionMutex;
        std::exception_ptr workerException;
        const auto worker = [&] {
            try {
                while (!stop.load(std::memory_order_relaxed)) {
                    const std::size_t index =
                        nextRepeat.fetch_add(1, std::memory_order_relaxed);
                    if (index >= repetitions) {
                        return;
                    }
                    const std::uint64_t repeatSeed = splitmix64(seed + index);
                    CetspOptions options;
                    options.numRepeats = 1;
                    options.seed = repeatSeed;
                    options.maxThreads = 1;
                    const std::vector<Point> tour = solveCetsp(circles, options);
                    if (!verifyTour(tour, circles)) {
                        throw std::logic_error(
                            "construction produced an invalid tour");
                    }
                    std::vector<SeedNode> nodes = fanOut(
                        tour, reduced.circles, 0,
                        splitmix64(repeatSeed));
                    const double expandedDistance = seedDistance(nodes);
                    const double sourceDistance = tourDistance(tour);
                    const double distanceTolerance =
                        fanOutValidationRelativeTolerance *
                        std::max(1.0, sourceDistance);
                    if (expandedDistance > sourceDistance + distanceTolerance) {
                        throw std::logic_error(
                            "fan-out increased the source tour length");
                    }
                    results[index] = Candidate{
                        expandedDistance, index, std::move(nodes)};
                }
            } catch (...) {
                std::lock_guard lock(exceptionMutex);
                if (!workerException) {
                    workerException = std::current_exception();
                }
                stop.store(true, std::memory_order_relaxed);
            }
        };

        std::vector<std::jthread> workers;
        workers.reserve(workerCount);
        for (std::size_t index = 0; index < workerCount; ++index) {
            workers.emplace_back(worker);
        }
        workers.clear();
        if (workerException) {
            std::rethrow_exception(workerException);
        }

        std::vector<Candidate> candidates;
        candidates.reserve(repetitions);
        for (auto& result : results) {
            if (result) {
                candidates.push_back(std::move(*result));
            }
        }
        std::vector<Candidate> selected = selectPool(
            std::move(candidates), std::min(poolSize, repetitions),
            minimumDistance);
        writeInstance(argv[2], reduced.circles, 0);
        writeSeeds(argv[3], selected, reduced.circles.size());

        const auto [minimum, maximum] = std::minmax_element(
            selected.begin(), selected.end(),
            [](const Candidate& left, const Candidate& right) {
                return left.distance < right.distance;
            });
        std::cout << std::fixed << std::setprecision(10)
                  << "generated=" << repetitions << '\n'
                  << "selected=" << selected.size() << '\n'
                  << "original_circles=" << circles.size() << '\n'
                  << "active_circles=" << reduced.circles.size() << '\n'
                  << "dominated_circles=" << reduced.removed << '\n'
                  << "best_distance=" << minimum->distance << '\n'
                  << "worst_selected_distance=" << maximum->distance << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "MA-CETSP seed export failed: " << error.what() << '\n';
        return 2;
    }
}
