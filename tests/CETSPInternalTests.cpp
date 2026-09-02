#include <cetsp/cetsp.hpp>

#include "circle_geometry.hpp"
#include "merge.hpp"
#include "reconstruct.hpp"
#include "tour.hpp"

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/core/access.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

namespace bg = boost::geometry;

int failures = 0;

void expect(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void expectNear(double actual, double expected, double tolerance, std::string_view message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr << "FAIL: " << message << " (actual=" << actual
                  << ", expected=" << expected << ")\n";
        ++failures;
    }
}

Circle makeCircle(double x, double y, double radius) {
    return {Point{x, y}, radius};
}

double pathLengthVia(
    const Point& edgeStart,
    const Point& point,
    const Point& edgeEnd) {
    return bg::distance(edgeStart, point) + bg::distance(point, edgeEnd);
}

Point closestPointOnSegmentForTest(
    const Point& point,
    const Point& segmentStart,
    const Point& segmentEnd) {
    const double startX = bg::get<0>(segmentStart);
    const double startY = bg::get<1>(segmentStart);
    const double dx = bg::get<0>(segmentEnd) - startX;
    const double dy = bg::get<1>(segmentEnd) - startY;
    const double lengthSquared = dx * dx + dy * dy;
    if (lengthSquared == 0.0) {
        return segmentStart;
    }

    const double offsetX = bg::get<0>(point) - startX;
    const double offsetY = bg::get<1>(point) - startY;
    const double fraction = std::clamp(
        (offsetX * dx + offsetY * dy) / lengthSquared,
        0.0,
        1.0);
    return Point{startX + fraction * dx, startY + fraction * dy};
}

double denseOptimalInsertionPathLength(
    const Point& center,
    double radius,
    const Point& edgeStart,
    const Point& edgeEnd) {
    const Point closestEdgePoint =
        closestPointOnSegmentForTest(center, edgeStart, edgeEnd);
    if (bg::distance(center, closestEdgePoint) <= radius) {
        return bg::distance(edgeStart, edgeEnd);
    }

    constexpr std::size_t sampleCount = 16384;
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < sampleCount; ++index) {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(index) /
            static_cast<double>(sampleCount);
        const Point point{
            bg::get<0>(center) + radius * std::cos(angle),
            bg::get<1>(center) + radius * std::sin(angle)};
        best = std::min(best, pathLengthVia(edgeStart, point, edgeEnd));
    }
    return best;
}

void testCoveringCircleReduction() {
    std::vector<Circle> nested{
        makeCircle(0.0, 0.0, 10.0),
        makeCircle(0.5, 0.0, 2.0),
        makeCircle(0.5, 0.0, 1.0),
    };

    removeCoveringCircles(nested);

    expect(nested.size() == 1, "nested neighborhoods reduce to the smallest circle");
    if (nested.size() == 1) {
        expectNear(nested.front().r, 1.0, 0.0,
                   "the larger, redundant neighborhoods are removed");
        expectNear(bg::get<0>(nested.front().center), 0.5, 0.0,
                   "covering reduction preserves the survivor's geometry");
    }

    std::vector<Circle> disjoint{
        makeCircle(-3.0, 0.0, 1.0),
        makeCircle(3.0, 0.0, 1.0),
    };
    removeCoveringCircles(disjoint);
    expect(disjoint.size() == 2, "disjoint neighborhoods are not removed");
}

void testMergeStateIsInternal() {
    std::vector<Circle> circles{
        makeCircle(0.0, 0.0, 5.0),
        makeCircle(20.0, 0.0, 1.0),
        makeCircle(0.0, 0.0, 1.0),
        makeCircle(40.0, 0.0, 5.0),
        makeCircle(40.0, 0.0, 1.0),
        makeCircle(60.0, 0.0, 1.0),
    };

    removeCoveringCircles(circles);
    expect(circles.size() == 4,
           "covering reduction can remove noncontiguous input positions");

    const std::vector<double> expectedCenters{20.0, 0.0, 40.0, 60.0};
    if (circles.size() == expectedCenters.size()) {
        for (std::size_t i = 0; i < circles.size(); ++i) {
            expectNear(bg::get<0>(circles[i].center), expectedCenters[i], 0.0,
                       "covering reduction preserves survivor order");
        }
    }

    const auto geometryBeforeMerge = circles;
    std::mt19937_64 randomEngine{0x6a09e667f3bcc909ULL};
    const auto tree = buildMergeTree(circles, randomEngine);

    expect(circles.size() == geometryBeforeMerge.size(),
           "building a merge tree does not resize caller-owned circles");
    for (std::size_t i = 0; i < circles.size(); ++i) {
        expect(bg::equals(circles[i].center, geometryBeforeMerge[i].center) &&
                   circles[i].r == geometryBeforeMerge[i].r,
               "building a merge tree does not mutate input geometry");
    }

    const std::size_t leafCount = circles.size();
    expect(tree.size() == 2 * leafCount - 1,
           "a merge tree has one internal node per merge");
    if (tree.size() == 2 * leafCount - 1) {
        for (std::size_t id = 0; id < leafCount; ++id) {
            expect(tree.isLeaf(id),
                   "surviving input circles occupy the merge-tree leaves");
            if (id < expectedCenters.size()) {
                expectNear(bg::get<0>(tree.neighborhood(id).center), expectedCenters[id], 1e-12,
                           "leaf IDs follow reduced survivor order");
            }
        }
        for (std::size_t id = leafCount; id < tree.size(); ++id) {
            const auto children = tree.children(id);
            expect(!tree.isLeaf(id) &&
                       children[0] < id && children[1] < id,
                   "internal merge-tree children precede their parent");
        }
    }

    const auto tour = reconstructTour(tree);
    expect(verifyTour(tour, circles),
           "a tree built from compacted survivors reconstructs a valid tour");
}

void testThreeLeafReconstruction() {
    // Exercise multiple unmerge levels using a history produced through the
    // same construction path as the solver.
    const std::vector<Circle> leaves{
        makeCircle(0.0, 0.0, 0.0),
        makeCircle(10.0, 0.0, 0.0),
        makeCircle(20.0, 0.0, 0.0),
    };
    std::mt19937_64 randomEngine{0x6a09e667f3bcc909ULL};
    const MergeTree tree = buildMergeTree(leaves, randomEngine);

    const auto tour = reconstructTour(tree);
    expect(tour.size() == leaves.size(),
           "a reinsertion cascade retains one tour point per point neighborhood");
    expect(verifyTour(tour, leaves),
           "a reinsertion cascade preserves every leaf neighborhood");
    expectNear(totalTourDistance(tour), 40.0, 1e-12,
               "a reinsertion cascade preserves the expected cyclic tour");
}

void testTourStructuralMaintenance() {
    Tour tour(3);
    expect(tour.empty(), "a new tour is empty");

    const TourNodeHandle first = tour.createFirstVisit(
        Point{0.0, 0.0}, 0, 3);
    tour.addAssignment(first, 1, 3);
    expect(tour.size() == 1,
           "multiple tree-node assignments can share one tour visit");
    expect(tour.previous(first) == first && tour.next(first) == first,
           "the first tour visit forms a singleton cycle");
    const auto singletonEdges = tour.nearestEdges(Point{1.0, 0.0}, 1);
    expect(singletonEdges.size() == 1 &&
               singletonEdges.front().start() == first &&
               singletonEdges.front().end() == first &&
               bg::equals(singletonEdges.front().startPoint(),
                          tour.point(first)) &&
               bg::equals(singletonEdges.front().endPoint(),
                          tour.point(first)),
           "edge queries expose complete singleton-edge semantics");
    expect(tour.visitFor(0) == first && tour.visitFor(1) == first,
           "the assignment index resolves every shared tree node");

    bool duplicateAssignmentRejected = false;
    try {
        tour.addAssignment(first, 1, 3);
    } catch (const std::logic_error&) {
        duplicateAssignmentRejected = true;
    }
    expect(duplicateAssignmentRejected,
           "a tree node cannot be assigned to the tour twice");

    bool invalidTreeNodeRejected = false;
    try {
        tour.addAssignment(first, 3, 3);
    } catch (const std::out_of_range&) {
        invalidTreeNodeRejected = true;
    }
    expect(invalidTreeNodeRejected,
           "assignment IDs must belong to the tour's merge tree");

    const TourNodeHandle second = tour.insertVisitBetween(
        Point{10.0, 0.0}, 2, first, first, 3);
    expect(tour.size() == 2, "inserting on the singleton edge grows the tour");
    expect(tour.next(first) == second && tour.previous(first) == second &&
               tour.next(second) == first && tour.previous(second) == first,
           "insertion updates both directions of the cyclic topology");
    tour.assertValid();

    tour.removeAssignment(0);
    expect(tour.size() == 2 && !tour.visitFor(0) &&
               tour.visitFor(1) == first,
           "removing one shared assignment retains its tour visit");

    const auto erasedAssignments = tour.eraseVisit(first);
    expect(erasedAssignments.size() == 1 && erasedAssignments.front() == 1,
           "erasing a visit returns and clears its remaining assignments");
    expect(tour.size() == 1 && !tour.visitFor(1),
           "erasing a visit updates the assignment index");
    expect(tour.previous(second) == second && tour.next(second) == second,
           "erasing from a two-visit tour restores a singleton cycle");

    bool staleHandleRejected = false;
    try {
        static_cast<void>(tour.point(first));
    } catch (const std::logic_error&) {
        staleHandleRejected = true;
    }
    expect(staleHandleRejected,
           "erased tour handles are rejected after their generation advances");

    tour.removeAssignment(2);
    expect(tour.empty(), "removing the final assignment empties the tour");
    tour.assertValid();
}

void testInsertionPointSelection() {
    const Point center{0.0, 0.0};

    const Point zeroRadiusResult = chooseInsertionPoint(
        Circle{center, 0.0}, Point{-2.0, 1.0}, Point{3.0, 4.0});
    expect(bg::equals(zeroRadiusResult, center),
           "a zero-radius neighborhood inserts its center");

    const Point crossingStart{-2.0, 0.0};
    const Point crossingEnd{2.0, 0.0};
    const Point crossingResult = chooseInsertionPoint(
        Circle{center, 1.0}, crossingStart, crossingEnd);
    expectNear(
        pathLengthVia(crossingStart, crossingResult, crossingEnd),
        bg::distance(crossingStart, crossingEnd),
        1e-12,
        "an edge crossing the neighborhood has zero insertion cost");

    const Point tangentStart{-2.0, 1.0};
    const Point tangentEnd{2.0, 1.0};
    const Point tangentResult = chooseInsertionPoint(
        Circle{center, 1.0}, tangentStart, tangentEnd);
    expectNear(
        pathLengthVia(tangentStart, tangentResult, tangentEnd),
        bg::distance(tangentStart, tangentEnd),
        1e-12,
        "a tangent edge has zero insertion cost");

    const Point repeatedEndpoint{2.0, 0.0};
    const Point repeatedResult = chooseInsertionPoint(
        Circle{center, 1.0}, repeatedEndpoint, repeatedEndpoint);
    expectNear(
        pathLengthVia(repeatedEndpoint, repeatedResult, repeatedEndpoint),
        2.0,
        1e-10,
        "a degenerate edge is minimized at the nearest boundary point");

    bool negativeRadiusRejected = false;
    try {
        static_cast<void>(chooseInsertionPoint(
            Circle{center, -1.0}, crossingStart, crossingEnd));
    } catch (const std::invalid_argument&) {
        negativeRadiusRejected = true;
    }
    expect(negativeRadiusRejected,
           "insertion-point geometry rejects a negative radius");

    const Point regressionStart{
        -0.6256443365812765,
        -0.812701981837634};
    const Point regressionEnd{
        -2.0757443438448773,
        -8.142396907134557};
    const Point regressionResult = chooseInsertionPoint(
        Circle{center, 1.0}, regressionStart, regressionEnd);
    const double regressionOracle = denseOptimalInsertionPathLength(
        center, 1.0, regressionStart, regressionEnd);
    expect(
        pathLengthVia(regressionStart, regressionResult, regressionEnd) <=
            regressionOracle + 1e-4,
        "asymmetric insertion stays close to a dense boundary-search oracle");

    const Point antipodalStart{-100.0, 1.01};
    const Point antipodalEnd{100.0, 1.01};
    const Point antipodalResult = chooseInsertionPoint(
        Circle{center, 1.0}, antipodalStart, antipodalEnd);
    expectNear(
        bg::get<0>(antipodalResult),
        0.0,
        1e-10,
        "nearly antipodal endpoints preserve the symmetric boundary point");
    expectNear(
        bg::get<1>(antipodalResult),
        1.0,
        1e-10,
        "nearly antipodal endpoints select the nearer boundary arc");

    std::mt19937_64 randomEngine{0x243f6a8885a308d3ULL};
    std::uniform_real_distribution<double> centerDistribution(-20.0, 20.0);
    std::uniform_real_distribution<double> radiusDistribution(0.1, 5.0);
    std::uniform_real_distribution<double> angleDistribution(
        0.0, 2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> distanceFactorDistribution(
        1.01, 10.0);
    double maximumRelativePathExcess = 0.0;
    for (std::size_t sample = 0; sample < 128; ++sample) {
        const Point randomCenter{
            centerDistribution(randomEngine),
            centerDistribution(randomEngine)};
        const double radius = radiusDistribution(randomEngine);
        const double firstAngle = angleDistribution(randomEngine);
        const double secondAngle = angleDistribution(randomEngine);
        const double firstDistance =
            radius * distanceFactorDistribution(randomEngine);
        const double secondDistance =
            radius * distanceFactorDistribution(randomEngine);
        const Point edgeStart{
            bg::get<0>(randomCenter) +
                firstDistance * std::cos(firstAngle),
            bg::get<1>(randomCenter) +
                firstDistance * std::sin(firstAngle)};
        const Point edgeEnd{
            bg::get<0>(randomCenter) +
                secondDistance * std::cos(secondAngle),
            bg::get<1>(randomCenter) +
                secondDistance * std::sin(secondAngle)};

        const Point result = chooseInsertionPoint(
            Circle{randomCenter, radius}, edgeStart, edgeEnd);
        const double resultDistance = bg::distance(result, randomCenter);
        const double actualPathLength =
            pathLengthVia(edgeStart, result, edgeEnd);
        const double oraclePathLength = denseOptimalInsertionPathLength(
            randomCenter, radius, edgeStart, edgeEnd);
        const double numericalTolerance =
            1e-9 * (1.0 + oraclePathLength);
        maximumRelativePathExcess = std::max(
            maximumRelativePathExcess,
            (actualPathLength - oraclePathLength) /
                (1.0 + oraclePathLength));

        expect(
            std::isfinite(bg::get<0>(result)) &&
                std::isfinite(bg::get<1>(result)),
            "randomized insertion geometry returns a finite point");
        expect(
            resultDistance <= radius + 1e-10 * (1.0 + radius),
            "randomized insertion geometry returns a point in the neighborhood");
        expect(
            actualPathLength + numericalTolerance >=
                bg::distance(edgeStart, edgeEnd),
            "randomized insertion geometry respects the triangle inequality");
    }
    expectNear(
        maximumRelativePathExcess,
        0.0,
        2e-3,
        "one-step insertion geometry stays close to a dense angular oracle");
}

void testCombinedCircleRadiusRange() {
    const Point firstCenter{0.0, 0.0};
    const Point secondCenter{1.1, 0.0};
    constexpr double firstRadius = 5.0;
    constexpr double secondRadius = 4.0;

    // Here overlapDepth is greater than the intersection half-chord. This used
    // to construct uniform_real_distribution with its bounds reversed.
    const double distance = bg::distance(firstCenter, secondCenter);
    const double overlapDepth = (firstRadius + secondRadius - distance) / 2.0;
    const double chordOffset =
        (firstRadius * firstRadius - secondRadius * secondRadius + distance * distance) /
        (2.0 * distance);
    const double halfChord =
        std::sqrt(firstRadius * firstRadius - chordOffset * chordOffset);
    expect(overlapDepth > halfChord,
           "the combined-circle regression case has reversed raw bounds");

    std::mt19937_64 randomEngine{0x4d595df4d0f33173ULL};
    for (int sample = 0; sample < 32; ++sample) {
        const Circle combined = makeCombinedCircle(
            Circle{firstCenter, firstRadius},
            Circle{secondCenter, secondRadius},
            randomEngine);
        expect(std::isfinite(combined.r), "the combined-circle radius is finite");
        expect(combined.r >= halfChord && combined.r <= overlapDepth,
               "the combined-circle radius lies between both geometric bounds");
    }
}

void testInternalBoundaryInputs() {
    std::mt19937_64 randomEngine{1234};
    expect(buildMergeTree({}, randomEngine).empty(),
           "an empty instance produces an empty merge tree");
    expect(reconstructTour(MergeTree{}).empty(),
           "an empty merge tree produces an empty tour");

    const std::vector<Circle> singleton{makeCircle(2.0, -4.0, 3.0)};
    const auto singletonTree = buildMergeTree(singleton, randomEngine);
    expect(singletonTree.size() == 1, "one circle produces one merge-tree leaf");
    if (singletonTree.size() == 1) {
        expect(singletonTree.isLeaf(*singletonTree.root()),
               "the singleton merge-tree node has leaf semantics");
    }
}

} // namespace

int main() {
    testCoveringCircleReduction();
    testMergeStateIsInternal();
    testThreeLeafReconstruction();
    testTourStructuralMaintenance();
    testInsertionPointSelection();
    testCombinedCircleRadiusRange();
    testInternalBoundaryInputs();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }

    std::cout << "All internal CETSP tests passed\n";
    return 0;
}
