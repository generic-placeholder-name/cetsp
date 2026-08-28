#include <cetsp/cetsp.hpp>

#include "circle_geometry.hpp"
#include "merge.hpp"
#include "reconstruct.hpp"
#include "tour.hpp"

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/core/access.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
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
            expect(tree[id].isLeaf(),
                   "surviving input circles occupy the merge-tree leaves");
            if (id < expectedCenters.size()) {
                expectNear(bg::get<0>(tree[id].center), expectedCenters[id], 1e-12,
                           "leaf IDs follow reduced survivor order");
            }
        }
        for (std::size_t id = leafCount; id < tree.size(); ++id) {
            expect(!tree[id].isLeaf() &&
                       tree[id].left < id && tree[id].right < id,
                   "internal merge-tree children precede their parent");
        }
    }

    const auto tour = reconstructTour(tree);
    expect(verifyTour(tour, circles),
           "a tree built from compacted survivors reconstructs a valid tour");
}

void testReconstructionReinsertionCascade() {
    // Unmerging the root first creates leaf 0 and branch 3 as two tour nodes.
    // Inserting branch 3 consumes two of leaf 0's three energy units. When
    // branch 3 is unmerged, inserting leaf 1 consumes the last unit and forces
    // leaf 0 through the explicit delete/reinsert work stack.
    const std::vector<TreeNode> tree{
        TreeNode::leaf(Point{0.0, 0.0}, 0.0),
        TreeNode::leaf(Point{10.0, 0.0}, 0.0),
        TreeNode::leaf(Point{20.0, 0.0}, 0.0),
        TreeNode::branch(1, 2, 1.0, Point{15.0, 0.0}, 5.0),
        TreeNode::branch(0, 3, 2.0, Point{10.0, 0.0}, 10.0),
    };
    const std::vector<Circle> leaves{
        makeCircle(0.0, 0.0, 0.0),
        makeCircle(10.0, 0.0, 0.0),
        makeCircle(20.0, 0.0, 0.0),
    };

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
        const auto [center, radius] =
            makeCombinedCircle(
                firstCenter,
                firstRadius,
                secondCenter,
                secondRadius,
                randomEngine);
        static_cast<void>(center);
        expect(std::isfinite(radius), "the combined-circle radius is finite");
        expect(radius >= halfChord && radius <= overlapDepth,
               "the combined-circle radius lies between both geometric bounds");
    }
}

void testInternalBoundaryInputs() {
    std::mt19937_64 randomEngine{1234};
    expect(buildMergeTree({}, randomEngine).empty(),
           "an empty instance produces an empty merge tree");
    expect(reconstructTour({}).empty(),
           "an empty merge tree produces an empty tour");

    const std::vector<Circle> singleton{makeCircle(2.0, -4.0, 3.0)};
    const auto singletonTree = buildMergeTree(singleton, randomEngine);
    expect(singletonTree.size() == 1, "one circle produces one merge-tree leaf");
    if (singletonTree.size() == 1) {
        expect(singletonTree.front().isLeaf(),
               "the singleton merge-tree node has leaf semantics");
    }
}

} // namespace

int main() {
    testCoveringCircleReduction();
    testMergeStateIsInternal();
    testReconstructionReinsertionCascade();
    testTourStructuralMaintenance();
    testCombinedCircleRadiusRange();
    testInternalBoundaryInputs();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }

    std::cout << "All internal CETSP tests passed\n";
    return 0;
}
