#include "reconstruct.hpp"

#include "circle_geometry.hpp"
#include "debug.hpp"
#include "tour.hpp"

#include <boost/geometry.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bg = boost::geometry;

namespace {

constexpr std::size_t energyPerInsertion = 3;
using NodeId = MergeTree::NodeId;

[[nodiscard]] constexpr bool reachedPowerOfTwoCheckpoint(
    std::size_t value) noexcept {
    return value >= 2 && (value & (value - 1)) == 0;
}

void sortByRadius(
    std::vector<NodeId>& ids,
    const MergeTree& mergeTree) {
    std::sort(ids.begin(), ids.end(), [&](NodeId a, NodeId b) {
        return mergeTree.neighborhood(a).r < mergeTree.neighborhood(b).r;
    });
}

class TourReconstructor {
public:
    explicit TourReconstructor(const MergeTree& mergeTree)
        : mergeTree_(mergeTree), tour_(mergeTree.size()) {
        insertionStack_.reserve(mergeTree.size());
    }

    [[nodiscard]] std::vector<Point> reconstruct();

private:
    enum class InsertionPhase {
        begin,
        afterPreviousNeighbor,
        afterFollowingNeighbor,
    };

    struct InsertionStep {
        NodeId treeNodeId;
        InsertionPhase phase;
    };

    [[nodiscard]] std::vector<NodeId> processNeighbor(
        TourNodeHandle neighborHandle);
    void scheduleAfterNeighbor(
        NodeId treeNodeId,
        InsertionPhase continuation,
        TourNodeHandle neighborHandle);
    void beginCircleInsertion(NodeId treeNodeId);
    void resumeCircleInsertionAfterPreviousNeighbor(NodeId treeNodeId);
    void finishCircleInsertion(NodeId treeNodeId);
    void insertCircle(NodeId treeNodeId);

    const MergeTree& mergeTree_;
    Tour tour_;
    std::vector<InsertionStep> insertionStack_;
};

std::vector<NodeId> TourReconstructor::processNeighbor(
    TourNodeHandle neighborHandle) {
    if (!tour_.consumeEnergy(neighborHandle)) {
        return {};
    }

    DBG("Neighbor visit has zero energy; removing and reinserting its assignments.");
    std::vector<NodeId> assignments = tour_.eraseVisit(neighborHandle);
    sortByRadius(assignments, mergeTree_);
    return assignments;
}

void TourReconstructor::scheduleAfterNeighbor(
    NodeId treeNodeId,
    InsertionPhase continuation,
    TourNodeHandle neighborHandle) {
    std::vector<NodeId> reinsertions = processNeighbor(neighborHandle);

    // This explicitly encodes the old recursive order. The continuation sits
    // below every reinsertion, and reverse scheduling makes the first
    // radius-sorted assignment execute first on the LIFO stack.
    insertionStack_.push_back({treeNodeId, continuation});
    for (auto it = reinsertions.rbegin(); it != reinsertions.rend(); ++it) {
        insertionStack_.push_back({*it, InsertionPhase::begin});
    }
}

void TourReconstructor::beginCircleInsertion(NodeId treeNodeId) {
    const Circle& neighborhood = mergeTree_.neighborhood(treeNodeId);
    const Point& center = neighborhood.center;
    const double radius = neighborhood.r;

    DBG(tour_.size() << " points and segments in the tour indexes.");
    DBG("Inserting circle id=" << treeNodeId << " at ("
        << bg::get<0>(center) << ", " << bg::get<1>(center)
        << "), r=" << radius);

    if (tour_.empty()) {
        static_cast<void>(tour_.createFirstVisit(
            center, treeNodeId, energyPerInsertion));
        DBG("Created first tour visit for tree node " << treeNodeId);
        return;
    }

    DBG("Trying to link to an existing point via R-tree...");
    const MaybeTourNodeHandle existingHandle = tour_.nearestVisit(center);
    if (existingHandle) {
        const double distance = bg::distance(
            tour_.point(*existingHandle), center);
        if (distance <= radius) {
            tour_.addAssignment(
                *existingHandle, treeNodeId, energyPerInsertion);
            DBG("Linked tree node " << treeNodeId
                << " to an existing tour visit");
            scheduleAfterNeighbor(
                treeNodeId,
                InsertionPhase::afterPreviousNeighbor,
                tour_.previous(*existingHandle));
            return;
        }
    }

    DBG("Finding best edge to insert via segment R-tree...");
    constexpr std::size_t candidateCount = 16;
    double bestAddedDistance = std::numeric_limits<double>::infinity();
    std::optional<TourEdge> bestEdge;
    Point bestPoint;

    for (const TourEdge& edge :
         tour_.nearestEdges(center, candidateCount)) {
        const Point& leftPoint = edge.startPoint();
        const Point& rightPoint = edge.endPoint();
        Point candidatePoint = chooseInsertionPoint(
            neighborhood, leftPoint, rightPoint);
        const double addedDistance =
            bg::distance(leftPoint, candidatePoint) +
            bg::distance(candidatePoint, rightPoint) -
            bg::distance(leftPoint, rightPoint);
        if (addedDistance < bestAddedDistance) {
            bestAddedDistance = addedDistance;
            bestEdge = edge;
            bestPoint = candidatePoint;
        }
    }

    if (!bestEdge) {
        throw std::logic_error(
            "failed to find a tour edge for circle insertion");
    }

    const TourNodeHandle leftHandle = bestEdge->start();
    const TourNodeHandle rightHandle = bestEdge->end();
    DBG("Best edge for insertion has point ("
        << bg::get<0>(bestPoint) << ", " << bg::get<1>(bestPoint)
        << "), addCost=" << bestAddedDistance);

    const TourNodeHandle newHandle = tour_.insertVisitBetween(
        bestPoint,
        treeNodeId,
        leftHandle,
        rightHandle,
        energyPerInsertion);
    DBG("Inserted a new tour visit for tree node " << treeNodeId);

    scheduleAfterNeighbor(
        treeNodeId,
        InsertionPhase::afterPreviousNeighbor,
        tour_.previous(newHandle));
}

void TourReconstructor::resumeCircleInsertionAfterPreviousNeighbor(
    NodeId treeNodeId) {
    const MaybeTourNodeHandle currentHandle = tour_.visitFor(treeNodeId);
    if (!currentHandle) {
        throw std::logic_error(
            "tree node lost its tour assignment during previous-neighbor processing");
    }

    scheduleAfterNeighbor(
        treeNodeId,
        InsertionPhase::afterFollowingNeighbor,
        tour_.next(*currentHandle));
}

void TourReconstructor::finishCircleInsertion(NodeId treeNodeId) {
    const MaybeTourNodeHandle currentHandle = tour_.visitFor(treeNodeId);
    if (!currentHandle) {
        throw std::logic_error(
            "tree node lost its tour assignment during following-neighbor processing");
    }

    if (reachedPowerOfTwoCheckpoint(
            tour_.insertionCount(*currentHandle))) {
        tour_.optimizeVisit(*currentHandle, mergeTree_);
    }
}

void TourReconstructor::insertCircle(NodeId treeNodeId) {
    if (!insertionStack_.empty()) {
        throw std::logic_error("nested tour insertion stack is already active");
    }

    insertionStack_.push_back({treeNodeId, InsertionPhase::begin});
    while (!insertionStack_.empty()) {
        const InsertionStep step = insertionStack_.back();
        insertionStack_.pop_back();

        switch (step.phase) {
        case InsertionPhase::begin:
            beginCircleInsertion(step.treeNodeId);
            break;
        case InsertionPhase::afterPreviousNeighbor:
            resumeCircleInsertionAfterPreviousNeighbor(step.treeNodeId);
            break;
        case InsertionPhase::afterFollowingNeighbor:
            finishCircleInsertion(step.treeNodeId);
            break;
        }
    }
}

std::vector<Point> TourReconstructor::reconstruct() {
    const std::size_t nodeCount = mergeTree_.size();
    const std::optional<NodeId> root = mergeTree_.root();
    if (!root) {
        return {};
    }

    DBG("Starting unmerge with " << nodeCount << " tree nodes.");
    std::priority_queue<std::pair<double, NodeId>> pendingBranches;
    if (!mergeTree_.isLeaf(*root)) {
        pendingBranches.push({mergeTree_.mergeGap(*root), *root});
    }

    DBG("Inserting root node id=" << *root);
    insertCircle(*root);
    tour_.assertValid();

    std::size_t nodesProcessed = 0;
    while (!pendingBranches.empty()) {
        const auto [mergeGap, treeNodeId] = pendingBranches.top();
        pendingBranches.pop();
        ++nodesProcessed;
        if (mergeTree_.isLeaf(treeNodeId)) continue;

        DBG("Unmerging node id=" << treeNodeId
            << " (mergeGap=" << mergeGap
            << ")");

        tour_.removeAssignment(treeNodeId);
        for (NodeId child : mergeTree_.children(treeNodeId)) {
            DBG("Inserting child node id=" << child);
            insertCircle(child);
            if (!mergeTree_.isLeaf(child)) {
                pendingBranches.push({mergeTree_.mergeGap(child), child});
            }
        }

        // Maintenance is weighted by assigned-tree-node multiplicity. Scanning
        // every assignment visits shared tour nodes once per assigned tree node.
        if (reachedPowerOfTwoCheckpoint(nodesProcessed) ||
            pendingBranches.empty()) {
            for (NodeId id = nodeCount; id-- > 0;) {
                const MaybeTourNodeHandle handle = tour_.visitFor(id);
                if (!handle) continue;
                if (reachedPowerOfTwoCheckpoint(
                        tour_.recordInsertion(*handle))) {
                    DBG("Optimizing tour visit at ("
                        << bg::get<0>(tour_.point(*handle)) << ", "
                        << bg::get<1>(tour_.point(*handle)) << ")");
                    tour_.optimizeVisit(*handle, mergeTree_);
                }
            }

            for (NodeId id = nodeCount; id-- > 0;) {
                DBG("Processing node id=" << id
                    << " for energy reduction.");
                const MaybeTourNodeHandle handle = tour_.visitFor(id);
                if (!handle || !tour_.consumeEnergy(*handle)) continue;

                std::vector<NodeId> assignments =
                    tour_.eraseVisit(*handle);
                sortByRadius(assignments, mergeTree_);
                for (NodeId assignment : assignments) {
                    insertCircle(assignment);
                }
            }

            for (NodeId id = nodeCount; id-- > 0;) {
                const MaybeTourNodeHandle handle = tour_.visitFor(id);
                if (!handle) continue;
                if (reachedPowerOfTwoCheckpoint(
                        tour_.recordInsertion(*handle))) {
                    DBG("Optimizing tour visit at ("
                        << bg::get<0>(tour_.point(*handle)) << ", "
                        << bg::get<1>(tour_.point(*handle)) << ")");
                    tour_.optimizeVisit(*handle, mergeTree_);
                }
            }
        }

        tour_.assertValid();
    }

    std::vector<Point> result = tour_.points();
    DBG("Tour reconstruction complete. Tour size: " << result.size());
    tour_.assertValid();
    return result;
}

} // namespace

std::vector<Point> reconstructTour(const MergeTree& mergeTree) {
    TourReconstructor reconstructor(mergeTree);
    return reconstructor.reconstruct();
}
