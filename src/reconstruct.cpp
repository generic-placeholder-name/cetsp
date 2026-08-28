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

[[nodiscard]] constexpr bool reachedPowerOfTwoCheckpoint(
    std::size_t value) noexcept {
    return value >= 2 && (value & (value - 1)) == 0;
}

void sortByRadius(
    std::vector<TreeNodeId>& ids,
    const std::vector<TreeNode>& treeNodes) {
    std::sort(ids.begin(), ids.end(), [&](TreeNodeId a, TreeNodeId b) {
        return treeNodes[a].r < treeNodes[b].r;
    });
}

std::string formatHandle(TourNodeHandle handle) {
    return std::to_string(handle.slot) + ":" +
           std::to_string(handle.generation);
}

class TourReconstructor {
public:
    explicit TourReconstructor(const std::vector<TreeNode>& treeNodes)
        : treeNodes_(treeNodes), tour_(treeNodes.size()) {
        insertionStack_.reserve(treeNodes.size());
    }

    [[nodiscard]] std::vector<Point> reconstruct();

private:
    enum class InsertionPhase {
        begin,
        afterPreviousNeighbor,
        afterFollowingNeighbor,
    };

    struct InsertionStep {
        TreeNodeId treeNodeId;
        InsertionPhase phase;
    };

    [[nodiscard]] std::vector<TreeNodeId> processNeighbor(
        TourNodeHandle neighborHandle);
    void scheduleAfterNeighbor(
        TreeNodeId treeNodeId,
        InsertionPhase continuation,
        TourNodeHandle neighborHandle);
    void beginCircleInsertion(TreeNodeId treeNodeId);
    void resumeCircleInsertionAfterPreviousNeighbor(TreeNodeId treeNodeId);
    void finishCircleInsertion(TreeNodeId treeNodeId);
    void insertCircle(TreeNodeId treeNodeId);

    const std::vector<TreeNode>& treeNodes_;
    Tour tour_;
    std::vector<InsertionStep> insertionStack_;
};

std::vector<TreeNodeId> TourReconstructor::processNeighbor(
    TourNodeHandle neighborHandle) {
    if (!tour_.consumeEnergy(neighborHandle)) {
        return {};
    }

    DBG("Neighbor handle=" << formatHandle(neighborHandle)
        << " has zero energy, removing and reinserting its assignments.");
    std::vector<TreeNodeId> assignments = tour_.eraseVisit(neighborHandle);
    sortByRadius(assignments, treeNodes_);
    return assignments;
}

void TourReconstructor::scheduleAfterNeighbor(
    TreeNodeId treeNodeId,
    InsertionPhase continuation,
    TourNodeHandle neighborHandle) {
    std::vector<TreeNodeId> reinsertions = processNeighbor(neighborHandle);

    // This explicitly encodes the old recursive order. The continuation sits
    // below every reinsertion, and reverse scheduling makes the first
    // radius-sorted assignment execute first on the LIFO stack.
    insertionStack_.push_back({treeNodeId, continuation});
    for (auto it = reinsertions.rbegin(); it != reinsertions.rend(); ++it) {
        insertionStack_.push_back({*it, InsertionPhase::begin});
    }
}

void TourReconstructor::beginCircleInsertion(TreeNodeId treeNodeId) {
    const TreeNode& treeNode = treeNodes_[treeNodeId];
    const Point& center = treeNode.center;
    const double radius = treeNode.r;

    DBG(tour_.size() << " points and segments in the tour indexes.");
    DBG("Inserting circle id=" << treeNodeId << " at ("
        << bg::get<0>(center) << ", " << bg::get<1>(center)
        << "), r=" << radius);

    if (tour_.empty()) {
        const TourNodeHandle handle = tour_.createFirstVisit(
            center, treeNodeId, energyPerInsertion);
        DBG("Created first TourNode handle=" << formatHandle(handle)
            << " for tree node " << treeNodeId);
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
                << " to existing handle " << formatHandle(*existingHandle));
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
    MaybeTourNodeHandle bestLeft;
    Point bestPoint;

    for (TourNodeHandle leftHandle :
         tour_.nearestEdgeStarts(center, candidateCount)) {
        const TourNodeHandle rightHandle = tour_.next(leftHandle);
        const Point leftPoint = tour_.point(leftHandle);
        const Point rightPoint = tour_.point(rightHandle);
        Point candidatePoint = chooseInsertionPoint(
            center, radius, leftPoint, rightPoint);
        const double addedDistance =
            bg::distance(leftPoint, candidatePoint) +
            bg::distance(candidatePoint, rightPoint) -
            bg::distance(leftPoint, rightPoint);
        if (addedDistance < bestAddedDistance) {
            bestAddedDistance = addedDistance;
            bestLeft = leftHandle;
            bestPoint = candidatePoint;
        }
    }

    if (!bestLeft) {
        throw std::logic_error(
            "failed to find a tour edge for circle insertion");
    }

    const TourNodeHandle leftHandle = *bestLeft;
    const TourNodeHandle rightHandle = tour_.next(leftHandle);
    DBG("Best edge for insertion: between handle="
        << formatHandle(leftHandle) << " and handle="
        << formatHandle(rightHandle) << " at point ("
        << bg::get<0>(bestPoint) << ", " << bg::get<1>(bestPoint)
        << "), addCost=" << bestAddedDistance);

    const TourNodeHandle newHandle = tour_.insertVisitBetween(
        bestPoint,
        treeNodeId,
        leftHandle,
        rightHandle,
        energyPerInsertion);
    DBG("Inserted new TourNode handle=" << formatHandle(newHandle)
        << " for tree node " << treeNodeId << " between handles "
        << formatHandle(leftHandle) << " and "
        << formatHandle(rightHandle));

    scheduleAfterNeighbor(
        treeNodeId,
        InsertionPhase::afterPreviousNeighbor,
        tour_.previous(newHandle));
}

void TourReconstructor::resumeCircleInsertionAfterPreviousNeighbor(
    TreeNodeId treeNodeId) {
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

void TourReconstructor::finishCircleInsertion(TreeNodeId treeNodeId) {
    const MaybeTourNodeHandle currentHandle = tour_.visitFor(treeNodeId);
    if (!currentHandle) {
        throw std::logic_error(
            "tree node lost its tour assignment during following-neighbor processing");
    }

    if (reachedPowerOfTwoCheckpoint(
            tour_.insertionCount(*currentHandle))) {
        tour_.optimizeVisit(*currentHandle, treeNodes_);
    }
}

void TourReconstructor::insertCircle(TreeNodeId treeNodeId) {
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
    const std::size_t nodeCount = treeNodes_.size();
    if (treeNodes_.empty()) {
        return {};
    }

    DBG("Starting unmerge with " << nodeCount << " tree nodes.");
    std::priority_queue<std::pair<double, TreeNodeId>> pendingBranches;
    pendingBranches.push(
        {treeNodes_[nodeCount - 1].mergeGap, nodeCount - 1});

    DBG("Inserting root node id=" << (nodeCount - 1));
    insertCircle(nodeCount - 1);
    tour_.assertValid();

    std::size_t nodesProcessed = 0;
    while (!pendingBranches.empty()) {
        const auto [mergeGap, treeNodeId] = pendingBranches.top();
        pendingBranches.pop();
        ++nodesProcessed;
        if (treeNodes_[treeNodeId].isLeaf()) continue;

        DBG("Unmerging node id=" << treeNodeId
            << " (mergeGap=" << mergeGap
            << "), left=" << treeNodes_[treeNodeId].left
            << ", right=" << treeNodes_[treeNodeId].right);

        tour_.removeAssignment(treeNodeId);
        for (TreeNodeId child : treeNodes_[treeNodeId].children()) {
            DBG("Inserting child node id=" << child);
            insertCircle(child);
            if (!treeNodes_[child].isLeaf()) {
                pendingBranches.push({treeNodes_[child].mergeGap, child});
            }
        }

        // Maintenance is weighted by assigned-tree-node multiplicity. Scanning
        // every assignment visits shared tour nodes once per assigned tree node.
        if (reachedPowerOfTwoCheckpoint(nodesProcessed) ||
            pendingBranches.empty()) {
            for (TreeNodeId id = nodeCount; id-- > 0;) {
                const MaybeTourNodeHandle handle = tour_.visitFor(id);
                if (!handle) continue;
                if (reachedPowerOfTwoCheckpoint(
                        tour_.recordInsertion(*handle))) {
                    DBG("Optimizing TourNode handle="
                        << formatHandle(*handle) << " at ("
                        << bg::get<0>(tour_.point(*handle)) << ", "
                        << bg::get<1>(tour_.point(*handle)) << ")");
                    tour_.optimizeVisit(*handle, treeNodes_);
                }
            }

            for (TreeNodeId id = nodeCount; id-- > 0;) {
                DBG("Processing node id=" << id
                    << " for energy reduction.");
                const MaybeTourNodeHandle handle = tour_.visitFor(id);
                if (!handle || !tour_.consumeEnergy(*handle)) continue;

                std::vector<TreeNodeId> assignments =
                    tour_.eraseVisit(*handle);
                sortByRadius(assignments, treeNodes_);
                for (TreeNodeId assignment : assignments) {
                    insertCircle(assignment);
                }
            }

            for (TreeNodeId id = nodeCount; id-- > 0;) {
                const MaybeTourNodeHandle handle = tour_.visitFor(id);
                if (!handle) continue;
                if (reachedPowerOfTwoCheckpoint(
                        tour_.recordInsertion(*handle))) {
                    DBG("Optimizing TourNode handle="
                        << formatHandle(*handle) << " at ("
                        << bg::get<0>(tour_.point(*handle)) << ", "
                        << bg::get<1>(tour_.point(*handle)) << ")");
                    tour_.optimizeVisit(*handle, treeNodes_);
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

std::vector<Point> reconstructTour(const std::vector<TreeNode>& treeNodes) {
    TourReconstructor reconstructor(treeNodes);
    return reconstructor.reconstruct();
}
