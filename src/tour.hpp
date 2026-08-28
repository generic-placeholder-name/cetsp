#pragma once

#include <cetsp/types.hpp>

#include "id_set.hpp"
#include "merge_tree.hpp"

#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

struct TourNodeHandle {
    std::size_t slot;
    std::size_t generation;

    constexpr TourNodeHandle(
        std::size_t slotIndex,
        std::size_t generationValue) noexcept
        : slot(slotIndex), generation(generationValue) {}

    friend constexpr bool operator==(
        const TourNodeHandle&,
        const TourNodeHandle&) = default;
};

using MaybeTourNodeHandle = std::optional<TourNodeHandle>;

// Owns the cyclic visit topology and every index that observes it. Mutations
// keep node storage, links, spatial indexes, and tree-node assignments in sync.
class Tour {
public:
    explicit Tour(std::size_t treeNodeCount);

    Tour(const Tour&) = delete;
    Tour& operator=(const Tour&) = delete;
    Tour(Tour&&) = delete;
    Tour& operator=(Tour&&) = delete;

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

    [[nodiscard]] MaybeTourNodeHandle visitFor(
        TreeNodeId treeNodeId) const;
    [[nodiscard]] Point point(TourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle previous(TourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle next(TourNodeHandle handle) const;

    [[nodiscard]] MaybeTourNodeHandle nearestVisit(
        const Point& point) const;
    [[nodiscard]] std::vector<TourNodeHandle> nearestEdgeStarts(
        const Point& point,
        std::size_t maximumCount) const;

    [[nodiscard]] TourNodeHandle createFirstVisit(
        const Point& point,
        TreeNodeId treeNodeId,
        std::size_t initialEnergy);
    [[nodiscard]] TourNodeHandle insertVisitBetween(
        const Point& point,
        TreeNodeId treeNodeId,
        TourNodeHandle previous,
        TourNodeHandle following,
        std::size_t initialEnergy);
    void addAssignment(
        TourNodeHandle handle,
        TreeNodeId treeNodeId,
        std::size_t energyIncrease);
    void removeAssignment(TreeNodeId treeNodeId);
    [[nodiscard]] std::vector<TreeNodeId> eraseVisit(
        TourNodeHandle handle);

    // Returns true when the visit has no energy remaining.
    [[nodiscard]] bool consumeEnergy(TourNodeHandle handle);
    [[nodiscard]] std::size_t recordInsertion(TourNodeHandle handle);
    [[nodiscard]] std::size_t insertionCount(
        TourNodeHandle handle) const;

    // Optimizes one visit while preserving all structural and spatial indexes.
    void optimizeVisit(
        TourNodeHandle handle,
        const std::vector<TreeNode>& treeNodes);

    [[nodiscard]] std::vector<Point> points() const;
    void assertValid() const;

private:
    using Segment = boost::geometry::model::segment<Point>;
    using PointValue = std::pair<Point, TourNodeHandle>;
    using SegmentValue = std::pair<Segment, TourNodeHandle>;
    using PointIndex = boost::geometry::index::rtree<
        PointValue,
        boost::geometry::index::rstar<16>>;
    using SegmentIndex = boost::geometry::index::rtree<
        SegmentValue,
        boost::geometry::index::rstar<16>>;

    struct Node {
        Point pos;
        IdSet<TreeNodeId> assignedTreeNodes;
        std::size_t energy;
        std::size_t insertions;
        TourNodeHandle prev;
        TourNodeHandle next;

        Node(
            const Point& point,
            TreeNodeId treeNodeId,
            std::size_t initialEnergy,
            TourNodeHandle previous,
            TourNodeHandle following);
    };

    struct NodeSlot {
        std::optional<Node> node;
        std::size_t generation = 1;
    };

    [[nodiscard]] Node* resolve(TourNodeHandle handle) noexcept;
    [[nodiscard]] const Node* resolve(TourNodeHandle handle) const noexcept;
    [[nodiscard]] Node& requireNode(TourNodeHandle handle);
    [[nodiscard]] const Node& requireNode(TourNodeHandle handle) const;
    void requireUnassignedTreeNode(TreeNodeId treeNodeId) const;
    [[nodiscard]] std::string formatHandle(
        MaybeTourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle createNode(
        const Point& point,
        TreeNodeId treeNodeId,
        std::size_t initialEnergy,
        MaybeTourNodeHandle previous = std::nullopt,
        MaybeTourNodeHandle following = std::nullopt);
    void destroyNode(TourNodeHandle handle);
    void deleteNode(TourNodeHandle handle);

    // Declared before observers and indexes so it is destroyed after them.
    std::vector<NodeSlot> nodeSlots;
    std::vector<std::size_t> freeSlots;
    std::size_t liveNodeCount = 0;
    MaybeTourNodeHandle head;
    PointIndex pointIndex;
    SegmentIndex segmentIndex;
    std::vector<MaybeTourNodeHandle> tourNodeForTreeNode;
};
