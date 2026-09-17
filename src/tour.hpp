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

class Tour;

class TourNodeHandle {
public:
    friend constexpr bool operator==(
        const TourNodeHandle&,
        const TourNodeHandle&) = default;

private:
    constexpr TourNodeHandle(
        std::size_t slot,
        std::size_t generation) noexcept
        : slot_(slot), generation_(generation) {}

    std::size_t slot_;
    std::size_t generation_;

    friend class Tour;
};

using MaybeTourNodeHandle = std::optional<TourNodeHandle>;

// A directed tour edge and its endpoint positions at query time. Mutating the
// tour may invalidate its handles and make the point snapshot stale.
class TourEdge {
public:
    [[nodiscard]] TourNodeHandle start() const noexcept;
    [[nodiscard]] TourNodeHandle end() const noexcept;
    [[nodiscard]] const Point& startPoint() const noexcept;
    [[nodiscard]] const Point& endPoint() const noexcept;

private:
    TourEdge(
        TourNodeHandle start,
        TourNodeHandle end,
        Point startPoint,
        Point endPoint);

    TourNodeHandle start_;
    TourNodeHandle end_;
    Point startPoint_;
    Point endPoint_;

    friend class Tour;
};

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
        MergeTree::NodeId treeNodeId) const;
    [[nodiscard]] Point point(TourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle previous(TourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle next(TourNodeHandle handle) const;

    [[nodiscard]] MaybeTourNodeHandle nearestVisit(
        const Point& point) const;
    [[nodiscard]] std::vector<TourEdge> nearestEdges(
        const Point& point,
        std::size_t maximumCount) const;

    [[nodiscard]] TourNodeHandle createFirstVisit(
        const Point& point,
        MergeTree::NodeId treeNodeId,
        std::size_t initialEnergy);
    [[nodiscard]] TourNodeHandle insertVisitBetween(
        const Point& point,
        MergeTree::NodeId treeNodeId,
        TourNodeHandle previous,
        TourNodeHandle following,
        std::size_t initialEnergy);
    void addAssignment(
        TourNodeHandle handle,
        MergeTree::NodeId treeNodeId,
        std::size_t energyIncrease);
    void removeAssignment(MergeTree::NodeId treeNodeId);
    [[nodiscard]] std::vector<MergeTree::NodeId> eraseVisit(
        TourNodeHandle handle);

    // Returns true when the visit has no energy remaining.
    [[nodiscard]] bool consumeEnergy(TourNodeHandle handle);
    [[nodiscard]] std::size_t recordInsertion(TourNodeHandle handle);
    [[nodiscard]] std::size_t insertionCount(
        TourNodeHandle handle) const;

    // Optimizes one visit while preserving all structural and spatial indexes.
    void optimizeVisit(
        TourNodeHandle handle,
        const MergeTree& mergeTree);

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
        IdSet<MergeTree::NodeId> assignedTreeNodes;
        std::size_t energy;
        std::size_t insertions;
        // Generation-checked links survive arena relocation and detect any
        // accidentally retained reference to a recycled slot.
        TourNodeHandle prev;
        TourNodeHandle next;

        Node(
            const Point& point,
            MergeTree::NodeId treeNodeId,
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
    void requireUnassignedTreeNode(MergeTree::NodeId treeNodeId) const;
    [[nodiscard]] std::string formatHandle(
        MaybeTourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle createNode(
        const Point& point,
        MergeTree::NodeId treeNodeId,
        std::size_t initialEnergy,
        MaybeTourNodeHandle previous = std::nullopt,
        MaybeTourNodeHandle following = std::nullopt);
    void destroyNode(TourNodeHandle handle);
    void deleteNode(TourNodeHandle handle);

    // Declared before observers and indexes so it is destroyed after them.
    std::vector<NodeSlot> nodeSlots_;
    std::vector<std::size_t> freeSlots_;
    std::size_t liveNodeCount_ = 0;
    MaybeTourNodeHandle head_;
    PointIndex pointIndex_;
    SegmentIndex segmentIndex_;
    std::vector<MaybeTourNodeHandle> tourNodeForTreeNode_;
};
