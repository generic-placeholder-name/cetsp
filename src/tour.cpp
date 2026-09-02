#include "tour.hpp"

#include "debug.hpp"

#include <boost/geometry.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using NodeId = MergeTree::NodeId;

TourEdge::TourEdge(
    TourNodeHandle start,
    TourNodeHandle end,
    Point startPoint,
    Point endPoint)
    : start_(start),
      end_(end),
      startPoint_(std::move(startPoint)),
      endPoint_(std::move(endPoint)) {}

TourNodeHandle TourEdge::start() const noexcept {
    return start_;
}

TourNodeHandle TourEdge::end() const noexcept {
    return end_;
}

const Point& TourEdge::startPoint() const noexcept {
    return startPoint_;
}

const Point& TourEdge::endPoint() const noexcept {
    return endPoint_;
}

Tour::Node::Node(
    const Point& point,
    NodeId treeNodeId,
    std::size_t initialEnergy,
    TourNodeHandle previous,
    TourNodeHandle following)
    : pos(point),
      energy(initialEnergy),
      insertions(1),
      prev(previous),
      next(following) {
    assignedTreeNodes.insert(treeNodeId);
}

Tour::Node* Tour::resolve(TourNodeHandle handle) noexcept {
    if (handle.slot_ >= nodeSlots.size()) {
        return nullptr;
    }
    NodeSlot& slot = nodeSlots[handle.slot_];
    if (slot.generation != handle.generation_) {
        return nullptr;
    }
    return slot.node ? &*slot.node : nullptr;
}

const Tour::Node* Tour::resolve(
    TourNodeHandle handle) const noexcept {
    if (handle.slot_ >= nodeSlots.size()) {
        return nullptr;
    }
    const NodeSlot& slot = nodeSlots[handle.slot_];
    if (slot.generation != handle.generation_) {
        return nullptr;
    }
    return slot.node ? &*slot.node : nullptr;
}

Tour::Node& Tour::requireNode(TourNodeHandle handle) {
    if (Node* node = resolve(handle)) {
        return *node;
    }
    throw std::logic_error("stale tour-node handle");
}

const Tour::Node& Tour::requireNode(
    TourNodeHandle handle) const {
    if (const Node* node = resolve(handle)) {
        return *node;
    }
    throw std::logic_error("stale tour-node handle");
}

void Tour::requireUnassignedTreeNode(NodeId treeNodeId) const {
    if (treeNodeId >= tourNodeForTreeNode.size()) {
        throw std::out_of_range("tree-node ID is outside this tour");
    }
    if (tourNodeForTreeNode[treeNodeId]) {
        throw std::logic_error("tree node already has a tour assignment");
    }
}

std::string Tour::formatHandle(MaybeTourNodeHandle handle) const {
    if (!handle) {
        return "<none>";
    }
    if (resolve(*handle) == nullptr) {
        return "<stale>";
    }
    return std::to_string(handle->slot_) + ":" +
           std::to_string(handle->generation_);
}

TourNodeHandle Tour::createNode(
    const Point& point,
    NodeId treeNodeId,
    std::size_t initialEnergy,
    MaybeTourNodeHandle previous,
    MaybeTourNodeHandle following) {
    if (previous.has_value() != following.has_value()) {
        throw std::logic_error(
            "tour node must have either two neighbors or neither");
    }
    if (previous &&
        (resolve(*previous) == nullptr || resolve(*following) == nullptr)) {
        throw std::logic_error("tour node neighbor handle is stale");
    }

    const bool reuseSlot = !freeSlots.empty();
    std::size_t slotIndex = 0;
    std::size_t generation = 1;
    if (reuseSlot) {
        slotIndex = freeSlots.back();
        const NodeSlot& slot = nodeSlots[slotIndex];
        assert(!slot.node);
        assert(slot.generation < std::numeric_limits<std::size_t>::max());
        generation = slot.generation;
    } else {
        slotIndex = nodeSlots.size();
    }

    const TourNodeHandle handle{slotIndex, generation};
    const TourNodeHandle previousHandle = previous.value_or(handle);
    const TourNodeHandle followingHandle = following.value_or(handle);
    std::optional<Node> owner{
        std::in_place,
        point,
        treeNodeId,
        initialEnergy,
        previousHandle,
        followingHandle};

    if (reuseSlot) {
        freeSlots.pop_back();
        nodeSlots[slotIndex].node = std::move(owner);
    } else {
        nodeSlots.push_back(NodeSlot{std::move(owner), generation});
    }
    ++liveNodeCount;
    return handle;
}

void Tour::destroyNode(TourNodeHandle handle) {
    if (handle.slot_ >= nodeSlots.size()) {
        throw std::logic_error("attempted to destroy an unknown tour node");
    }
    NodeSlot& slot = nodeSlots[handle.slot_];
    if (!slot.node || slot.generation != handle.generation_) {
        throw std::logic_error("attempted to destroy a stale tour node");
    }

    slot.node.reset();
    --liveNodeCount;

    // Generation values never wrap; exhausted slots remain retired.
    if (slot.generation == std::numeric_limits<std::size_t>::max() - 1) {
        slot.generation = std::numeric_limits<std::size_t>::max();
        return;
    }
    ++slot.generation;
    freeSlots.push_back(handle.slot_);
}

void Tour::deleteNode(TourNodeHandle handle) {
    Node& node = requireNode(handle);
    const TourNodeHandle leftHandle = node.prev;
    const TourNodeHandle rightHandle = node.next;
    if ((leftHandle == handle) != (rightHandle == handle)) {
        throw std::logic_error("tour node has a half-self-linked ring");
    }
    Node& left = requireNode(leftHandle);
    Node& right = requireNode(rightHandle);

    DBG("Deleting TourNode handle=" << formatHandle(handle)
        << " at (" << bg::get<0>(node.pos) << ", "
        << bg::get<1>(node.pos) << ")");

    pointIndex.remove({node.pos, handle});
    segmentIndex.remove({Segment(node.pos, right.pos), handle});

    if (rightHandle != handle) {
        left.next = rightHandle;
        right.prev = leftHandle;

        segmentIndex.remove({Segment(left.pos, node.pos), leftHandle});
        segmentIndex.insert({Segment(left.pos, right.pos), leftHandle});

        DBG("Updated segment for left handle="
            << formatHandle(leftHandle) << " to ("
            << bg::get<0>(left.pos) << ", " << bg::get<1>(left.pos)
            << ") -> (" << bg::get<0>(right.pos) << ", "
            << bg::get<1>(right.pos) << ")");
    }

    if (head == handle) {
        head = rightHandle != handle
            ? MaybeTourNodeHandle{rightHandle}
            : std::nullopt;
        DBG("Head updated to handle=" << formatHandle(head));
    }

    for (NodeId treeNodeId : node.assignedTreeNodes) {
        tourNodeForTreeNode[treeNodeId].reset();
    }

    DBG("About to delete node handle=" << formatHandle(handle));
    DBG("pos=(" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");
    DBG("prev=" << formatHandle(node.prev)
        << " next=" << formatHandle(node.next));
    DBG("assignedTreeNodes.size() = " << node.assignedTreeNodes.size());

    destroyNode(handle);
    DBG("TourNode deleted.");
}

Tour::Tour(std::size_t treeNodeCount)
    : tourNodeForTreeNode(treeNodeCount) {
    nodeSlots.reserve(treeNodeCount);
    freeSlots.reserve(treeNodeCount);
}

bool Tour::empty() const noexcept {
    return liveNodeCount == 0;
}

std::size_t Tour::size() const noexcept {
    return liveNodeCount;
}

MaybeTourNodeHandle Tour::visitFor(NodeId treeNodeId) const {
    if (treeNodeId >= tourNodeForTreeNode.size()) {
        throw std::out_of_range("tree-node ID is outside the tour assignment map");
    }
    return tourNodeForTreeNode[treeNodeId];
}

Point Tour::point(TourNodeHandle handle) const {
    return requireNode(handle).pos;
}

TourNodeHandle Tour::previous(TourNodeHandle handle) const {
    return requireNode(handle).prev;
}

TourNodeHandle Tour::next(TourNodeHandle handle) const {
    return requireNode(handle).next;
}

MaybeTourNodeHandle Tour::nearestVisit(const Point& point) const {
    const auto candidate = pointIndex.qbegin(bgi::nearest(point, 1));
    if (candidate == pointIndex.qend()) {
        return std::nullopt;
    }
    return candidate->second;
}

std::vector<TourEdge> Tour::nearestEdges(
    const Point& point,
    std::size_t maximumCount) const {
    std::vector<TourEdge> edges;
    edges.reserve(maximumCount);
    for (auto candidate = segmentIndex.qbegin(
             bgi::nearest(point, maximumCount));
         candidate != segmentIndex.qend(); ++candidate) {
        const TourNodeHandle startHandle = candidate->second;
        const Node& start = requireNode(startHandle);
        const TourNodeHandle endHandle = start.next;
        const Node& end = requireNode(endHandle);
        edges.push_back(TourEdge(
            startHandle, endHandle, start.pos, end.pos));
    }
    return edges;
}

TourNodeHandle Tour::createFirstVisit(
    const Point& point,
    NodeId treeNodeId,
    std::size_t initialEnergy) {
    if (head) {
        throw std::logic_error("cannot create a first visit in a nonempty tour");
    }
    requireUnassignedTreeNode(treeNodeId);
    const TourNodeHandle handle = createNode(
        point, treeNodeId, initialEnergy);
    head = handle;
    tourNodeForTreeNode[treeNodeId] = handle;
    pointIndex.insert({point, handle});
    segmentIndex.insert({Segment(point, point), handle});
    return handle;
}

TourNodeHandle Tour::insertVisitBetween(
    const Point& point,
    NodeId treeNodeId,
    TourNodeHandle previousHandle,
    TourNodeHandle followingHandle,
    std::size_t initialEnergy) {
    requireUnassignedTreeNode(treeNodeId);
    const Node& previousBefore = requireNode(previousHandle);
    const Node& followingBefore = requireNode(followingHandle);
    if (previousBefore.next != followingHandle ||
        followingBefore.prev != previousHandle) {
        throw std::logic_error("tour insertion handles do not name one edge");
    }

    const TourNodeHandle handle = createNode(
        point,
        treeNodeId,
        initialEnergy,
        previousHandle,
        followingHandle);
    Node& previous = requireNode(previousHandle);
    Node& following = requireNode(followingHandle);
    Node& node = requireNode(handle);

    segmentIndex.remove(
        {Segment(previous.pos, following.pos), previousHandle});
    pointIndex.insert({point, handle});
    segmentIndex.insert(
        {Segment(previous.pos, node.pos), previousHandle});
    segmentIndex.insert(
        {Segment(node.pos, following.pos), handle});

    previous.next = handle;
    following.prev = handle;
    tourNodeForTreeNode[treeNodeId] = handle;
    return handle;
}

void Tour::addAssignment(
    TourNodeHandle handle,
    NodeId treeNodeId,
    std::size_t energyIncrease) {
    requireUnassignedTreeNode(treeNodeId);
    Node& node = requireNode(handle);
    node.assignedTreeNodes.insert(treeNodeId);
    node.energy += energyIncrease;
    ++node.insertions;
    tourNodeForTreeNode[treeNodeId] = handle;
}

void Tour::removeAssignment(NodeId treeNodeId) {
    const MaybeTourNodeHandle handle = visitFor(treeNodeId);
    if (!handle) {
        return;
    }

    Node& node = requireNode(*handle);
    node.assignedTreeNodes.erase(treeNodeId);
    tourNodeForTreeNode[treeNodeId].reset();
    if (node.assignedTreeNodes.empty()) {
        deleteNode(*handle);
    }
}

std::vector<NodeId> Tour::eraseVisit(TourNodeHandle handle) {
    const Node& node = requireNode(handle);
    std::vector<NodeId> assignments(
        node.assignedTreeNodes.begin(), node.assignedTreeNodes.end());
    deleteNode(handle);
    return assignments;
}

bool Tour::consumeEnergy(TourNodeHandle handle) {
    Node& node = requireNode(handle);
    if (node.energy > 0) {
        --node.energy;
    }
    return node.energy == 0;
}

std::size_t Tour::recordInsertion(TourNodeHandle handle) {
    return ++requireNode(handle).insertions;
}

std::size_t Tour::insertionCount(TourNodeHandle handle) const {
    return requireNode(handle).insertions;
}

void Tour::optimizeVisit(
    TourNodeHandle handle,
    const MergeTree& mergeTree) {
    Node& node = requireNode(handle);
    if (node.prev == handle) return;
    constexpr double EPS = 1e-12;

    Node& previous = requireNode(node.prev);
    Node& following = requireNode(node.next);

    DBG("[optimizeTourPoint] Optimizing handle="
        << formatHandle(handle) << " at (" << bg::get<0>(node.pos)
        << ", " << bg::get<1>(node.pos) << "), assignedTreeNodes="
        << node.assignedTreeNodes.size());

    pointIndex.remove({node.pos, handle});
    segmentIndex.remove(
        {Segment(previous.pos, node.pos), node.prev});
    segmentIndex.remove(
        {Segment(node.pos, following.pos), handle});

    DBG("[optimizeTourPoint] Handle " << formatHandle(handle)
        << " assigned tree nodes: ");
    for (NodeId treeNodeId : node.assignedTreeNodes) {
        DBG_NOENDL(treeNodeId << " ");
    }
    DBG("");
    DBG("[optimizeTourPoint] Handle " << formatHandle(handle)
        << " prev=" << formatHandle(node.prev)
        << " pos=(" << bg::get<0>(previous.pos) << ", "
        << bg::get<1>(previous.pos) << ")"
        << " next=" << formatHandle(node.next)
        << " pos=(" << bg::get<0>(following.pos) << ", "
        << bg::get<1>(following.pos) << ")");

    for (auto assignedIt = node.assignedTreeNodes.begin();
         assignedIt != node.assignedTreeNodes.end();) {
        const auto currentAssignment = assignedIt++;
        const NodeId treeNodeId = *currentAssignment;
        const Circle& neighborhood = mergeTree.neighborhood(treeNodeId);
        const auto candidate = pointIndex.qbegin(
            bgi::nearest(neighborhood.center, 1));
        if (candidate != pointIndex.qend()) {
            const auto& [candidatePoint, otherHandle] = *candidate;
            const double distance = bg::distance(
                candidatePoint, neighborhood.center);
            if (otherHandle != handle &&
                distance <= neighborhood.r) {
                Node& other = requireNode(otherHandle);
                node.assignedTreeNodes.erase(currentAssignment);
                other.assignedTreeNodes.insert(treeNodeId);
                tourNodeForTreeNode[treeNodeId] = otherHandle;
                DBG("[optimizeTourPoint] Moved tree node " << treeNodeId
                    << " from handle " << formatHandle(handle)
                    << " to handle " << formatHandle(otherHandle));
            }
        }
    }

    if (node.assignedTreeNodes.empty()) {
        DBG("[optimizeTourPoint] No assignments left in handle "
            << formatHandle(handle) << ", deleting node.");
        deleteNode(handle);
        return;
    }

    Point a = previous.pos;
    Point b = following.pos;
    double ax = bg::get<0>(a), ay = bg::get<1>(a);
    double bx = bg::get<0>(b), by = bg::get<1>(b);
    double dx = bx - ax, dy = by - ay;
    double T_min = 0.0, T_max = 1.0;

    if (bg::distance(a, b) < EPS) {
        for (NodeId treeNodeId : node.assignedTreeNodes) {
            const Circle& neighborhood = mergeTree.neighborhood(treeNodeId);
            const auto& center = neighborhood.center;
            const double radius = neighborhood.r;
            if (bg::distance(a, center) > radius) {
                T_min = 1.0;
                T_max = 0.0;
                break;
            }
        }
    } else {
        for (NodeId treeNodeId : node.assignedTreeNodes) {
            const Circle& neighborhood = mergeTree.neighborhood(treeNodeId);
            const auto& center = neighborhood.center;
            const double cx = bg::get<0>(center);
            const double cy = bg::get<1>(center);
            const double radius = neighborhood.r;
            const double fx = ax - cx;
            const double fy = ay - cy;
            const double quadraticA = dx * dx + dy * dy;
            const double quadraticB = 2 * (fx * dx + fy * dy);
            const double quadraticC =
                fx * fx + fy * fy - radius * radius;

            const double discriminant =
                quadraticB * quadraticB -
                4 * quadraticA * quadraticC;
            if (discriminant < 0) {
                T_min = 1.0;
                T_max = 0.0;
                break;
            }

            const double root = std::sqrt(discriminant);
            const double firstRoot =
                (-quadraticB - root) / (2 * quadraticA);
            const double secondRoot =
                (-quadraticB + root) / (2 * quadraticA);
            double low = std::min(firstRoot, secondRoot);
            double high = std::max(firstRoot, secondRoot);
            low = std::max(low, 0.0);
            high = std::min(high, 1.0);
            if (low > high) {
                T_min = 1.0;
                T_max = 0.0;
                break;
            }

            T_min = std::max(T_min, low);
            T_max = std::min(T_max, high);
            if (T_min > T_max) break;
        }
    }

    if (T_min <= T_max) {
        const double t = 0.5 * (T_min + T_max);
        node.pos = Point(ax + t * dx, ay + t * dy);
        DBG("[optimizeTourPoint] Found common-intersection point at t=" << t);
    } else {
        double gradientX = 0;
        double gradientY = 0;
        const double x = bg::get<0>(node.pos);
        const double y = bg::get<1>(node.pos);
        {
            double offsetX = bg::get<0>(a) - x;
            double offsetY = bg::get<1>(a) - y;
            double distance = std::hypot(offsetX, offsetY);
            if (distance > EPS) {
                gradientX += offsetX / distance;
                gradientY += offsetY / distance;
            }
            offsetX = bg::get<0>(b) - x;
            offsetY = bg::get<1>(b) - y;
            distance = std::hypot(offsetX, offsetY);
            if (distance > EPS) {
                gradientX += offsetX / distance;
                gradientY += offsetY / distance;
            }
        }
        double norm = std::hypot(gradientX, gradientY);
        if (norm < EPS) norm = 1.0;
        gradientX /= norm;
        gradientY /= norm;

        double maximumStep = std::numeric_limits<double>::infinity();
        for (NodeId treeNodeId : node.assignedTreeNodes) {
            const Circle& circle = mergeTree.neighborhood(treeNodeId);
            const double centerX = bg::get<0>(circle.center);
            const double centerY = bg::get<1>(circle.center);
            const double offsetX = x - centerX;
            const double offsetY = y - centerY;
            const double quadraticA =
                gradientX * gradientX + gradientY * gradientY;
            const double quadraticB =
                2 * (offsetX * gradientX + offsetY * gradientY);
            const double quadraticC =
                offsetX * offsetX + offsetY * offsetY - circle.r * circle.r;
            const double discriminant = std::max(
                quadraticB * quadraticB -
                    4 * quadraticA * quadraticC,
                0.0);
            const double firstRoot =
                (-quadraticB + std::sqrt(discriminant)) /
                (2 * quadraticA);
            const double secondRoot =
                (-quadraticB - std::sqrt(discriminant)) /
                (2 * quadraticA);
            const double step = std::max({firstRoot, secondRoot, 0.0});
            if (step < maximumStep) maximumStep = step;
        }
        if (maximumStep < EPS ||
            maximumStep == std::numeric_limits<double>::infinity()) {
            DBG("[optimizeTourPoint] No valid step found for handle "
                << formatHandle(handle)
                << ", skipping optimization.");
        } else {
            DBG("[optimizeTourPoint] Gradient step for handle "
                << formatHandle(handle)
                << ": gx=" << gradientX
                << " gy=" << gradientY
                << " step=" << maximumStep);
            node.pos = Point(
                x + gradientX * maximumStep,
                y + gradientY * maximumStep);
        }
    }

    pointIndex.insert({node.pos, handle});
    segmentIndex.insert(
        {Segment(node.pos, following.pos), handle});
    segmentIndex.insert(
        {Segment(previous.pos, node.pos), node.prev});

    DBG("[optimizeTourPoint] Optimization complete for handle="
        << formatHandle(handle) << " at (" << bg::get<0>(node.pos)
        << ", " << bg::get<1>(node.pos) << ")");
}

std::vector<Point> Tour::points() const {
    std::vector<Point> result;
    result.reserve(liveNodeCount);
    if (!head) {
        return result;
    }

    TourNodeHandle current = *head;
    do {
        const Node& node = requireNode(current);
        result.push_back(node.pos);
        current = node.next;
    } while (current != *head);
    return result;
}

void Tour::assertValid() const {
#ifndef NDEBUG
    assert(head.has_value() == (liveNodeCount != 0));
    assert(pointIndex.size() == liveNodeCount);
    assert(segmentIndex.size() == liveNodeCount);

    std::vector<bool> isFree(nodeSlots.size(), false);
    for (std::size_t slotIndex : freeSlots) {
        assert(slotIndex < nodeSlots.size());
        assert(!isFree[slotIndex]);
        isFree[slotIndex] = true;
        assert(!nodeSlots[slotIndex].node);
        assert(nodeSlots[slotIndex].generation <
               std::numeric_limits<std::size_t>::max());
    }

    std::size_t occupiedSlots = 0;
    for (std::size_t index = 0; index < nodeSlots.size(); ++index) {
        const NodeSlot& slot = nodeSlots[index];
        if (!slot.node) {
            assert(isFree[index] ||
                   slot.generation ==
                       std::numeric_limits<std::size_t>::max());
            continue;
        }

        ++occupiedSlots;
        assert(!isFree[index]);
        assert(slot.generation < std::numeric_limits<std::size_t>::max());
        const TourNodeHandle handle{index, slot.generation};
        const Node& node = *slot.node;
        assert(!node.assignedTreeNodes.empty());
        const Node* previous = resolve(node.prev);
        const Node* following = resolve(node.next);
        assert(previous != nullptr);
        assert(following != nullptr);
        assert(previous->next == handle);
        assert(following->prev == handle);
        for (NodeId treeNodeId : node.assignedTreeNodes) {
            assert(treeNodeId < tourNodeForTreeNode.size());
            assert(tourNodeForTreeNode[treeNodeId] == handle);
        }
    }
    assert(occupiedSlots == liveNodeCount);

    for (NodeId treeNodeId = 0;
         treeNodeId < tourNodeForTreeNode.size(); ++treeNodeId) {
        const MaybeTourNodeHandle handle = tourNodeForTreeNode[treeNodeId];
        if (handle) {
            const Node* node = resolve(*handle);
            assert(node != nullptr);
            assert(node->assignedTreeNodes.contains(treeNodeId));
        }
    }

    std::vector<bool> pointIndexed(nodeSlots.size(), false);
    for (const PointValue& value : pointIndex) {
        const auto& [indexedPoint, handle] = value;
        const Node* node = resolve(handle);
        assert(node != nullptr);
        assert(bg::equals(indexedPoint, node->pos));
        assert(!pointIndexed[handle.slot_]);
        pointIndexed[handle.slot_] = true;
    }

    std::vector<bool> segmentIndexed(nodeSlots.size(), false);
    for (const SegmentValue& value : segmentIndex) {
        const auto& [segment, handle] = value;
        const Node* node = resolve(handle);
        assert(node != nullptr);
        const Node* following = resolve(node->next);
        assert(following != nullptr);
        assert(bg::equals(segment, Segment(node->pos, following->pos)));
        assert(!segmentIndexed[handle.slot_]);
        segmentIndexed[handle.slot_] = true;
    }

    for (std::size_t index = 0; index < nodeSlots.size(); ++index) {
        if (nodeSlots[index].node) {
            assert(pointIndexed[index]);
            assert(segmentIndexed[index]);
        }
    }

    if (head) {
        TourNodeHandle current = *head;
        std::size_t visited = 0;
        do {
            const Node* node = resolve(current);
            assert(node != nullptr);
            ++visited;
            assert(visited <= liveNodeCount);
            current = node->next;
        } while (current != *head);
        assert(visited == liveNodeCount);
    }
#endif
}
