#include "reconstruct.hpp"

#include "debug.hpp"
#include "id_set.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {

using Segment = bg::model::segment<Point>;

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
using PointValue = std::pair<Point, TourNodeHandle>;
using SegmentValue = std::pair<Segment, TourNodeHandle>;

constexpr std::size_t energyPerInsertion = 3;

[[nodiscard]] constexpr bool reachedPowerOfTwoCheckpoint(
    std::size_t value) noexcept {
    return value >= 2 && (value & (value - 1)) == 0;
}

// Tour has sole ownership of every live node. Persistent relationships use
// generation-checked handles; raw pointers returned by resolve() are transient.
struct TourNode {
    Point pos;
    IdSet<TreeNodeId> assignedTreeNodes;
    std::size_t energy;
    std::size_t insertions;
    TourNodeHandle prev;
    TourNodeHandle next;

    TourNode(
        const Point& point,
        TreeNodeId treeNodeId,
        TourNodeHandle previous,
        TourNodeHandle following)
        : pos(point),
          energy(energyPerInsertion),
          insertions(1),
          prev(previous),
          next(following) {
        assignedTreeNodes.insert(treeNodeId);
    }
};

void sortByRadius(
    std::vector<TreeNodeId>& ids,
    const std::vector<TreeNode>& treeNodes) {
    std::sort(ids.begin(), ids.end(), [&](TreeNodeId a, TreeNodeId b) {
        return treeNodes[a].r < treeNodes[b].r;
    });
}

class Tour {
public:
    explicit Tour(const std::vector<TreeNode>& treeNodes);
    Tour(const Tour&) = delete;
    Tour& operator=(const Tour&) = delete;
    Tour(Tour&&) = delete;
    Tour& operator=(Tour&&) = delete;

    [[nodiscard]] std::vector<Point> reconstruct();

private:
    using PointIndex = bgi::rtree<PointValue, bgi::rstar<16>>;
    using SegmentIndex = bgi::rtree<SegmentValue, bgi::rstar<16>>;

    enum class InsertionPhase {
        begin,
        afterPreviousNeighbor,
        afterFollowingNeighbor,
    };

    struct InsertionStep {
        TreeNodeId treeNodeId;
        InsertionPhase phase;
    };

    struct NodeSlot {
        std::optional<TourNode> node;
        std::size_t generation = 1;
    };

    [[nodiscard]] TourNode* resolve(TourNodeHandle handle) noexcept;
    [[nodiscard]] const TourNode* resolve(
        TourNodeHandle handle) const noexcept;
    [[nodiscard]] TourNode& requireNode(TourNodeHandle handle);
    [[nodiscard]] std::string formatNodeHandle(
        MaybeTourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle createNode(
        const Point& point,
        TreeNodeId treeNodeId,
        MaybeTourNodeHandle previous = std::nullopt,
        MaybeTourNodeHandle following = std::nullopt);
    void destroyNode(TourNodeHandle handle);
    void deleteNode(TourNodeHandle handle);
    void optimizePoint(TourNodeHandle handle);
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
    void assertValid() const;

    const std::vector<TreeNode>& treeNodes_;
    // Declared before observers and indexes so it is destroyed after them.
    std::vector<NodeSlot> nodeSlots_;
    std::vector<std::size_t> freeSlots_;
    std::size_t liveNodeCount_ = 0;
    MaybeTourNodeHandle head_;
    PointIndex pointIndex_;
    SegmentIndex segmentIndex_;
    std::vector<MaybeTourNodeHandle> tourNodeForTreeNode_;
    std::vector<InsertionStep> insertionStack_;
};

// Given a circle (center, radius) and an edge [p1,p2], returns the point
// on or within the circle that minimizes insertion‐cost between p1 and p2.
// Note: neither p1 nor p2 should be inside the circle.
Point findOptimalPoint(const Point& center, double radius,
                       const Point& p1, const Point& p2)
{
    const double EPS = 1e-12;

    // If radius is (nearly) zero, just return center
    if (radius < EPS) {
        return center;
    }

    // Try projecting the center onto the segment p1→p2
    {
        double dx = bg::get<0>(p2) - bg::get<0>(p1);
        double dy = bg::get<1>(p2) - bg::get<1>(p1);
        double len2 = dx*dx + dy*dy;
        if (len2 > EPS) {
            double vx = bg::get<0>(center) - bg::get<0>(p1);
            double vy = bg::get<1>(center) - bg::get<1>(p1);
            double t = (dx*vx + dy*vy) / len2;
            if (t >= 0.0 && t <= 1.0) {
                Point proj{
                    bg::get<0>(p1) + t * dx,
                    bg::get<1>(p1) + t * dy
                };
                if (bg::distance(proj, center) <= radius) {
                    return proj;
                }
            }
        }
    }

    // Otherwise approximate the boundary reflection point with one Newton step (Alhazen’s problem)

    // compute initial half-angle guess phi = (alpha + beta) / 2
    double alpha = std::atan2(
        bg::get<1>(p1) - bg::get<1>(center),
        bg::get<0>(p1) - bg::get<0>(center)
    );
    double beta = std::atan2(
        bg::get<1>(p2) - bg::get<1>(center),
        bg::get<0>(p2) - bg::get<0>(center)
    );
    double ax = std::cos(alpha), ay = std::sin(alpha);
    double bx = std::cos(beta),  by = std::sin(beta);
    double mx = ax + bx,         my = ay + by;
    double phi = std::atan2(my, mx);

    // define f(phi) = (p1−X)·T − (p2−X)·T
    auto f = [&](double phi) {
        double x = bg::get<0>(center) + radius * std::cos(phi);
        double y = bg::get<1>(center) + radius * std::sin(phi);
        double tx = -std::sin(phi);
        double ty =  std::cos(phi);
        double dot1 = (bg::get<0>(p1) - x) * tx
                    + (bg::get<1>(p1) - y) * ty;
        double dot2 = (bg::get<0>(p2) - x) * tx
                    + (bg::get<1>(p2) - y) * ty;
        return dot1 - dot2;
    };

    // derivative f'(phi)
    auto fprime = [&](double phi) {
        double x = bg::get<0>(center) + radius * std::cos(phi);
        double y = bg::get<1>(center) + radius * std::sin(phi);
        double ux = bg::get<0>(p1) - x;
        double uy = bg::get<1>(p1) - y;
        double vx = bg::get<0>(p2) - x;
        double vy = bg::get<1>(p2) - y;
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        double d1 = -(ux * cphi + uy * sphi)
                    + radius * (ux * (-sphi) + uy * cphi);
        double d2 = -(vx * cphi + vy * sphi)
                    + radius * (vx * (-sphi) + vy * cphi);
        return d1 - d2;
    };

    // one Newton–Raphson iteration
    double value = f(phi);
    double deriv = fprime(phi);
    if (std::abs(deriv) > EPS) {
        double new_phi = phi - value / deriv;
        Point oldPoint{
            bg::get<0>(center) + radius * std::cos(phi),
            bg::get<1>(center) + radius * std::sin(phi)
        };
        Point newPoint{
            bg::get<0>(center) + radius * std::cos(new_phi),
            bg::get<1>(center) + radius * std::sin(new_phi)
        };

        if (bg::distance(p1, oldPoint) + bg::distance(p2, oldPoint) <
            bg::distance(p1, newPoint) + bg::distance(p2, newPoint)) {
            return oldPoint; // If the new point is worse, return the old one
        } else {
            return newPoint; // Otherwise, use the new point
        }
    } else {
        return Point(
            bg::get<0>(center) + radius * std::cos(phi),
            bg::get<1>(center) + radius * std::sin(phi)
        ); // If derivative is zero, return the point on the circle
    }
}

Tour::Tour(const std::vector<TreeNode>& treeNodes)
    : treeNodes_(treeNodes), tourNodeForTreeNode_(treeNodes.size()) {
    nodeSlots_.reserve(treeNodes.size());
    freeSlots_.reserve(treeNodes.size());
    insertionStack_.reserve(treeNodes.size());
}

TourNode* Tour::resolve(TourNodeHandle handle) noexcept {
    if (handle.slot >= nodeSlots_.size()) {
        return nullptr;
    }
    NodeSlot& slot = nodeSlots_[handle.slot];
    if (slot.generation != handle.generation) {
        return nullptr;
    }
    return slot.node ? &*slot.node : nullptr;
}

const TourNode* Tour::resolve(TourNodeHandle handle) const noexcept {
    if (handle.slot >= nodeSlots_.size()) {
        return nullptr;
    }
    const NodeSlot& slot = nodeSlots_[handle.slot];
    if (slot.generation != handle.generation) {
        return nullptr;
    }
    return slot.node ? &*slot.node : nullptr;
}

TourNode& Tour::requireNode(TourNodeHandle handle) {
    if (TourNode* node = resolve(handle)) {
        return *node;
    }
    throw std::logic_error("stale tour-node handle");
}

std::string Tour::formatNodeHandle(MaybeTourNodeHandle handle) const {
    if (!handle) {
        return "<none>";
    }
    const TourNode* node = resolve(*handle);
    if (node == nullptr) {
        return "<stale>";
    }
    return std::to_string(handle->slot) + ":" +
           std::to_string(handle->generation);
}

void Tour::assertValid() const {
#ifndef NDEBUG
    assert(head_.has_value() == (liveNodeCount_ != 0));
    assert(pointIndex_.size() == liveNodeCount_);
    assert(segmentIndex_.size() == liveNodeCount_);

    std::vector<bool> isFree(nodeSlots_.size(), false);
    for (std::size_t slotIndex : freeSlots_) {
        assert(slotIndex < nodeSlots_.size());
        assert(!isFree[slotIndex]);
        isFree[slotIndex] = true;
        assert(!nodeSlots_[slotIndex].node);
        assert(nodeSlots_[slotIndex].generation <
               std::numeric_limits<std::size_t>::max());
    }

    std::size_t occupiedSlots = 0;
    for (std::size_t index = 0; index < nodeSlots_.size(); ++index) {
        const NodeSlot& slot = nodeSlots_[index];
        if (!slot.node) {
            assert(isFree[index] ||
                   slot.generation ==
                       std::numeric_limits<std::size_t>::max());
            continue;
        }

        ++occupiedSlots;
        assert(!isFree[index]);
        assert(slot.generation <
               std::numeric_limits<std::size_t>::max());
        const TourNodeHandle handle{index, slot.generation};
        const TourNode& node = *slot.node;
        assert(!node.assignedTreeNodes.empty());
        const TourNode* previous = resolve(node.prev);
        const TourNode* following = resolve(node.next);
        assert(previous != nullptr);
        assert(following != nullptr);
        assert(previous->next == handle);
        assert(following->prev == handle);
        for (TreeNodeId treeNodeId : node.assignedTreeNodes) {
            assert(treeNodeId < tourNodeForTreeNode_.size());
            assert(tourNodeForTreeNode_[treeNodeId] == handle);
        }
    }
    assert(occupiedSlots == liveNodeCount_);

    for (TreeNodeId treeNodeId = 0;
         treeNodeId < tourNodeForTreeNode_.size(); ++treeNodeId) {
        const MaybeTourNodeHandle handle = tourNodeForTreeNode_[treeNodeId];
        if (handle) {
            const TourNode* const node = resolve(*handle);
            assert(node != nullptr);
            assert(node->assignedTreeNodes.contains(treeNodeId));
        }
    }

    std::vector<bool> pointIndexed(nodeSlots_.size(), false);
    for (const PointValue& value : pointIndex_) {
        const auto& [point, handle] = value;
        const TourNode* node = resolve(handle);
        assert(node != nullptr);
        assert(bg::equals(point, node->pos));
        assert(!pointIndexed[handle.slot]);
        pointIndexed[handle.slot] = true;
    }

    std::vector<bool> segmentIndexed(nodeSlots_.size(), false);
    for (const SegmentValue& value : segmentIndex_) {
        const auto& [segment, handle] = value;
        const TourNode* node = resolve(handle);
        assert(node != nullptr);
        const TourNode* following = resolve(node->next);
        assert(following != nullptr);
        assert(bg::equals(segment, Segment(node->pos, following->pos)));
        assert(!segmentIndexed[handle.slot]);
        segmentIndexed[handle.slot] = true;
    }

    for (std::size_t index = 0; index < nodeSlots_.size(); ++index) {
        if (nodeSlots_[index].node) {
            assert(pointIndexed[index]);
            assert(segmentIndexed[index]);
        }
    }

    if (head_) {
        TourNodeHandle current = *head_;
        std::size_t visited = 0;
        do {
            const TourNode* node = resolve(current);
            assert(node != nullptr);
            ++visited;
            assert(visited <= liveNodeCount_);
            current = node->next;
        } while (current != *head_);
        assert(visited == liveNodeCount_);
    }
#endif
}

TourNodeHandle Tour::createNode(
    const Point& point,
    TreeNodeId treeNodeId,
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

    const bool reuseSlot = !freeSlots_.empty();
    std::size_t slotIndex = 0;
    std::size_t generation = 1;
    if (reuseSlot) {
        slotIndex = freeSlots_.back();
        const NodeSlot& slot = nodeSlots_[slotIndex];
        assert(!slot.node);
        assert(slot.generation <
               std::numeric_limits<std::size_t>::max());
        generation = slot.generation;
    } else {
        slotIndex = nodeSlots_.size();
    }

    const TourNodeHandle handle{slotIndex, generation};
    const TourNodeHandle previousHandle = previous.value_or(handle);
    const TourNodeHandle followingHandle = following.value_or(handle);
    std::optional<TourNode> owner{
        std::in_place,
        point,
        treeNodeId,
        previousHandle,
        followingHandle};

    if (reuseSlot) {
        freeSlots_.pop_back();
        nodeSlots_[slotIndex].node = std::move(owner);
    } else {
        nodeSlots_.push_back(NodeSlot{std::move(owner), generation});
    }
    ++liveNodeCount_;
    return handle;
}

void Tour::destroyNode(TourNodeHandle handle) {
    if (handle.slot >= nodeSlots_.size()) {
        throw std::logic_error("attempted to destroy an unknown tour node");
    }
    NodeSlot& slot = nodeSlots_[handle.slot];
    if (!slot.node || slot.generation != handle.generation) {
        throw std::logic_error("attempted to destroy a stale tour node");
    }

    slot.node.reset();
    --liveNodeCount_;

    // Generation values never wrap; exhausted slots remain retired.
    if (slot.generation ==
        std::numeric_limits<std::size_t>::max() - 1) {
        slot.generation = std::numeric_limits<std::size_t>::max();
        return;
    }
    ++slot.generation;
    freeSlots_.push_back(handle.slot);
}

// Delete a TourNode from the cyclic list & any other data structures.
void Tour::deleteNode(TourNodeHandle handle) {
    TourNode& node = requireNode(handle);
    const TourNodeHandle leftHandle = node.prev;
    const TourNodeHandle rightHandle = node.next;
    if ((leftHandle == handle) != (rightHandle == handle)) {
        throw std::logic_error("tour node has a half-self-linked ring");
    }
    TourNode& left = requireNode(leftHandle);
    TourNode& right = requireNode(rightHandle);

    DBG("Deleting TourNode handle=" << formatNodeHandle(handle)
        << " at (" << bg::get<0>(node.pos) << ", "
        << bg::get<1>(node.pos) << ")");

    // Remove point and segment belonging to 'node'
    pointIndex_.remove({node.pos, handle});
    segmentIndex_.remove({Segment(node.pos, right.pos), handle});

    if (rightHandle != handle) { // Singleton nodes link to themselves.
        // Splice out 'node': link L->R
        left.next = rightHandle;
        right.prev = leftHandle;

        // Update L's segment in R-tree
        segmentIndex_.remove({Segment(left.pos, node.pos), leftHandle});
        segmentIndex_.insert({Segment(left.pos, right.pos), leftHandle});

        DBG("Updated segment for left handle="
            << formatNodeHandle(leftHandle) << " to ("
            << bg::get<0>(left.pos) << ", " << bg::get<1>(left.pos)
            << ") -> (" << bg::get<0>(right.pos) << ", "
            << bg::get<1>(right.pos) << ")");
    }

    // Adjust head if needed
    if (head_ == handle) {
        head_ = rightHandle != handle
            ? MaybeTourNodeHandle{rightHandle}
            : std::nullopt;
        DBG("Head updated to handle=" << formatNodeHandle(head_));
    }

    // Clear the tree-node assignments.
    for (TreeNodeId treeNodeId : node.assignedTreeNodes) {
        tourNodeForTreeNode_[treeNodeId].reset();
    }

    DBG("About to delete node handle=" << formatNodeHandle(handle));
    DBG("pos=(" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");
    DBG("prev=" << formatNodeHandle(node.prev)
        << " next=" << formatNodeHandle(node.next));
    DBG("assignedTreeNodes.size() = " << node.assignedTreeNodes.size());

    destroyNode(handle);
    DBG("TourNode deleted.");
}

// Optimize the location of a point in the tour, possibly moving its circles to other points
void Tour::optimizePoint(TourNodeHandle handle) {
    TourNode& node = requireNode(handle);
    if (node.prev == handle) return;
    constexpr double EPS = 1e-12;

    TourNode& previous = requireNode(node.prev);
    TourNode& following = requireNode(node.next);

    DBG("[optimizeTourPoint] Optimizing handle="
        << formatNodeHandle(handle) << " at (" << bg::get<0>(node.pos)
        << ", " << bg::get<1>(node.pos) << "), assignedTreeNodes="
        << node.assignedTreeNodes.size());

    // Remove from R-trees
    pointIndex_.remove({node.pos, handle});
    segmentIndex_.remove({Segment(previous.pos, node.pos), node.prev});
    segmentIndex_.remove({Segment(node.pos, following.pos), handle});

    DBG("[optimizeTourPoint] Handle " << formatNodeHandle(handle)
        << " assigned tree nodes: ");
    for (TreeNodeId cid : node.assignedTreeNodes) {
        DBG_NOENDL(cid << " ");
    }
    DBG("");
    DBG("[optimizeTourPoint] Handle " << formatNodeHandle(handle)
        << " prev=" << formatNodeHandle(node.prev)
        << " pos=(" << bg::get<0>(previous.pos) << ", "
        << bg::get<1>(previous.pos) << ")"
        << " next=" << formatNodeHandle(node.next)
        << " pos=(" << bg::get<0>(following.pos) << ", "
        << bg::get<1>(following.pos) << ")"
    );

    // Try to move circles to other points
    for (auto assignedIt = node.assignedTreeNodes.begin();
         assignedIt != node.assignedTreeNodes.end();) {
        const auto currentAssignment = assignedIt++;
        const TreeNodeId cid = *currentAssignment;
        // The source node is absent from rtreeP, so the nearest point is the
        // only point we need to test.
        const auto candidate = pointIndex_.qbegin(
            bgi::nearest(treeNodes_[cid].center, 1));
        if (candidate != pointIndex_.qend()) {
            const auto& [point, otherHandle] = *candidate;
            const double distance = bg::distance(point, treeNodes_[cid].center);
            if (otherHandle != handle &&
                distance <= treeNodes_[cid].r) {
                TourNode& other = requireNode(otherHandle);
                node.assignedTreeNodes.erase(currentAssignment);
                other.assignedTreeNodes.insert(cid);
                tourNodeForTreeNode_[cid] = otherHandle;
                DBG("[optimizeTourPoint] Moved tree node " << cid
                    << " from handle " << formatNodeHandle(handle)
                    << " to handle " << formatNodeHandle(otherHandle));
            }
        }
    }

    // If no assignments remain, delete the tour node and return.
    if (node.assignedTreeNodes.empty()) {
        DBG("[optimizeTourPoint] No assignments left in handle "
            << formatNodeHandle(handle) << ", deleting node.");
        deleteNode(handle);
        return;
    }

    // Get previous and next points
    Point a = previous.pos, b = following.pos;

    // Check if [a, b] intersects the intersection of all circles
    double ax = bg::get<0>(a), ay = bg::get<1>(a);
    double bx = bg::get<0>(b), by = bg::get<1>(b);
    double dx = bx - ax, dy = by - ay;

    // [Tmin, Tmax] ∩= each circle‐interval
    double T_min = 0.0, T_max = 1.0;

    if (bg::distance(a, b) < EPS) {
        // If a and b are the same point, the quadratic equation degenerates
        // Therefore, we have to check it separately
        for (TreeNodeId cid : node.assignedTreeNodes) {
            const auto& C = treeNodes_[cid].center;
            double r = treeNodes_[cid].r;

            // Check if the circle contains point a (or b, since a == b)
            if (bg::distance(a, C) > r)  {
                T_min = 1.0; T_max = 0.0; // No intersection
                break;
            }
        }
    } else {
        for (TreeNodeId cid : node.assignedTreeNodes) {
            const auto& C  = treeNodes_[cid].center;
            double  cx     = bg::get<0>(C),
                    cy     = bg::get<1>(C),
                    r      = treeNodes_[cid].r;

            // Quadratic: A t^2 + B t + C0 <= 0
            double fx = ax - cx, fy = ay - cy;
            double Aq = dx*dx + dy*dy;
            double Bq = 2*(fx*dx + fy*dy);
            double Cq = fx*fx + fy*fy - r*r;

            double disc = Bq*Bq - 4*Aq*Cq;
            if (disc < 0) {
                // No real intersection → segment AB misses this circle
                T_min = 1.0; T_max = 0.0;
                break;
            }

            double s  = std::sqrt(disc);
            double t1 = (-Bq - s) / (2*Aq);
            double t2 = (-Bq + s) / (2*Aq);
            double lo = std::min(t1, t2);
            double hi = std::max(t1, t2);

            // Intersect with [0,1]
            lo = std::max(lo, 0.0);
            hi = std::min(hi, 1.0);
            if (lo > hi) {
                T_min = 1.0; T_max = 0.0;
                break;
            }

            T_min = std::max(T_min, lo);
            T_max = std::min(T_max, hi);
            if (T_min > T_max) break;
        }
    }

    if (T_min <= T_max) {
        // There *is* a nonempty intersection
        double t_star = 0.5 * (T_min + T_max);
        Point P(ax + t_star * dx, ay + t_star * dy);
        DBG("[optimizeTourPoint] Found common-intersection point at t=" << t_star);
        node.pos = P;
    } else {
        // Projected gradient descent step
        // Gradient: sum of (negative) gradients of distances to a and b
        double gx = 0, gy = 0;
        double x = bg::get<0>(node.pos), y = bg::get<1>(node.pos);
        {
            double dx = bg::get<0>(a) - x, dy = bg::get<1>(a) - y;
            double d = std::hypot(dx, dy);
            if (d > EPS) { gx += dx / d; gy += dy / d; }
            dx = bg::get<0>(b) - x; dy = bg::get<1>(b) - y;
            d = std::hypot(dx, dy);
            if (d > EPS) { gx += dx / d; gy += dy / d; }
        }
        double norm = std::hypot(gx, gy);
        if (norm < EPS) norm = 1.0;
        gx /= norm; gy /= norm;

        // Find max step size before leaving any circle
        double maxStep = std::numeric_limits<double>::infinity();
        for (TreeNodeId cid : node.assignedTreeNodes) {
            const TreeNode& c = treeNodes_[cid];
            double cx = bg::get<0>(c.center), cy = bg::get<1>(c.center);
            double dx = x - cx, dy = y - cy;
            double a2 = gx*gx + gy*gy;
            double b2 = 2 * (dx*gx + dy*gy);
            double c2 = dx*dx + dy*dy - c.r * c.r;
            // Solve a2*t^2 + b2*t + c2 = 0 for t
            double disc = std::max(b2*b2 - 4*a2*c2, 0.);
            double t1 = (-b2 + std::sqrt(disc)) / (2*a2);
            double t2 = (-b2 - std::sqrt(disc)) / (2*a2);
            double tmax = std::max({t1, t2, 0.});
            if (tmax < maxStep) maxStep = tmax;
        }
        if (maxStep < EPS || maxStep == std::numeric_limits<double>::infinity()) {
            DBG("[optimizeTourPoint] No valid step found for handle "
                << formatNodeHandle(handle) << ", skipping optimization.");
        }
        else {
            DBG("[optimizeTourPoint] Gradient step for handle "
                << formatNodeHandle(handle) << ": gx=" << gx
                << " gy=" << gy << " step=" << maxStep);
            node.pos = Point(x + gx * maxStep, y + gy * maxStep);
        }
    }

    // Reinsert into R-trees
    pointIndex_.insert({node.pos, handle});
    segmentIndex_.insert({Segment(node.pos, following.pos), handle});
    segmentIndex_.insert({Segment(previous.pos, node.pos), node.prev});

    DBG("[optimizeTourPoint] Optimization complete for handle="
        << formatNodeHandle(handle) << " at (" << bg::get<0>(node.pos)
        << ", " << bg::get<1>(node.pos) << ")");
}

std::vector<TreeNodeId> Tour::processNeighbor(
    TourNodeHandle neighborHandle) {
    TourNode& neighbor = requireNode(neighborHandle);
    if (neighbor.energy > 0) {
        --neighbor.energy;
    }
    if (neighbor.energy != 0) {
        return {};
    }

    DBG("Neighbor handle=" << formatNodeHandle(neighborHandle)
        << " has zero energy, removing and reinserting its assignments.");
    std::vector<TreeNodeId> linked(
        neighbor.assignedTreeNodes.begin(), neighbor.assignedTreeNodes.end());
    deleteNode(neighborHandle);
    sortByRadius(linked, treeNodes_);
    return linked;
}

void Tour::scheduleAfterNeighbor(
    TreeNodeId treeNodeId,
    InsertionPhase continuation,
    TourNodeHandle neighborHandle) {
    std::vector<TreeNodeId> reinsertions = processNeighbor(neighborHandle);

    // This is an explicit encoding of the previous recursive call order. The
    // continuation sits below every reinsertion, and reverse scheduling makes
    // the first radius-sorted assignment execute first on the LIFO stack.
    insertionStack_.push_back({treeNodeId, continuation});
    for (auto it = reinsertions.rbegin(); it != reinsertions.rend(); ++it) {
        insertionStack_.push_back({*it, InsertionPhase::begin});
    }
}

// Insert a merge-tree circle into the cyclic tour, updating spatial indexes
// and triggering any energy-based delete/reinsert cascades.
void Tour::beginCircleInsertion(TreeNodeId cid) {
    const TreeNode& tn = treeNodes_[cid];
    const Point& center = tn.center;
    double rad = tn.r;

    DBG(pointIndex_.size() << " points in R-tree, "
        << segmentIndex_.size() << " segments in R-tree.");
    DBG("Inserting circle id=" << cid << " at (" << bg::get<0>(center) << ", " << bg::get<1>(center) << "), r=" << rad);

    // Empty tour: create first node
    if (!head_) {
        const TourNodeHandle nodeHandle = createNode(center, cid);
        head_ = nodeHandle;
        tourNodeForTreeNode_[cid] = nodeHandle;
        pointIndex_.insert({center, nodeHandle});
        // For a single point, use a degenerate segment (point, point)
        segmentIndex_.insert({Segment(center, center), nodeHandle});
        DBG("Created first TourNode handle="
            << formatNodeHandle(nodeHandle) << " for tree node " << cid);
        return;
    }

    // 1) Try to link to an existing point via R-tree
    DBG("Trying to link to an existing point via R-tree...");
    const auto nearestPoint = pointIndex_.qbegin(bgi::nearest(center, 1));
    if (nearestPoint != pointIndex_.qend()) {
        const auto& [point, existingHandle] = *nearestPoint;
        const double distance = bg::distance(point, center);
        if (distance <= rad) {
            TourNode& node = requireNode(existingHandle);
            node.assignedTreeNodes.insert(cid);
            node.energy += energyPerInsertion;
            node.insertions++;
            tourNodeForTreeNode_[cid] = existingHandle;

            DBG("Linked tree node " << cid << " to existing handle "
                << formatNodeHandle(existingHandle));

            scheduleAfterNeighbor(
                cid,
                InsertionPhase::afterPreviousNeighbor,
                node.prev);
            return;
        }
    }

    // 2) Find best edge to insert via segment R-tree
    DBG("Finding best edge to insert via segment R-tree...");
    constexpr int K = 16;
    double bestAdd = std::numeric_limits<double>::infinity();
    MaybeTourNodeHandle bestU;
    Point bestPt;

    for (auto candidate = segmentIndex_.qbegin(bgi::nearest(center, K));
         candidate != segmentIndex_.qend(); ++candidate) {
        const TourNodeHandle uHandle = candidate->second;
        const TourNode& u = requireNode(uHandle);
        const TourNode& v = requireNode(u.next);
        Point cand = findOptimalPoint(center, rad, u.pos, v.pos);
        double addCost = bg::distance(u.pos, cand) +
            bg::distance(cand, v.pos) - bg::distance(u.pos, v.pos);
        if (addCost < bestAdd) {
            bestAdd = addCost;
            bestU = uHandle;
            bestPt = cand;
        }
    }

    if (!bestU) {
        throw std::logic_error("failed to find a tour edge for circle insertion");
    }

    const TourNodeHandle leftHandle = *bestU;
    const TourNodeHandle rightHandle = requireNode(leftHandle).next;
    DBG("Best edge for insertion: between handle="
        << formatNodeHandle(leftHandle) << " and handle="
        << formatNodeHandle(rightHandle)
        << " at point (" << bg::get<0>(bestPt) << ", " << bg::get<1>(bestPt) << "), addCost=" << bestAdd);

    // 3) Insert new TourNode between bestU and bestU->next
    const TourNodeHandle nodeHandle = createNode(
        bestPt, cid, leftHandle, rightHandle);
    TourNode& left = requireNode(leftHandle);
    TourNode& right = requireNode(rightHandle);
    TourNode& node = requireNode(nodeHandle);

    // Remove old segment and insert new segments in R-tree
    segmentIndex_.remove({Segment(left.pos, right.pos), leftHandle});
    pointIndex_.insert({bestPt, nodeHandle});
    segmentIndex_.insert({Segment(left.pos, node.pos), leftHandle});
    segmentIndex_.insert({Segment(node.pos, right.pos), nodeHandle});

    // Link in the new node
    left.next = nodeHandle;
    right.prev = nodeHandle;

    tourNodeForTreeNode_[cid] = nodeHandle;

    DBG("Inserted new TourNode handle=" << formatNodeHandle(nodeHandle)
        << " for tree node " << cid << " between handles "
        << formatNodeHandle(leftHandle) << " and "
        << formatNodeHandle(rightHandle));

    scheduleAfterNeighbor(
        cid,
        InsertionPhase::afterPreviousNeighbor,
        node.prev);
}

void Tour::resumeCircleInsertionAfterPreviousNeighbor(
    TreeNodeId treeNodeId) {
    const MaybeTourNodeHandle currentHandle =
        tourNodeForTreeNode_[treeNodeId];
    if (!currentHandle) {
        throw std::logic_error(
            "tree node lost its tour assignment during previous-neighbor processing");
    }

    scheduleAfterNeighbor(
        treeNodeId,
        InsertionPhase::afterFollowingNeighbor,
        requireNode(*currentHandle).next);
}

void Tour::finishCircleInsertion(TreeNodeId treeNodeId) {
    const MaybeTourNodeHandle currentHandle =
        tourNodeForTreeNode_[treeNodeId];
    if (!currentHandle) {
        throw std::logic_error(
            "tree node lost its tour assignment during following-neighbor processing");
    }

    const TourNode& updatedNode = requireNode(*currentHandle);
    if (reachedPowerOfTwoCheckpoint(updatedNode.insertions)) {
        optimizePoint(*currentHandle);
    }
}

void Tour::insertCircle(TreeNodeId treeNodeId) {
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

std::vector<Point> Tour::reconstruct() {
    const std::size_t N = treeNodes_.size();

    if (treeNodes_.empty()) {
        return {};
    }

    DBG("Starting unmerge with " << N << " tree nodes.");

    // Priority queue to process nodes in descending order of merge gap.
    std::priority_queue<std::pair<double, TreeNodeId>> pq;
    pq.push({treeNodes_[N - 1].mergeGap, N - 1}); // Start from the root

    // Insert the root node's circle into the tour
    DBG("Inserting root node id=" << (N - 1));
    insertCircle(N - 1);
    assertValid();

    // Process nodes in order of decreasing merge gap.
    size_t nodesProcessed = 0;
    while (!pq.empty()) {
        auto [mergeGap, idx] = pq.top();
        pq.pop();
        nodesProcessed++;
        if (treeNodes_[idx].isLeaf()) continue;

        DBG("Unmerging node id=" << idx << " (mergeGap=" << mergeGap
            << "), left=" << treeNodes_[idx].left
        << ", right=" << treeNodes_[idx].right);

        // Remove the current node's circle from the tour
        const MaybeTourNodeHandle tourNodeHandle = tourNodeForTreeNode_[idx];
        if (tourNodeHandle) {
            TourNode& tourNode = requireNode(*tourNodeHandle);
            auto& assignments = tourNode.assignedTreeNodes;
            assignments.erase(idx);
            tourNodeForTreeNode_[idx].reset();

            if (assignments.empty()) {
                DBG("Deleting TourNode id=" << idx << " during unmerge.");
                deleteNode(*tourNodeHandle);
            }
        }
        // Insert the left and right children into the tour
        for (TreeNodeId c : treeNodes_[idx].children()) {
            DBG("Inserting child node id=" << c);
            insertCircle(c);
            // If the child is not a leaf, add it to the queue for further unmerging
            if (!treeNodes_[c].isLeaf()) {
                pq.push({treeNodes_[c].mergeGap, c});
            }
        }

        // Maintenance is weighted by assigned-tree-node multiplicity. Scanning
        // tourNodeForTreeNode_ visits a tour node once per tree node assigned
        // to it, so
        // shared nodes accumulate optimization visits and consume energy faster.
        if (reachedPowerOfTwoCheckpoint(nodesProcessed) || pq.empty()) {
            // Optimize tour nodes with multiplicity determined by assignment count.
            for (TreeNodeId nid = N; nid-- > 0;) {
                const MaybeTourNodeHandle handle = tourNodeForTreeNode_[nid];
                if (!handle) continue;
                TourNode& node = requireNode(*handle);
                node.insertions++;
                if (reachedPowerOfTwoCheckpoint(node.insertions)) {
                    DBG("Optimizing TourNode handle="
                        << formatNodeHandle(*handle) << " at ("
                        << bg::get<0>(node.pos) << ", "
                        << bg::get<1>(node.pos) << ")");
                    optimizePoint(*handle);
                }
            }

            // Reduce energy once per associated circle ID and rebuild nodes
            // whose energy reaches zero.
            for (TreeNodeId nid = N; nid-- > 0;) {
                DBG("Processing node id=" << nid << " for energy reduction.");
                const MaybeTourNodeHandle handle = tourNodeForTreeNode_[nid];
                if (!handle) continue;
                TourNode& node = requireNode(*handle);
                if (node.energy > 0) {
                    --node.energy;
                }
                if (node.energy == 0) {
                    std::vector<TreeNodeId> linked(
                        node.assignedTreeNodes.begin(),
                        node.assignedTreeNodes.end());
                    deleteNode(*handle);
                    sortByRadius(linked, treeNodes_);
                    for (TreeNodeId treeNodeId : linked) {
                        insertCircle(treeNodeId);
                    }
                }
            }

            // Repeat the same assignment-weighted optimization after rebuilding.
            for (TreeNodeId nid = N; nid-- > 0;) {
                const MaybeTourNodeHandle handle = tourNodeForTreeNode_[nid];
                if (!handle) continue;
                TourNode& node = requireNode(*handle);
                node.insertions++;
                if (reachedPowerOfTwoCheckpoint(node.insertions)) {
                    DBG("Optimizing TourNode handle="
                        << formatNodeHandle(*handle) << " at ("
                        << bg::get<0>(node.pos) << ", "
                        << bg::get<1>(node.pos) << ")");
                    optimizePoint(*handle);
                }
            }
        }

        // Full-tour invariant check at the phase boundary.
        assertValid();
    }

    // Collect the tour points in order, starting from head
    std::vector<Point> tour;
    if (head_) {
        TourNodeHandle current = *head_;
        do {
            const TourNode& node = requireNode(current);
            tour.push_back(node.pos);
            current = node.next;
        } while (current != *head_);
    }

    DBG("Tour reconstruction complete. Tour size: " << tour.size());
    assertValid();

    return tour;
}

} // namespace

std::vector<Point> reconstructTour(const std::vector<TreeNode>& treeNodes) {
    Tour tour(treeNodes);
    return tour.reconstruct();
}
