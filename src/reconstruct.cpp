#include "reconstruct.hpp"

#include "debug.hpp"
#include "hash_set.hpp"

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
using PointValue = std::pair<Point, TreeNodeId>;
using SegmentValue = std::pair<Segment, TreeNodeId>;

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

// Tour has sole ownership of every live node. Persistent relationships use
// generation-checked handles; raw pointers returned by resolve() are transient.
struct TourNode {
    Point pos;
    HashSet<TreeNodeId> circles;
    std::size_t energy;
    std::size_t insertions;
    TourNodeHandle prev;
    TourNodeHandle next;
    TreeNodeId id;

    TourNode(
        const Point& point,
        TreeNodeId circleId,
        TourNodeHandle previous,
        TourNodeHandle following)
        : pos(point),
          energy(3),
          insertions(1),
          prev(previous),
          next(following),
          id(circleId) {
        circles.insert(circleId);
    }
};

void sortByRadius(
    std::vector<size_t>& ids,
    const std::vector<TreeNode>& treeNodes) {
    std::sort(ids.begin(), ids.end(), [&](size_t a, size_t b) {
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

    struct NodeSlot {
        std::optional<TourNode> node;
        std::size_t generation = 1;
    };

    [[nodiscard]] TourNode* resolve(TourNodeHandle handle) noexcept;
    [[nodiscard]] const TourNode* resolve(
        TourNodeHandle handle) const noexcept;
    [[nodiscard]] TourNode& requireNode(TourNodeHandle handle);
    [[nodiscard]] std::string formatNodeId(
        MaybeTourNodeHandle handle) const;
    [[nodiscard]] TourNodeHandle createNode(
        const Point& point,
        TreeNodeId circleId,
        MaybeTourNodeHandle previous = std::nullopt,
        MaybeTourNodeHandle following = std::nullopt);
    void destroyNode(TourNodeHandle handle);
    void deleteNode(TourNodeHandle handle);
    void optimizePoint(TourNodeHandle handle);
    void processNeighbor(TourNodeHandle neighborHandle);
    void insertCircle(TreeNodeId circleId);
    void assertValid() const;

    const std::vector<TreeNode>& treeNodes_;
    // Declared before observers and indexes so it is destroyed after them.
    std::vector<NodeSlot> nodeSlots_;
    std::vector<std::size_t> freeSlots_;
    std::size_t liveNodeCount_ = 0;
    MaybeTourNodeHandle head_;
    PointIndex pointIndex_;
    SegmentIndex segmentIndex_;
    std::vector<MaybeTourNodeHandle> nodeForCircle_;
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
    : treeNodes_(treeNodes), nodeForCircle_(treeNodes.size()) {
    nodeSlots_.reserve(treeNodes.size());
    freeSlots_.reserve(treeNodes.size());
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

std::string Tour::formatNodeId(MaybeTourNodeHandle handle) const {
    if (!handle) {
        return "<none>";
    }
    const TourNode* node = resolve(*handle);
    return node == nullptr ? "<stale>" : std::to_string(node->id);
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
        assert(!node.circles.empty());
        assert(node.circles.contains(node.id));
        const TourNode* previous = resolve(node.prev);
        const TourNode* following = resolve(node.next);
        assert(previous != nullptr);
        assert(following != nullptr);
        assert(previous->next == handle);
        assert(following->prev == handle);
        for (TreeNodeId circleId : node.circles) {
            assert(circleId < nodeForCircle_.size());
            assert(nodeForCircle_[circleId] == handle);
        }
    }
    assert(occupiedSlots == liveNodeCount_);

    for (TreeNodeId circleId = 0;
         circleId < nodeForCircle_.size(); ++circleId) {
        const MaybeTourNodeHandle handle = nodeForCircle_[circleId];
        if (handle) {
            const TourNode* const node = resolve(*handle);
            assert(node != nullptr);
            assert(node->circles.contains(circleId));
        }
    }

    std::vector<bool> pointIndexed(nodeSlots_.size(), false);
    for (const PointValue& value : pointIndex_) {
        const auto& [point, representativeId] = value;
        assert(representativeId < nodeForCircle_.size());
        const MaybeTourNodeHandle handle =
            nodeForCircle_[representativeId];
        assert(handle.has_value());
        const TourNode* node = resolve(*handle);
        assert(node != nullptr);
        assert(node->id == representativeId);
        assert(bg::equals(point, node->pos));
        assert(!pointIndexed[handle->slot]);
        pointIndexed[handle->slot] = true;
    }

    std::vector<bool> segmentIndexed(nodeSlots_.size(), false);
    for (const SegmentValue& value : segmentIndex_) {
        const auto& [segment, representativeId] = value;
        assert(representativeId < nodeForCircle_.size());
        const MaybeTourNodeHandle handle =
            nodeForCircle_[representativeId];
        assert(handle.has_value());
        const TourNode* node = resolve(*handle);
        assert(node != nullptr);
        assert(node->id == representativeId);
        const TourNode* following = resolve(node->next);
        assert(following != nullptr);
        assert(bg::equals(segment, Segment(node->pos, following->pos)));
        assert(!segmentIndexed[handle->slot]);
        segmentIndexed[handle->slot] = true;
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
    TreeNodeId circleId,
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
        circleId,
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

    DBG("Deleting TourNode id=" << node.id << " at (" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");

    // Remove point and segment belonging to 'node'
    pointIndex_.remove({node.pos, node.id});
    segmentIndex_.remove({Segment(node.pos, right.pos), node.id});

    if (rightHandle != handle) { // Singleton nodes link to themselves.
        // Splice out 'node': link L->R
        left.next = rightHandle;
        right.prev = leftHandle;

        // Update L's segment in R-tree
        segmentIndex_.remove({Segment(left.pos, node.pos), left.id});
        segmentIndex_.insert({Segment(left.pos, right.pos), left.id});

        DBG("Updated segment for L id=" << left.id << " to (" << bg::get<0>(left.pos) << ", " << bg::get<1>(left.pos) << ") -> (" << bg::get<0>(right.pos) << ", " << bg::get<1>(right.pos) << ")");
    }

    // Adjust head if needed
    if (head_ == handle) {
        head_ = rightHandle != handle
            ? MaybeTourNodeHandle{rightHandle}
            : std::nullopt;
        DBG("Head updated to id=" << formatNodeId(head_));
    }

    // Clear the circles list and set map[cid] to nullptr
    for (TreeNodeId cid : node.circles) {
        nodeForCircle_[cid].reset();
    }

    DBG("About to delete node id=" << node.id);
    DBG("pos=(" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");
    DBG("prev=" << formatNodeId(node.prev)
        << " next=" << formatNodeId(node.next));
    DBG("circles.size() = " << node.circles.size());

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

    DBG("[optimizeTourPoint] Optimizing node id=" << node.id << " at (" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << "), circles=" << node.circles.size());

    // Remove from R-trees
    pointIndex_.remove({node.pos, node.id});
    segmentIndex_.remove({Segment(previous.pos, node.pos), previous.id});
    segmentIndex_.remove({Segment(node.pos, following.pos), node.id});

    DBG("[optimizeTourPoint] Node " << node.id << " circles: ");
    for (size_t cid : node.circles) {
        DBG_NOENDL(cid << " ");
    }
    DBG("");
    DBG("[optimizeTourPoint] Node " << node.id
        << " prev: id=" << formatNodeId(node.prev)
        << " pos=(" << bg::get<0>(previous.pos) << ", "
        << bg::get<1>(previous.pos) << ")"
        << " next: id=" << formatNodeId(node.next)
        << " pos=(" << bg::get<0>(following.pos) << ", "
        << bg::get<1>(following.pos) << ")"
    );

    // Try to move circles to other points
    for (auto circleIt = node.circles.begin();
         circleIt != node.circles.end();) {
        const auto currentCircle = circleIt++;
        const TreeNodeId cid = *currentCircle;
        // The source node is absent from rtreeP, so the nearest point is the
        // only point we need to test.
        const auto candidate = pointIndex_.qbegin(
            bgi::nearest(treeNodes_[cid].center, 1));
        if (candidate != pointIndex_.qend()) {
            const auto& [point, pointId] = *candidate;
            const MaybeTourNodeHandle otherHandle = nodeForCircle_[pointId];
            const double distance = bg::distance(point, treeNodes_[cid].center);
            if (otherHandle && *otherHandle != handle &&
                distance <= treeNodes_[cid].r) {
                TourNode& other = requireNode(*otherHandle);
                node.circles.erase(currentCircle);
                other.circles.insert(cid);
                nodeForCircle_[cid] = *otherHandle;
                DBG("[optimizeTourPoint] Moved circle " << cid << " from node " << node.id << " to node " << other.id);
            }
        }
    }

    // If no circles left, delete node and return
    if (node.circles.empty()) {
        DBG("[optimizeTourPoint] No circles left in node " << node.id << ", deleting node.");
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
        for (TreeNodeId cid : node.circles) {
            const auto& C = treeNodes_[cid].center;
            double r = treeNodes_[cid].r;

            // Check if the circle contains point a (or b, since a == b)
            if (bg::distance(a, C) > r)  {
                T_min = 1.0; T_max = 0.0; // No intersection
                break;
            }
        }
    } else {
        for (TreeNodeId cid : node.circles) {
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
        for (TreeNodeId cid : node.circles) {
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
            DBG("[optimizeTourPoint] No valid step found for node " << node.id << ", skipping optimization.");
        }
        else {
            DBG("[optimizeTourPoint] Gradient step for node " << node.id << ": gx=" << gx << " gy=" << gy << " step=" << maxStep);
            node.pos = Point(x + gx * maxStep, y + gy * maxStep);
        }
    }

    // If id is not in the circle anymore, reset it
    if (!node.circles.contains(node.id)) {
        TreeNodeId oldId = node.id;
        node.id = *node.circles.begin();
        DBG("[optimizeTourPoint] Changed node id from " << oldId << " to " << node.id);
    }

    // Reinsert into R-trees
    pointIndex_.insert({node.pos, node.id});
    segmentIndex_.insert({Segment(node.pos, following.pos), node.id});
    segmentIndex_.insert({Segment(previous.pos, node.pos), previous.id});

    DBG("[optimizeTourPoint] Optimization complete for node id=" << node.id << " at (" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");
}

void Tour::processNeighbor(TourNodeHandle neighborHandle) {
    TourNode& neighbor = requireNode(neighborHandle);
    if (neighbor.energy > 0) {
        --neighbor.energy;
    }
    if (neighbor.energy != 0) {
        return;
    }

    DBG("Neighbor id=" << neighbor.id
        << " has zero energy, removing and reinserting its circles.");
    std::vector<TreeNodeId> linked(
        neighbor.circles.begin(), neighbor.circles.end());
    deleteNode(neighborHandle);
    sortByRadius(linked, treeNodes_);
    for (TreeNodeId circleId : linked) {
        insertCircle(circleId);
    }
}

// Insert a merge-tree circle into the cyclic tour, updating spatial indexes
// and triggering any energy-based delete/reinsert cascades.
void Tour::insertCircle(TreeNodeId cid) {
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
        nodeForCircle_[cid] = nodeHandle;
        pointIndex_.insert({center, cid});
        // For a single point, use a degenerate segment (point, point)
        segmentIndex_.insert({Segment(center, center), cid});
        DBG("Created first TourNode id=" << cid);
        return;
    }

    // 1) Try to link to an existing point via R-tree
    DBG("Trying to link to an existing point via R-tree...");
    const auto nearestPoint = pointIndex_.qbegin(bgi::nearest(center, 1));
    if (nearestPoint != pointIndex_.qend()) {
        const auto& [point, pointId] = *nearestPoint;
        const double distance = bg::distance(point, center);
        MaybeTourNodeHandle nodeHandle = nodeForCircle_[pointId];
        if (nodeHandle && distance <= rad) {
            TourNode& node = requireNode(*nodeHandle);
            node.circles.insert(cid);
            node.energy += 3;
            node.insertions++;
            nodeForCircle_[cid] = *nodeHandle;

            DBG("Linked circle id=" << cid << " to existing node id=" << pointId);

            // Process the neighbors
            processNeighbor(node.prev);
            nodeHandle = nodeForCircle_[cid];
            if (!nodeHandle) {
                throw std::logic_error(
                    "circle lost its tour node during neighbor processing");
            }
            processNeighbor(requireNode(*nodeHandle).next);

            nodeHandle = nodeForCircle_[cid];
            if (!nodeHandle) {
                throw std::logic_error(
                    "circle lost its tour node during neighbor processing");
            }
            const TourNode& updatedNode = requireNode(*nodeHandle);
            if (updatedNode.insertions >= 2 &&
                (updatedNode.insertions & (updatedNode.insertions - 1)) == 0) {
                optimizePoint(*nodeHandle);
            }

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
        TreeNodeId uid = candidate->second;
        const MaybeTourNodeHandle uHandle = nodeForCircle_[uid];
        if (!uHandle) {
            continue;
        }
        const TourNode& u = requireNode(*uHandle);
        const TourNode& v = requireNode(u.next);
        Point cand = findOptimalPoint(center, rad, u.pos, v.pos);
        double addCost = bg::distance(u.pos, cand) +
            bg::distance(cand, v.pos) - bg::distance(u.pos, v.pos);
        if (addCost < bestAdd) {
            bestAdd = addCost;
            bestU = *uHandle;
            bestPt = cand;
        }
    }

    if (!bestU) {
        throw std::logic_error("failed to find a tour edge for circle insertion");
    }

    const TourNodeHandle leftHandle = *bestU;
    const TourNodeHandle rightHandle = requireNode(leftHandle).next;
    DBG("Best edge for insertion: between id=" << requireNode(leftHandle).id << " and id=" << requireNode(rightHandle).id
        << " at point (" << bg::get<0>(bestPt) << ", " << bg::get<1>(bestPt) << "), addCost=" << bestAdd);

    // 3) Insert new TourNode between bestU and bestU->next
    const TourNodeHandle nodeHandle = createNode(
        bestPt, cid, leftHandle, rightHandle);
    TourNode& left = requireNode(leftHandle);
    TourNode& right = requireNode(rightHandle);
    TourNode& node = requireNode(nodeHandle);

    // Remove old segment and insert new segments in R-tree
    segmentIndex_.remove({Segment(left.pos, right.pos), left.id});
    pointIndex_.insert({bestPt, cid});
    segmentIndex_.insert({Segment(left.pos, node.pos), left.id});
    segmentIndex_.insert({Segment(node.pos, right.pos), node.id});

    // Link in the new node
    left.next = nodeHandle;
    right.prev = nodeHandle;

    nodeForCircle_[cid] = nodeHandle;

    DBG("Inserted new TourNode id=" << cid << " between id=" << left.id << " and id=" << right.id);

    // Process the neighbors
    processNeighbor(node.prev);
    MaybeTourNodeHandle currentHandle = nodeForCircle_[cid];
    if (!currentHandle) {
        throw std::logic_error(
            "circle lost its tour node during neighbor processing");
    }
    processNeighbor(requireNode(*currentHandle).next);

    currentHandle = nodeForCircle_[cid];
    if (!currentHandle) {
        throw std::logic_error(
            "circle lost its tour node during neighbor processing");
    }
    const TourNode& updatedNode = requireNode(*currentHandle);
    if (updatedNode.insertions >= 2 &&
        (updatedNode.insertions & (updatedNode.insertions - 1)) == 0) {
        optimizePoint(*currentHandle);
    }
}

std::vector<Point> Tour::reconstruct() {
    const std::size_t N = treeNodes_.size();

    if (treeNodes_.empty()) {
        return {};
    }

    DBG("Starting unmerge with " << N << " tree nodes.");

    // Priority queue to process nodes in descending order of weight (merge cost)
    std::priority_queue<std::pair<double, TreeNodeId>> pq;
    pq.push({treeNodes_[N - 1].weight, N - 1}); // Start from the root

    // Insert the root node's circle into the tour
    DBG("Inserting root node id=" << (N - 1));
    insertCircle(N - 1);
    assertValid();

    // Process nodes in order of decreasing weight
    size_t nodesProcessed = 0;
    while (!pq.empty()) {
        auto [w, idx] = pq.top();
        pq.pop();
        nodesProcessed++;
        if (treeNodes_[idx].isLeaf()) continue;

        DBG("Unmerging node id=" << idx << " (weight=" << w
            << "), left=" << treeNodes_[idx].left
        << ", right=" << treeNodes_[idx].right);

        // Remove the current node's circle from the tour
        const MaybeTourNodeHandle tourNodeHandle = nodeForCircle_[idx];
        if (tourNodeHandle) {
            TourNode& tourNode = requireNode(*tourNodeHandle);
            auto& vc = tourNode.circles;
            // Remove this node's id from the circle list
            vc.erase(idx);
            nodeForCircle_[idx].reset();

            if (vc.empty()) {
                DBG("Deleting TourNode id=" << idx << " during unmerge.");
                deleteNode(*tourNodeHandle);
            } else if (idx == tourNode.id) {
                // We removed the main id; need to update id and R-trees
                TreeNodeId oldId = tourNode.id;
                TreeNodeId newId = *vc.begin();

                // Update id
                tourNode.id = newId;

                // Update rtreeP: remove old point, insert new point with new id
                pointIndex_.remove({tourNode.pos, oldId});
                pointIndex_.insert({tourNode.pos, newId});

                // Update rtreeS: remove old segments keyed by oldId, reinsert keyed by newId
                const TourNode& following = requireNode(tourNode.next);
                segmentIndex_.remove(
                    {Segment(tourNode.pos, following.pos), oldId});
                segmentIndex_.insert(
                    {Segment(tourNode.pos, following.pos), newId});

                // Also update previous node's segment if needed (segment keyed by previous node id)
                if (tourNode.prev != *tourNodeHandle) {
                    const TourNode& previous = requireNode(tourNode.prev);
                    segmentIndex_.remove(
                        {Segment(previous.pos, tourNode.pos), previous.id});
                    segmentIndex_.insert(
                        {Segment(previous.pos, tourNode.pos), previous.id});
                }

                DBG("Updated TourNode id from " << oldId << " to " << newId << " and updated R-trees.");
            }
        }
        // Insert the left and right children into the tour
        for (TreeNodeId c : treeNodes_[idx].children()) {
            DBG("Inserting child node id=" << c);
            insertCircle(c);
            // If the child is not a leaf, add it to the queue for further unmerging
            if (!treeNodes_[c].isLeaf()) {
                pq.push({treeNodes_[c].weight, c});
            }
        }

        // Maintenance is weighted by assigned-circle multiplicity. Scanning
        // nodeForCircle_ visits a tour node once per circle assigned to it, so
        // shared nodes accumulate optimization visits and consume energy faster.
        if ((nodesProcessed >= 2 && (nodesProcessed & (nodesProcessed - 1)) == 0) || pq.empty()) {
            // Optimize tour nodes with multiplicity determined by circle count.
            for (TreeNodeId nid = N; nid-- > 0;) {
                const MaybeTourNodeHandle handle = nodeForCircle_[nid];
                if (!handle) continue;
                TourNode& node = requireNode(*handle);
                node.insertions++;
                if (node.insertions >= 2 &&
                    (node.insertions & (node.insertions - 1)) == 0) {
                    DBG("Optimizing TourNode id=" << node.id << " at (" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");
                    optimizePoint(*handle);
                }
            }

            // Reduce energy once per associated circle ID and rebuild nodes
            // whose energy reaches zero.
            for (TreeNodeId nid = N; nid-- > 0;) {
                DBG("Processing node id=" << nid << " for energy reduction.");
                const MaybeTourNodeHandle handle = nodeForCircle_[nid];
                if (!handle) continue;
                TourNode& node = requireNode(*handle);
                if (node.energy > 0) {
                    --node.energy;
                }
                if (node.energy == 0) {
                    std::vector<TreeNodeId> linked(
                        node.circles.begin(), node.circles.end());
                    deleteNode(*handle);
                    sortByRadius(linked, treeNodes_);
                    for (TreeNodeId circleId : linked) {
                        insertCircle(circleId);
                    }
                }
            }

            // Repeat the same circle-weighted optimization after rebuilding.
            for (TreeNodeId nid = N; nid-- > 0;) {
                const MaybeTourNodeHandle handle = nodeForCircle_[nid];
                if (!handle) continue;
                TourNode& node = requireNode(*handle);
                node.insertions++;
                if (node.insertions >= 2 &&
                    (node.insertions & (node.insertions - 1)) == 0) {
                    DBG("Optimizing TourNode id=" << node.id << " at (" << bg::get<0>(node.pos) << ", " << bg::get<1>(node.pos) << ")");
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
