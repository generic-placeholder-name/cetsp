#include "merge.hpp"

#include "circle_geometry.hpp"
#include "debug.hpp"
#include "id_set.hpp"

#include <boost/container/small_vector.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace {

using BoxValue = std::pair<Box, TreeNodeId>;
using HeapEntry = std::pair<double, TreeNodeId>;

template<std::size_t Capacity>
using BoxValues = boost::container::small_vector<BoxValue, Capacity>;

constexpr TreeNodeId noNeighbor = std::numeric_limits<TreeNodeId>::max();

std::string formatNeighbor(TreeNodeId id) {
    return id == noNeighbor ? "<none>" : std::to_string(id);
}

struct ActiveCluster {
    Circle circle;
    TreeNodeId nearestNeighbor = noNeighbor;
    IdSet<TreeNodeId> reverseNeighbors;

    explicit ActiveCluster(Circle value) : circle(std::move(value)) {}
};

} // namespace

void removeCoveringCircles(std::vector<Circle>& circles) {
    size_t n = circles.size();

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return circles[a].r > circles[b].r;
    });

    bgi::rtree<BoxValue, bgi::rstar<16>> rtree_pr;
    std::vector<bool> removed(n, false);
    for (size_t id : order) {
        if (removed[id]) continue;
        BoxValues<16> candidates;
        rtree_pr.query(bgi::covers(circles[id].center), std::back_inserter(candidates));
        for (auto& v : candidates) {
            Circle& cj = circles[v.second];
            if (bg::distance(cj.center, circles[id].center) + circles[id].r <= cj.r) {
                removed[v.second] = true;
                rtree_pr.remove(std::make_pair(circleBox(cj), v.second));
            }
        }
        rtree_pr.insert(std::make_pair(circleBox(circles[id]), id));
    }
    
    // Filter out removed circles
    std::vector<Circle> tmp;
    tmp.reserve(n);
    for (size_t i = 0; i < circles.size(); ++i) {
        if (!removed[i]) tmp.push_back(std::move(circles[i]));
    }
    circles = std::move(tmp);

    DBG("Done removing");
    DBG("Remaining circles: " << circles.size());
    for (size_t i = 0; i < circles.size(); ++i) {
        const Circle& c = circles[i];
        DBG("Circle ID: " << i
            << ", Center: (" << bg::get<0>(c.center) << ", " << bg::get<1>(c.center) << ")"
            << ", Radius: " << c.r);
    }
}

std::vector<TreeNode> buildMergeTree(
    const std::vector<Circle>& circles,
    std::mt19937_64& randomEngine) {
    size_t n = circles.size();

    if (circles.empty()) {
        return {};
    }

    if (circles.size() == 1) {
        return {TreeNode::leaf(circles.front().center, circles.front().r)};
    }

    std::vector<ActiveCluster> activeClusters;
    activeClusters.reserve(2 * n - 1);
    for (const Circle& circle : circles) {
        activeClusters.emplace_back(circle);
    }

    // Generate a random angle theta for rotation
    std::uniform_real_distribution<> dis(0.0, 2 * std::numbers::pi);
    double theta = dis(randomEngine);

    const double cosine = std::cos(theta);
    const double sine = std::sin(theta);

    auto rotatePoint = [](Point& p, double cosAngle, double sinAngle) {
        double x = bg::get<0>(p);
        double y = bg::get<1>(p);
        bg::set<0>(p, x * cosAngle - y * sinAngle);
        bg::set<1>(p, x * sinAngle + y * cosAngle);
    };

    for (auto& active : activeClusters) {
        rotatePoint(active.circle.center, cosine, sine);
    }

    // Initialize combination-tree nodes for leaves
    std::vector<TreeNode> treeNodes;
    treeNodes.reserve(2 * n);
    for (size_t i = 0; i < n; ++i) {
        const Circle& circle = activeClusters[i].circle;
        treeNodes.push_back(TreeNode::leaf(circle.center, circle.r));
    }

    // Build R*-tree for merge-phase nearest-neighbor queries
    bgi::rtree<BoxValue, bgi::rstar<16>> rtree_kd;
    for (size_t id = 0; id < n; ++id) {
        rtree_kd.insert({circleBox(activeClusters[id].circle), id});
    }

    // 2) Dynamic NN maintenance: each active cluster stores its nearest
    // neighbor and reverse links, while the heap stores the corresponding gap.
    std::set<HeapEntry> heap;
    // An engaged entry records the current heap key for that circle. A
    // disengaged entry means any matching key still in the heap is stale.
    std::vector<std::optional<double>> currentHeapGap(2 * n);
    for (size_t id = 0; id < n; ++id) {
        ActiveCluster& active = activeClusters[id];
        const Circle& circle = active.circle;

        // Find nearest neighbor via a bounded R-tree query.
        constexpr int K = 16;
        BoxValues<16> res;
        rtree_kd.query(bgi::nearest(circle.center, K), std::back_inserter(res));
        TreeNodeId nearestNeighbor = noNeighbor;
        double best = std::numeric_limits<double>::infinity();
        for (auto& v : res) {
            size_t j = v.second;
            if (j == id) continue;
            const Circle& candidate = activeClusters[j].circle;
            double d = gapDist(
                circle.center, candidate.center, circle.r, candidate.r);
            if (d < best) {
                best = d;
                nearestNeighbor = j;
            }
        }
        active.nearestNeighbor = nearestNeighbor;
        if (nearestNeighbor == noNeighbor) {
            throw std::logic_error("failed to find a nearest neighbor while building merge tree");
        }
        heap.insert({best, id});
        currentHeapGap[id] = best;
        activeClusters[nearestNeighbor].reverseNeighbors.insert(id);
    }

    DBG("Done initializing NNs");
    DBG("Nearest neighbors of each node:");
    for (size_t id = 0; id < n; ++id) {
        DBG("Node " << id
            << " -> NN: " << activeClusters[id].nearestNeighbor
            << " (gap: " << *currentHeapGap[id] << ")");
    }

    size_t mergeCount = 0;
    // Combine loop
    for (size_t itr = 0; itr < n - 1; ++itr) {
        // Extract global minimum gap circle
        TreeNodeId id1;
        double bestDist;
        do {
            if (heap.empty()) {
                throw std::logic_error(
                    "merge queue became empty before the merge tree was complete");
            }
            auto it = heap.begin();
            bestDist = it->first;
            id1 = it->second;
            heap.erase(it);
        } while (!currentHeapGap[id1] || *currentHeapGap[id1] != bestDist);
        currentHeapGap[id1].reset();

        ActiveCluster& first = activeClusters[id1];
        if (first.nearestNeighbor == noNeighbor) {
            throw std::logic_error(
                "active merge cluster has no nearest neighbor");
        }
        TreeNodeId id2 = first.nearestNeighbor;
        if (id2 >= activeClusters.size() || !currentHeapGap[id2]) {
            throw std::logic_error(
                "selected nearest neighbor is not active");
        }
        ActiveCluster& second = activeClusters[id2];

        // Capture reverse-neighbor lists
        auto reverseNeighbors = first.reverseNeighbors;
        reverseNeighbors.insert(
            second.reverseNeighbors.begin(), second.reverseNeighbors.end());
        reverseNeighbors.erase(id1); // remove self-reference
        reverseNeighbors.erase(id2); // remove self-reference

        TreeNodeId newId = activeClusters.size();
        ++mergeCount;
        DBG("Merge #" << mergeCount
              << ": " << id1 << " + " << id2
              << " (gap=" << bestDist << ") -> newId=" << newId);

        // Build combined circle
        const Circle& firstCircle = first.circle;
        const Circle& secondCircle = second.circle;
        auto combined =
            makeCombinedCircle(
                firstCircle.center, firstCircle.r,
                secondCircle.center, secondCircle.r,
                randomEngine);
        Circle newCircle{combined.first, combined.second};

        treeNodes.push_back(TreeNode::branch(
            id1, id2, bestDist, newCircle.center, newCircle.r));

        // Remove old circles from tree & data structures
        currentHeapGap[id2].reset();
        if (second.nearestNeighbor == noNeighbor ||
            second.nearestNeighbor >= activeClusters.size()) {
            throw std::logic_error(
                "selected merge cluster has no active nearest neighbor");
        }
        activeClusters[second.nearestNeighbor].reverseNeighbors.erase(id2);
        rtree_kd.remove({circleBox(firstCircle), id1});
        rtree_kd.remove({circleBox(secondCircle), id2});
        DBG("Removed old circles");

        // Insert combined circle
        activeClusters.emplace_back(newCircle);
        rtree_kd.insert({circleBox(newCircle), newId});

        // Compute NN for new circle
        BoxValues<2> resn;
        rtree_kd.query(bgi::nearest(newCircle.center, 2),
                        std::back_inserter(resn));
        TreeNodeId nn3 = noNeighbor;
        double b3 = std::numeric_limits<double>::infinity();
        for (auto& v : resn) {
            size_t j = v.second;
            if (j == newId) continue;
            const Circle& candidate = activeClusters[j].circle;
            double d3 = gapDist(
                newCircle.center, candidate.center,
                newCircle.r, candidate.r);
            if (d3 < b3) {
                b3 = d3;
                nn3 = j;
            }
        }
        activeClusters[newId].nearestNeighbor = nn3;
        if (nn3 == noNeighbor && itr != n - 2) {
            throw std::logic_error(
                "merged circle has no neighbor before the final merge");
        }
        if (nn3 != noNeighbor) {
            heap.insert({b3, newId});
            currentHeapGap[newId] = b3;
            activeClusters[nn3].reverseNeighbors.insert(newId);
        }
        DBG("  New circle " << newId
            << " nn=" << formatNeighbor(nn3)
            << " gap=" << b3);

        // Recompute NN for circles that pointed to id1 or id2
        auto recompute = [&](TreeNodeId cid) {
            if (cid >= activeClusters.size() || !currentHeapGap[cid]) {
                throw std::logic_error(
                    "reverse-neighbor link does not name an active circle");
            }
            const std::size_t erased =
                heap.erase({*currentHeapGap[cid], cid});
            if (erased != 1) {
                throw std::logic_error(
                    "active circle's heap entry could not be removed");
            }
            currentHeapGap[cid].reset();

            ActiveCluster& active = activeClusters[cid];
            const Circle& circle = active.circle;
            BoxValues<2> resc;
            rtree_kd.query(bgi::nearest(circle.center, 2),
                           std::back_inserter(resc));
            TreeNodeId nn2 = noNeighbor;
            double b2 = std::numeric_limits<double>::infinity();
            for (auto& v : resc) {
                size_t j = v.second;
                if (j == cid) continue;
                const Circle& candidate = activeClusters[j].circle;
                double d2 = gapDist(
                    circle.center, candidate.center,
                    circle.r, candidate.r);
                if (d2 < b2) {
                    b2 = d2;
                    nn2 = j;
                }
            }
            if (nn2 == noNeighbor) {
                throw std::logic_error(
                    "failed to recompute an active circle's nearest neighbor");
            }
            active.nearestNeighbor = nn2;
            heap.insert({b2, cid});
            currentHeapGap[cid] = b2;
            activeClusters[nn2].reverseNeighbors.insert(cid);
            DBG("  Update NN for cluster " << cid
                << " -> nn=" << formatNeighbor(nn2)
                << " gap=" << b2);
        };
        for (TreeNodeId cid : reverseNeighbors) recompute(cid);
    }

    // Rotate all TreeNode points back by -theta
    for (auto& node : treeNodes) {
        rotatePoint(node.center, cosine, -sine);
    }

    DBG("Combine complete. Total tree nodes=" << treeNodes.size());

    return treeNodes;
}
