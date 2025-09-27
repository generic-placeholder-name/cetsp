#include "merge.hpp"
#include "debug.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <vector>

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
        std::vector<BoxValue> candidates;
        rtree_pr.query(bgi::covers(circles[id].center), std::back_inserter(candidates));
        for (auto& v : candidates) {
            Circle& cj = circles[v.second];
            if (bg::distance(cj.center, circles[id].center) + circles[id].r < cj.r) {
                removed[v.second] = true;
                rtree_pr.remove(std::make_pair(circleBox(cj), v.second));
            }
        }
        rtree_pr.insert(std::make_pair(circleBox(circles[id]), id));
    }
    
    // Filter out removed circles
    std::vector<Circle> tmp;
    tmp.reserve(n);
    for (auto& c : circles) {
        if (!removed[c.nodeId]) tmp.push_back(c);
    }
    circles = std::move(tmp);

    // Reset the node IDs
    for (size_t i = 0; i < circles.size(); ++i) {
        circles[i].nodeId = i;
    }

    DBG("Done removing");
    DBG("Remaining circles: " << circles.size());
    for (const auto& c : circles) {
        DBG("Circle ID: " << c.nodeId 
            << ", Center: (" << bg::get<0>(c.center) << ", " << bg::get<1>(c.center) << ")"
            << ", Radius: " << c.r);
    }
}

std::vector<TreeNode> buildMergeTree(std::vector<Circle> circles) {
    size_t n = circles.size();

    // Generate a random angle theta for rotation
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 2 * std::numbers::pi);
    double theta = dis(gen);

    // Rotate all points by theta
    auto rotatePoint = [&](Point& p, double angle) {
        double x = bg::get<0>(p);
        double y = bg::get<1>(p);
        bg::set<0>(p, x * std::cos(angle) - y * std::sin(angle));
        bg::set<1>(p, x * std::sin(angle) + y * std::cos(angle));
    };

    for (auto& c : circles) {
        rotatePoint(c.center, theta);
    }

    // Initialize combination-tree nodes for leaves
    std::vector<TreeNode> treeNodes;
    treeNodes.reserve(2 * n);
    for (size_t i = 0; i < circles.size(); ++i) {
        treeNodes.emplace_back(static_cast<size_t>(-1), static_cast<size_t>(-1), 0.0, circles[i].center, circles[i].r);
    }
    size_t nextNodeId = n;

    // Build R*-tree for merge-phase nearest-neighbor queries
    bgi::rtree<BoxValue, bgi::rstar<16>> rtree_kd;
    for (auto& c : circles) {
        rtree_kd.insert(std::make_pair(circleBox(c), c.nodeId));
    }

    // 2) Dynamic NN maintenance: each circle stores its nn and gap, and we keep a global min-heap
    std::set<std::pair<double, size_t>> heap;
    std::vector<std::pair<double, size_t>> handles(2 * n, {-1, -1});
    for (auto& c : circles) {
        // Find nearest neighbor via two-NN query
        constexpr int K = 16;
        std::vector<BoxValue> res;
        rtree_kd.query(bgi::nearest(c.center, K), std::back_inserter(res));
        size_t nn = static_cast<size_t>(-1);
        double best = std::numeric_limits<double>::infinity();
        for (auto& v : res) {
            size_t j = v.second;
            if (j == c.nodeId) continue;
            Circle& cj = circles[j];
            double d = gapDist(c.center, cj.center, c.r, cj.r);
            if (d < best) {
                best = d;
                nn = j;
            }
        }
        c.nn = nn;
        c.gap = best;
        auto it = heap.insert({best, c.nodeId}).first;
        handles[c.nodeId] = *it;
        circles[nn].rev.insert(c.nodeId);
    }

    DBG("Done initializing NNs");
    DBG("Nearest neighbors of each node:");
    for (const auto& c : circles) {
        DBG("Node " << c.nodeId << " -> NN: " << c.nn 
                  << " (gap: " << c.gap << ")");
    }

    size_t mergeCount = 0;
    // Combine loop
    for (size_t itr = 0; itr < n - 1; ++itr) {
        // Extract global minimum gap circle
        size_t id1;
        double bestDist;
        do {
            auto it = heap.begin();
            bestDist = it->first;
            id1 = it->second;
            heap.erase(it);
        } while (handles[id1].second == static_cast<size_t>(-1));
        Circle& c1 = circles[id1];
        size_t id2 = c1.nn;
        Circle& c2 = circles[id2];

        // Capture reverse-neighbor lists
        auto revs = c1.rev;
        revs.insert(c2.rev.begin(), c2.rev.end());
        revs.erase(id1); // remove self-reference
        revs.erase(id2); // remove self-reference

        size_t newId = nextNodeId++;
        ++mergeCount;
        DBG("Merge #" << mergeCount
              << ": " << id1 << " + " << id2
              << " (gap=" << bestDist << ") -> newId=" << newId);

        // Build combined circle
        Circle newC;
        auto combined = makeCombinedCircle(c1.center, c1.r, c2.center, c2.r);
        newC.center = combined.first;
        newC.r = combined.second;
        newC.nodeId = newId;
        newC.nn = static_cast<size_t>(-1);
        newC.gap = 0.0;
        newC.rev.clear();

        treeNodes.emplace_back(id1, id2, bestDist, newC.center, newC.r);

        // Remove old circles from tree & data structures
        handles[id1] = {-1, static_cast<size_t>(-1)};
        handles[id2] = {-1, static_cast<size_t>(-1)};
        circles[c2.nn].rev.erase(id2); // remove reverse link
        rtree_kd.remove(std::make_pair(circleBox(c1), id1));
        rtree_kd.remove(std::make_pair(circleBox(c2), id2));
        DBG("Removed old circles");

        // Insert combined circle
        circles.push_back(newC);
        rtree_kd.insert(std::make_pair(circleBox(newC), newId));

        // Compute NN for new circle
        std::vector<BoxValue> resn;
        rtree_kd.query(bgi::nearest(newC.center, 2),
                        std::back_inserter(resn));
        size_t nn3 = static_cast<size_t>(-1);
        double b3 = std::numeric_limits<double>::infinity();
        for (auto& v : resn) {
            size_t j = v.second;
            if (j == newId) continue;
            Circle& cj = circles[j];
            double d3 = gapDist(newC.center,
                                cj.center,
                                newC.r,
                                cj.r);
            if (d3 < b3) {
                b3 = d3;
                nn3 = j;
            }
        }
        circles[newId].nn = nn3;
        circles[newId].gap = b3;
        auto it3 = heap.insert({b3, newId}).first;
        handles[newId] = *it3;
        if (nn3 != static_cast<size_t>(-1)) {
            circles[nn3].rev.insert(newId);
        }
        DBG("  New circle " << newId
            << " nn=" << nn3
            << " gap=" << b3);

        // Recompute NN for circles that pointed to id1 or id2
        auto recompute = [&](size_t cid) {
            if (cid == static_cast<size_t>(-1)) return;
            heap.erase(handles[cid]);
            handles[cid] = {-1, static_cast<size_t>(-1)};

            std::vector<BoxValue> resc;
            rtree_kd.query(bgi::nearest(circles[cid].center, 2),
                           std::back_inserter(resc));
            size_t nn2 = static_cast<size_t>(-1);
            double b2 = std::numeric_limits<double>::infinity();
            for (auto& v : resc) {
                size_t j = v.second;
                if (j == cid) continue;
                Circle& cj = circles[j];
                double d2 = gapDist(circles[cid].center,
                                    cj.center,
                                    circles[cid].r,
                                    cj.r);
                if (d2 < b2) {
                    b2 = d2;
                    nn2 = j;
                }
            }
            circles[cid].nn = nn2;
            circles[cid].gap = b2;
            auto it2 = heap.insert({b2, cid}).first;
            handles[cid] = *it2;
            if (nn2 != static_cast<size_t>(-1)) {
                circles[nn2].rev.insert(cid);
            }
            DBG("  Update NN for circle " << cid
                << " -> nn=" << nn2
                << " gap=" << b2);
        };
        for (size_t cid : revs) recompute(cid);
    }

    // Rotate all TreeNode points back by -theta
    for (auto& node : treeNodes) {
        rotatePoint(node.center, -theta);
    }

    DBG("Combine complete. Total tree nodes=" << treeNodes.size());

    return treeNodes;
}