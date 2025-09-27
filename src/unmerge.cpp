#include "unmerge.hpp"
#include "debug.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <vector>

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

// Delete a TourNode from the cyclic list & any other data structures.
void deleteTourNode(
    TourNode*& head,
    TourNode* node,
    bgi::rtree<PointValue, bgi::rstar<16>>& rtreeP,
    bgi::rtree<SegmentValue, bgi::rstar<16>>& rtreeS,
    std::vector<TourNode*>& map
) {
    // Sanity check
    if (!node) return;

    TourNode* L = node->prev;
    TourNode* R = node->next;

    assert(L != nullptr && "Previous node should not be null");
    assert(R != nullptr && "Next node should not be null");

    DBG("Deleting TourNode id=" << node->id << " at (" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << ")");

    // Remove point and segment belonging to 'node'
    rtreeP.remove({node->pos, node->id});
    rtreeS.remove({Segment(node->pos, R->pos), node->id});

    if (R != node) { // If the list is a singleton, we can simply delete and leave
        // Splice out 'node': link L->R
        L->next = R;
        R->prev = L;

        // Update L's segment in R-tree
        rtreeS.remove({Segment(L->pos, node->pos), L->id});
        rtreeS.insert({Segment(L->pos, R->pos), L->id});

        DBG("Updated segment for L id=" << L->id << " to (" << bg::get<0>(L->pos) << ", " << bg::get<1>(L->pos) << ") -> (" << bg::get<0>(R->pos) << ", " << bg::get<1>(R->pos) << ")");
    }

    // Adjust head if needed
    if (head == node) {
        head = (R != node) ? R : nullptr;
        DBG("Head updated to id=" << (head ? head->id : -1));
    }

    // Clear the circles list and set map[cid] to nullptr
    for (size_t cid : node->circles) {
        map[cid] = nullptr;
    }

    DBG("About to delete node id=" << node->id);
    DBG("pos=(" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << ")");
    DBG("prev=" << (node->prev ? node->prev->id : -1)
        << " next=" << (node->next ? node->next->id : -1));
    DBG("circles.size() = " << node->circles.size());

    delete node;
    DBG("TourNode deleted.");
}

// Optimize the location of a point in the tour, possibly moving its circles to other points
void optimizeTourPoint(
    TourNode*& node,
    TourNode*& head,
    bgi::rtree<PointValue, bgi::rstar<16>>& rtreeP,
    bgi::rtree<SegmentValue, bgi::rstar<16>>& rtreeS,
    std::vector<TourNode*>& map,
    const std::vector<TreeNode>& treeNodes
) {
    if (!node || node == node->prev) return;
    constexpr double EPS = 1e-12;

    DBG("[optimizeTourPoint] Optimizing node id=" << node->id << " at (" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << "), circles=" << node->circles.size());

    // Remove from R-trees
    rtreeP.remove({node->pos, node->id});
    rtreeS.remove({Segment(node->prev->pos, node->pos), node->prev->id});
    rtreeS.remove({Segment(node->pos, node->next->pos), node->id});
 
    DBG("[optimizeTourPoint] Node " << node->id << " circles: ");
    for (size_t cid : node->circles) {
        DBG_NOENDL(cid << " ");
    }
    DBG("");
    DBG("[optimizeTourPoint] Node " << node->id
        << " prev: id=" << (node->prev ? node->prev->id : -1)
        << " pos=(" << (node->prev ? std::to_string(bg::get<0>(node->prev->pos)) : "null") << ", "
        << (node->prev ? std::to_string(bg::get<1>(node->prev->pos)) : "null") << ")"
        << " next: id=" << (node->next ? node->next->id : -1)
        << " pos=(" << (node->next ? std::to_string(bg::get<0>(node->next->pos)) : "null") << ", "
        << (node->next ? std::to_string(bg::get<1>(node->next->pos)) : "null") << ")"
    );

    // Try to move circles to other points
    HashSet<size_t> stillHere;
    for (size_t cid : node->circles) {
        bool moved = false;
        // Use rtreeP to search for a point within the circle
        std::vector<PointValue> candidates;
        rtreeP.query(bgi::nearest(treeNodes[cid].center, 1), std::back_inserter(candidates));
        for (const auto& [pt, pid] : candidates) {
            if (bg::distance(pt, treeNodes[cid].center) <= treeNodes[cid].r) {
                TourNode* other = map[pid];
                assert(other != node && "Cannot move circle to itself");
                assert(other != nullptr && "Other node should not be null");
                other->circles.insert(cid);
                map[cid] = other;
                moved = true;
                DBG("[optimizeTourPoint] Moved circle " << cid << " from node " << node->id << " to node " << other->id);
                break;
            }
        }
        if (!moved) stillHere.insert(cid);
    }
    node->circles = std::move(stillHere);

    // If no circles left, delete node and return
    if (node->circles.empty()) {
        DBG("[optimizeTourPoint] No circles left in node " << node->id << ", deleting node.");
        deleteTourNode(head, node, rtreeP, rtreeS, map);
        return;
    }

    // Get previous and next points
    TourNode* aNode = node->prev;
    TourNode* bNode = node->next;
    Point a = aNode->pos, b = bNode->pos;
    
    // Check if [a, b] intersects the intersection of all circles
    double ax = bg::get<0>(a), ay = bg::get<1>(a);
    double bx = bg::get<0>(b), by = bg::get<1>(b);
    double dx = bx - ax, dy = by - ay;

    // [Tmin, Tmax] ∩= each circle‐interval
    double T_min = 0.0, T_max = 1.0;

    if (bg::distance(a, b) < EPS) {
        // If a and b are the same point, the quadratic equation degenerates
        // Therefore, we have to check it separately
        for (size_t cid : node->circles) {
            const auto& C = treeNodes[cid].center;
            double r = treeNodes[cid].r;

            // Check if the circle contains point a (or b, since a == b)
            if (bg::distance(a, C) > r)  {
                T_min = 1.0; T_max = 0.0; // No intersection
                break;
            }
        }
    } else {
        for (size_t cid : node->circles) {
            const auto& C  = treeNodes[cid].center;
            double  cx     = bg::get<0>(C),
                    cy     = bg::get<1>(C),
                    r      = treeNodes[cid].r;

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
        node->pos = P;
    } else {
        // Projected gradient descent step
        // Gradient: sum of (negative) gradients of distances to a and b
        double gx = 0, gy = 0;
        double x = bg::get<0>(node->pos), y = bg::get<1>(node->pos);
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
        for (size_t cid : node->circles) {
            const TreeNode& c = treeNodes[cid];
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
            DBG("[optimizeTourPoint] No valid step found for node " << node->id << ", skipping optimization.");
        }
        else {
            DBG("[optimizeTourPoint] Gradient step for node " << node->id << ": gx=" << gx << " gy=" << gy << " step=" << maxStep);
            node->pos = Point(x + gx * maxStep, y + gy * maxStep);
        }
    }

    // If id is not in the circle anymore, reset it
    if (!node->circles.contains(node->id)) {
        size_t oldId = node->id;
        node->id = *node->circles.begin();
        DBG("[optimizeTourPoint] Changed node id from " << oldId << " to " << node->id);
    }

    // Reinsert into R-trees
    rtreeP.insert({node->pos, node->id});
    rtreeS.insert({Segment(node->pos, node->next->pos), node->id});
    rtreeS.insert({Segment(node->prev->pos, node->pos), node->prev->id});

    DBG("[optimizeTourPoint] Optimization complete for node id=" << node->id << " at (" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << ")");
}

// Insert the circle `cid` (with data in `tn`) into the current tour starting at `head`.
// Updates energies, and if any neighbor hits energy=0, recursively removes it and
// reinserts its linked circles in ascending‐radius order.
// `map[cid]` will be set to wherever the circle ultimately lives.
// Insert the circle `cid` into the current tour, updating R-trees and energies.
void insertCircle(
    size_t cid,
    const TreeNode& tn,
    TourNode*& head,
    bgi::rtree<PointValue, bgi::rstar<16>>& rtreeP,
    bgi::rtree<SegmentValue, bgi::rstar<16>>& rtreeS,
    std::vector<TourNode*>& map,
    const std::vector<TreeNode>& allNodes
) {
    const Point& center = tn.center;
    double rad = tn.r;

    DBG(rtreeP.size() << " points in R-tree, " << rtreeS.size() << " segments in R-tree.");
    DBG("Inserting circle id=" << cid << " at (" << bg::get<0>(center) << ", " << bg::get<1>(center) << "), r=" << rad);

    // Empty tour: create first node
    if (!head) {
        TourNode* node = new TourNode(center, cid);
        node->prev = node;
        node->next = node;
        head = node;
        map[cid] = node;
        rtreeP.insert({center, cid});
        // For a single point, use a degenerate segment (point, point)
        rtreeS.insert({Segment(center, center), cid});
        DBG("Created first TourNode id=" << cid);
        return;
    }

    // Utility function: reduce the energy of a node and trigger a reinsert if it hits zero
    auto processNeighbor = [&](TourNode*& neighborPtr) {
        neighborPtr->energy--;
        if (neighborPtr->energy == 0) {
            DBG("Neighbor id=" << neighborPtr->id << " has zero energy, removing and reinserting its circles.");
            std::vector<size_t> linked(neighborPtr->circles.begin(), neighborPtr->circles.end());
            deleteTourNode(head, neighborPtr, rtreeP, rtreeS, map);
            std::sort(linked.begin(), linked.end(), [&](size_t a, size_t b) { return allNodes[a].r < allNodes[b].r; });
            for (size_t oc : linked) {
                insertCircle(oc, allNodes[oc], head, rtreeP, rtreeS, map, allNodes);
            }
        }
    };

    // 1) Try to link to an existing point via R-tree
    DBG("Trying to link to an existing point via R-tree...");
    std::vector<PointValue> qres;
    rtreeP.query(bgi::nearest(center, 1), std::back_inserter(qres));
    if (!qres.empty()) {
        auto [pt, pid] = qres[0];
        if (bg::distance(pt, center) <= rad) {
            TourNode* node = map[pid];
            assert(node != nullptr);
            node->circles.insert(cid);
            node->energy += 3; 
            node->insertions++;
            map[cid] = node;

            DBG("Linked circle id=" << cid << " to existing node id=" << pid);

            // Process the neighbors
            processNeighbor(node->prev);
            node = map[cid]; // The delete-reinsert cascade might have invalidated this pointer
            processNeighbor(node->next);
            
            node = map[cid]; // Re-fetch the node after potential deletion
            if (node->insertions >= 2 && (node->insertions & (node->insertions - 1)) == 0) {
                optimizeTourPoint(node, head, rtreeP, rtreeS, map, allNodes);
            }

            return;
        }
    }

    // 2) Find best edge to insert via segment R-tree
    DBG("Finding best edge to insert via segment R-tree...");
    constexpr int K = 16;
    std::vector<SegmentValue> sres;
    rtreeS.query(bgi::nearest(center, K), std::back_inserter(sres));
    double bestAdd = std::numeric_limits<double>::infinity();
    TourNode* bestU = nullptr;
    Point bestPt;

    for (auto& sv : sres) {
        size_t uid = sv.second;
        TourNode* u = map[uid];
        TourNode* v = u->next;
        Point cand = findOptimalPoint(center, rad, u->pos, v->pos);
        double addCost = bg::distance(u->pos, cand) + bg::distance(cand, v->pos) - bg::distance(u->pos, v->pos);
        if (addCost < bestAdd) {
            bestAdd = addCost;
            bestU = u;
            bestPt = cand;
        }
    }

    DBG("Best edge for insertion: between id=" << bestU->id << " and id=" << bestU->next->id
        << " at point (" << bg::get<0>(bestPt) << ", " << bg::get<1>(bestPt) << "), addCost=" << bestAdd);

    // 3) Insert new TourNode between bestU and bestU->next
    TourNode* L = bestU;
    TourNode* R = bestU->next;
    TourNode* node = new TourNode(bestPt, cid);

    // Remove old segment and insert new segments in R-tree
    rtreeS.remove({Segment(L->pos, R->pos), L->id});
    rtreeP.insert({bestPt, cid});
    rtreeS.insert({Segment(L->pos, node->pos), L->id});
    rtreeS.insert({Segment(node->pos, R->pos), node->id});

    // Link in the new node
    L->next = node;
    node->prev = L;
    node->next = R;
    R->prev = node;

    map[cid] = node;
    node->energy += 3;
    node->insertions++;

    DBG("Inserted new TourNode id=" << cid << " between id=" << L->id << " and id=" << R->id);

    // Process the neighbors
    processNeighbor(node->prev);
    node = map[cid]; // The delete-reinsert cascade might have invalidated this pointer
    processNeighbor(node->next);
    
    node = map[cid]; // Re-fetch the node after potential deletion
    if (node->insertions >= 2 && (node->insertions & (node->insertions - 1)) == 0) {
        optimizeTourPoint(node, head, rtreeP, rtreeS, map, allNodes);
    }
}

std::vector<Point> unmerge(std::vector<TreeNode> treeNodes) {
    size_t N = treeNodes.size();

    DBG("Starting unmerge with " << N << " tree nodes.");

    // Priority queue to process nodes in descending order of weight (merge cost)
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({treeNodes[N - 1].weight, N - 1}); // Start from the root

    TourNode* head = nullptr; // Head of the cyclic tour
    bgi::rtree<PointValue, bgi::rstar<16>> rtreeP; // R-tree for tour points
    bgi::rtree<SegmentValue, bgi::rstar<16>> rtreeS; // R-tree for tour segments
    std::vector<TourNode*> map(N, nullptr); // Map from node id to TourNode*

    // Insert the root node's circle into the tour
    DBG("Inserting root node id=" << (N - 1));
    insertCircle(N - 1, treeNodes[N - 1], head, rtreeP, rtreeS, map, treeNodes);

    // Process nodes in order of decreasing weight
    size_t nodesProcessed = 0;
    while (!pq.empty()) {
        auto [w, idx] = pq.top();
        pq.pop();
        nodesProcessed++;
        if (treeNodes[idx].left == static_cast<size_t>(-1)) continue; // Skip leaves

        DBG("Unmerging node id=" << idx << " (weight=" << w << "), left=" << treeNodes[idx].left << ", right=" << treeNodes[idx].right);

        // Remove the current node's circle from the tour
        TourNode* tnode = map[idx];
        if (tnode) {
            auto& vc = tnode->circles;
            // Remove this node's id from the circle list
            vc.erase(idx);
            map[idx] = nullptr;

            if (vc.empty()) {
                DBG("Deleting TourNode id=" << idx << " during unmerge.");
                deleteTourNode(head, tnode, rtreeP, rtreeS, map);
            } else if (idx == tnode->id) {
                // We removed the main id; need to update id and R-trees
                size_t oldId = tnode->id;
                size_t newId = *vc.begin();

                // Update id
                tnode->id = newId;

                // Update rtreeP: remove old point, insert new point with new id
                rtreeP.remove({tnode->pos, oldId});
                rtreeP.insert({tnode->pos, newId});

                // Update rtreeS: remove old segments keyed by oldId, reinsert keyed by newId
                rtreeS.remove({Segment(tnode->pos, tnode->next->pos), oldId});
                rtreeS.insert({Segment(tnode->pos, tnode->next->pos), newId});

                // Also update previous node's segment if needed (segment keyed by previous node id)
                if (tnode->prev && tnode->prev != tnode) {
                    rtreeS.remove({Segment(tnode->prev->pos, tnode->pos), tnode->prev->id});
                    rtreeS.insert({Segment(tnode->prev->pos, tnode->pos), tnode->prev->id});
                }

                DBG("Updated TourNode id from " << oldId << " to " << newId << " and updated R-trees.");
            }
        }

        // Insert the left and right children into the tour
        for (size_t c : {treeNodes[idx].left, treeNodes[idx].right}) {
            DBG("Inserting child node id=" << c);
            insertCircle(c, treeNodes[c], head, rtreeP, rtreeS, map, treeNodes);
            // If the child is not a leaf, add it to the queue for further unmerging
            if (treeNodes[c].left != static_cast<size_t>(-1)) {
                pq.push({treeNodes[c].weight, c});
            }
        }

        // Once in a while, optimize all nodes in the tour
        if ((nodesProcessed >= 2 && (nodesProcessed & (nodesProcessed - 1)) == 0) || pq.empty()) {
            // Optimize all points in the tour
            for (size_t nid = N; nid-- > 0;) {
                if (map[nid] == nullptr) continue; // Skip if not in the tour
                TourNode* t = map[nid];
                if (t) {
                    t->insertions++;
                    if (t->insertions >= 2 && (t->insertions & (t->insertions - 1)) == 0) {
                        DBG("Optimizing TourNode id=" << t->id << " at (" << bg::get<0>(t->pos) << ", " << bg::get<1>(t->pos) << ")");
                        optimizeTourPoint(t, head, rtreeP, rtreeS, map, treeNodes);
                    }
                }
            }
            
            // Interleave optimization with energy reduction
            // For every node in the pq, decrease energy and handle zero-energy nodes
            for (size_t nid = N; nid-- > 0;) {
                DBG("Processing node id=" << nid << " for energy reduction.");
                if (map[nid] == nullptr) continue; // Skip if not in the tour
                TourNode* t = map[nid];
                if (t) {
                    t->energy--;
                    if (t->energy == 0) {
                        std::vector<size_t> linked(t->circles.begin(), t->circles.end());
                        deleteTourNode(head, t, rtreeP, rtreeS, map);
                        std::sort(linked.begin(), linked.end(), [&](size_t a, size_t b) { return treeNodes[a].r < treeNodes[b].r; });
                        for (size_t oc : linked) {
                            insertCircle(oc, treeNodes[oc], head, rtreeP, rtreeS, map, treeNodes);
                        }
                    }
                }
            }

            // Optimize all points in the tour
            for (size_t nid = N; nid-- > 0;) {
                if (map[nid] == nullptr) continue; // Skip if not in the tour
                TourNode* t = map[nid];
                if (t) {
                    t->insertions++;
                    if (t->insertions >= 2 && (t->insertions & (t->insertions - 1)) == 0) {
                        DBG("Optimizing TourNode id=" << t->id << " at (" << bg::get<0>(t->pos) << ", " << bg::get<1>(t->pos) << ")");
                        optimizeTourPoint(t, head, rtreeP, rtreeS, map, treeNodes);
                    }
                }
            }
        }
    }

    // Collect the tour points in order, starting from head
    std::vector<Point> tour;
    if (head) {
        TourNode* cur = head;
        do {
            tour.push_back(cur->pos);
            cur = cur->next;
        } while (cur != head);
    }

    DBG("Tour reconstruction complete. Tour size: " << tour.size());

    // Clean up: delete all TourNodes to prevent memory leaks
    if (head) {
        TourNode* cur = head->next;
        while (cur != head) {
            TourNode* nxt = cur->next;
            delete cur;
            cur = nxt;
        }
        delete head;
    }

    return tour;
}