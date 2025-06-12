// This is the algorithm
// Code works and is compilable
// Pending refactor

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
typedef bg::model::segment<Point> Segment;
typedef bg::model::box<Point> Box;

typedef std::pair<Point, int> PointValue;
typedef std::pair<Box, int> BoxValue;
typedef std::pair<Segment, int> SegmentValue;

struct TreeNode { 
    int left, right; 
    double weight; 
    Point center;
    double r;

    TreeNode(int l, int r, double w, Point c = Point(), double ra = 0.0) 
        : left(l), right(r), weight(w), center(c), r(ra) {}
};

struct Circle {
    Point center;
    double r;
    int nodeId;
    int nn;
    double gap;
    std::unordered_set<int> rev;
};

// Create bounding box around a circle
Box circleBox(const Circle& c) {
    double adjr = c.r + 1e-12; // Add a small epsilon to avoid precision issues
    return Box(
        Point(bg::get<0>(c.center) - adjr, bg::get<1>(c.center) - adjr),
        Point(bg::get<0>(c.center) + adjr, bg::get<1>(c.center) + adjr)
    );
}

// Compute "gap" distance (center-dist minus sum of radii)
double gapDist(const Point& a, const Point& b, double ra, double rb) {
    return bg::distance(a, b) - (ra + rb);
}

// Compute a combined "representative" circle from two circles. 
std::pair<Point, double> makeCombinedCircle(const Point& p1, double r1, const Point& p2, double r2) {
    double x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
    double x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
    double dx = x2 - x1, dy = y2 - y1;
    double d = std::hypot(dx, dy);

    // If one circle is completely inside the other
    if (d + std::min(r1, r2) <= std::max(r1, r2)) {
        return (r1 < r2) ? std::make_pair(p1, r1) : std::make_pair(p2, r2);
    }

    // Normalize direction vector
    double ux = dx / d;
    double uy = dy / d;

    // Points on the edges of the circles along the line connecting centers
    double xEdge1 = x1 + ux * r1;
    double yEdge1 = y1 + uy * r1;
    double xEdge2 = x2 - ux * r2;
    double yEdge2 = y2 - uy * r2;

    // Midpoint of these edge points is the new center
    Point center((xEdge1 + xEdge2) / 2.0, (yEdge1 + yEdge2) / 2.0);

    // Handle non-overlapping case (set radius = 0)
    if (d >= r1 + r2) {
        return std::make_pair(center, 0.0);
    }

    // Circles intersect — compute potential radius range
    double overlapDepth = (r1 + r2 - d) / 2.0;

    // Compute distance from center of first circle to chord midpoint
    double a = (r1*r1 - r2*r2 + d*d) / (2*d);
    double h = std::sqrt(std::max(0.0, r1*r1 - a*a));  // half chord length

    // Radius range: [overlapDepth, h]
    double minRadius = std::min(overlapDepth, h);
    double maxRadius = std::max(overlapDepth, h);
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minRadius, maxRadius);
    double radius = dis(gen);

    return std::make_pair(center, radius);
}

void removeCoveringCircles(std::vector<Circle>& circles) {
    // Pre-remove in O(n log n) via R*-tree
    int n = circles.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return circles[a].r > circles[b].r;
    });

    bgi::rtree<BoxValue, bgi::rstar<16>> rtree_pr;
    std::vector<bool> removed(n, false);
    for (int id : order) {
        if (removed[id]) continue;
        std::vector<BoxValue> candidates;
        rtree_pr.query(bgi::covers(circles[id].center), std::back_inserter(candidates));
        for (auto& v : candidates) {
            Circle& cj = circles[v.second];
            if (bg::distance(cj.center, circles[id].center) + circles[id].r < cj.r) {
                // Remove the larger circle
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

    std::cerr << "Done removing" << std::endl;
    std::cerr << "Remaining circles: " << circles.size() << std::endl;
    for (const auto& c : circles) {
        std::cerr << "Circle ID: " << c.nodeId 
                  << ", Center: (" << bg::get<0>(c.center) << ", " << bg::get<1>(c.center) << ")"
                  << ", Radius: " << c.r << std::endl;
    }
}

// Build a merge tree from the circles.
// This function constructs a tree where each node represents a merge operation
// between two circles, storing the gap distance and the IDs of the merged circles.
std::vector<TreeNode> buildMergeTree(std::vector<Circle> circles) {
    int n = circles.size();

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
        treeNodes.emplace_back(-1, -1, 0.0, circles[i].center, circles[i].r);
    }
    int nextNodeId = n;

    // Build R*-tree for merge-phase nearest-neighbor queries
    bgi::rtree<BoxValue, bgi::rstar<16>> rtree_kd;
    for (auto& c : circles) {
        rtree_kd.insert(std::make_pair(circleBox(c), c.nodeId));
    }

    // 2) Dynamic NN maintenance: each circle stores its nn and gap, and we keep a global min-heap
    std::set<std::pair<double, int>> heap;
    std::vector<std::pair<double, int>> handles(2 * n, {-1, -1});
    for (auto& c : circles) {
        // find nearest neighbor via two-NN query
        constexpr int K = 16;
        std::vector<BoxValue> res;
        rtree_kd.query(bgi::nearest(c.center, K), std::back_inserter(res));
        int nn = -1;
        double best = std::numeric_limits<double>::infinity();
        for (auto& v : res) {
            int j = v.second;
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

    std::cerr << "Done initializing NNs" << std::endl;
    std::cerr << "Nearest neighbors of each node:" << std::endl;
    for (const auto& c : circles) {
        std::cerr << "Node " << c.nodeId << " -> NN: " << c.nn 
                  << " (gap: " << c.gap << ")" << std::endl;
    }

    int mergeCount = 0;
    // Combine loop
    for (int itr = 0; itr < n - 1; ++itr) {
        // Extract global minimum gap circle
        int id1;
        double bestDist;
        do {
            auto it = heap.begin();
            bestDist = it->first;
            id1 = it->second;
            heap.erase(it);
        } while (handles[id1].second == -1);
        Circle& c1 = circles[id1];
        int id2 = c1.nn;
        Circle& c2 = circles[id2];

        // Capture reverse-neighbor lists
        auto revs = c1.rev;
        revs.insert(c2.rev.begin(), c2.rev.end());
        revs.erase(id1); // remove self-reference
        revs.erase(id2); // remove self-reference

        int newId = nextNodeId++;
        ++mergeCount;
        std::cerr << "Merge #" << mergeCount
                  << ": " << id1 << " + " << id2
                  << " (gap=" << bestDist << ") -> newId=" << newId << std::endl;

        // Build combined circle
        Circle newC;
        auto combined = makeCombinedCircle(c1.center, c1.r, c2.center, c2.r);
        newC.center = combined.first;
        newC.r = combined.second;
        newC.nodeId = newId;
        newC.nn = -1;
        newC.gap = 0.0;
        newC.rev.clear();

        treeNodes.emplace_back(id1, id2, bestDist, newC.center, newC.r);

        // Remove old circles from tree & data structures
        handles[id1] = {-1, -1};
        handles[id2] = {-1, -1};
        circles[c2.nn].rev.erase(id2); // remove reverse link
        rtree_kd.remove(std::make_pair(circleBox(c1), id1));
        rtree_kd.remove(std::make_pair(circleBox(c2), id2));
        std::cerr << "Removed old circles" << std::endl;

        // Insert combined circle
        circles.push_back(newC);
        rtree_kd.insert(std::make_pair(circleBox(newC), newId));

        // Compute NN for new circle
        std::vector<BoxValue> resn;
        rtree_kd.query(bgi::nearest(newC.center, 2),
                        std::back_inserter(resn));
        int nn3 = -1;
        double b3 = std::numeric_limits<double>::infinity();
        for (auto& v : resn) {
            int j = v.second;
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
        if (nn3 != -1) {
            circles[nn3].rev.insert(newId);
        }
        std::cerr << "  New circle " << newId
                    << " nn=" << nn3
                    << " gap=" << b3 << std::endl;

        // Recompute NN for circles that pointed to id1 or id2
        auto recompute = [&](int cid) {
            if (cid == -1) return;
            heap.erase(handles[cid]);
            handles[cid] = {-1, -1};

            std::vector<BoxValue> resc;
            rtree_kd.query(bgi::nearest(circles[cid].center, 2),
                           std::back_inserter(resc));
            int nn2 = -1;
            double b2 = std::numeric_limits<double>::infinity();
            for (auto& v : resc) {
                int j = v.second;
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
            if (nn2 != -1) {
                circles[nn2].rev.insert(cid);
            }
            std::cerr << "  Update NN for circle " << cid
                      << " -> nn=" << nn2
                      << " gap=" << b2 << std::endl;
        };
        for (int cid : revs) recompute(cid);
    }

    // Rotate all TreeNode points back by -theta
    for (auto& node : treeNodes) {
        rotatePoint(node.center, -theta);
    }

    // Print the root for debugging
    int root = circles.back().nodeId;
    std::cerr << "Combine complete. Total tree nodes=" << treeNodes.size() << std::endl;

    return treeNodes;
}

// Node in the cyclic tour
// Contains a point, list of circles attached, energy counter, and pointers
struct TourNode {
    Point pos;
    std::vector<int> circles;
    int energy;
    int insertions;
    TourNode* prev;
    TourNode* next;
    int id; // circle index
    
    TourNode(const Point& p, int id_ = -1)
        : pos(p), energy(0), insertions(0), prev(nullptr), next(nullptr), id(id_), circles(1, id_) {}
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

    std::cerr << "Deleting TourNode id=" << node->id << " at (" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << ")" << std::endl;

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

        std::cerr << "Updated segment for L id=" << L->id << " to (" << bg::get<0>(L->pos) << ", " << bg::get<1>(L->pos) << ") -> (" << bg::get<0>(R->pos) << ", " << bg::get<1>(R->pos) << ")" << std::endl;
    }

    // Adjust head if needed
    if (head == node) {
        head = (R != node) ? R : nullptr;
        std::cerr << "Head updated to id=" << (head ? head->id : -1) << std::endl;
    }

    // Clear the circles list and set map[cid] to nullptr
    for (int cid : node->circles) {
        map[cid] = nullptr;
    }

    std::cerr << "About to delete node id=" << node->id << std::endl;
    std::cerr << "pos=(" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << ")" << std::endl;
    std::cerr << "prev=" << (node->prev ? node->prev->id : -1)
            << " next=" << (node->next ? node->next->id : -1) << std::endl;
    std::cerr << "circles.size() = " << node->circles.size() << std::endl;

    delete node;
    std::cerr << "TourNode deleted." << std::endl;
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

    std::cerr << "[optimizeTourPoint] Optimizing node id=" << node->id << " at (" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << "), circles=" << node->circles.size() << std::endl;

    // Remove from R-trees
    rtreeP.remove({node->pos, node->id});
    rtreeS.remove({Segment(node->prev->pos, node->pos), node->prev->id});
    rtreeS.remove({Segment(node->pos, node->next->pos), node->id});
 
    std::cerr << "[optimizeTourPoint] Node " << node->id << " circles: ";
    for (int cid : node->circles) {
        std::cerr << cid << " ";
    }
    std::cerr << std::endl;
    std::cerr << "[optimizeTourPoint] Node " << node->id << " prev: id=" << (node->prev ? node->prev->id : -1)
              << " pos=(" << (node->prev ? std::to_string(bg::get<0>(node->prev->pos)) : "null") << ", "
              << (node->prev ? std::to_string(bg::get<1>(node->prev->pos)) : "null") << ")"
              << " next: id=" << (node->next ? node->next->id : -1)
              << " pos=(" << (node->next ? std::to_string(bg::get<0>(node->next->pos)) : "null") << ", "
              << (node->next ? std::to_string(bg::get<1>(node->next->pos)) : "null") << ")"
              << std::endl;

    // Try to move circles to other points
    std::vector<int> stillHere;
    for (int cid : node->circles) {
        bool moved = false;
        // Use rtreeP to search for a point within the circle
        std::vector<PointValue> candidates;
        rtreeP.query(bgi::nearest(treeNodes[cid].center, 1), std::back_inserter(candidates));
        for (const auto& [pt, pid] : candidates) {
            if (bg::distance(pt, treeNodes[cid].center) <= treeNodes[cid].r) {
                TourNode* other = map[pid];
                assert(other != node && "Cannot move circle to itself");
                assert(other != nullptr && "Other node should not be null");
                other->circles.push_back(cid);
                map[cid] = other;
                moved = true;
                std::cerr << "[optimizeTourPoint] Moved circle " << cid << " from node " << node->id << " to node " << other->id << std::endl;
                break;
            }
        }
        if (!moved) stillHere.push_back(cid);
    }
    node->circles = stillHere;

    // If no circles left, delete node and return
    if (node->circles.empty()) {
        std::cerr << "[optimizeTourPoint] No circles left in node " << node->id << ", deleting node." << std::endl;
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
        for (int cid : node->circles) {
            const auto& C = treeNodes[cid].center;
            double r = treeNodes[cid].r;

            // Check if the circle contains point a (or b, since a == b)
            if (bg::distance(a, C) > r)  {
                T_min = 1.0; T_max = 0.0; // No intersection
                break;
            }
        }
    } else {
        for (int cid : node->circles) {
            const auto& C  = treeNodes[cid].center;
            double   cx     = bg::get<0>(C),
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
        std::cerr << "[optimizeTourPoint] Found common-intersection point at t=" 
                  << t_star << "\n";
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
        for (int cid : node->circles) {
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
            std::cerr << "[optimizeTourPoint] No valid step found for node " << node->id << ", skipping optimization." << std::endl;
        }
        else {
            std::cerr << "[optimizeTourPoint] Gradient step for node " << node->id << ": gx=" << gx << " gy=" << gy << " step=" << maxStep << std::endl;
            node->pos = Point(x + gx * maxStep, y + gy * maxStep);
        }
    }

    // If id is not the smallest circle id, reset it
    if (std::find(node->circles.begin(), node->circles.end(), node->id) == node->circles.end()) {
        int oldId = node->id;
        node->id = node->circles.front();
        std::cerr << "[optimizeTourPoint] Changed node id from " << oldId << " to " << node->id << std::endl;
    }

    // Reinsert into R-trees
    rtreeP.insert({node->pos, node->id});
    rtreeS.insert({Segment(node->pos, node->next->pos), node->id});
    rtreeS.insert({Segment(node->prev->pos, node->pos), node->prev->id});

    std::cerr << "[optimizeTourPoint] Optimization complete for node id=" << node->id << " at (" << bg::get<0>(node->pos) << ", " << bg::get<1>(node->pos) << ")" << std::endl;
}

int totInserts = 0; // Global counter for debugging
// Insert the circle `cid` (with data in `tn`) into the current tour starting at `head`.
// Updates energies, and if any neighbor hits energy=0, recursively removes it and
// reinserts its linked circles in ascending‐radius order.
// `map[cid]` will be set to wherever the circle ultimately lives.
// Insert the circle `cid` into the current tour, updating R-trees and energies.
void insertCircle(
    int cid,
    const TreeNode& tn,
    TourNode*& head,
    bgi::rtree<PointValue, bgi::rstar<16>>& rtreeP,
    bgi::rtree<SegmentValue, bgi::rstar<16>>& rtreeS,
    std::vector<TourNode*>& map,
    const std::vector<TreeNode>& allNodes
) {
    totInserts++;
    const Point& center = tn.center;
    double rad = tn.r;

    std::cerr << rtreeP.size() << " points in R-tree, " << rtreeS.size() << " segments in R-tree." << std::endl;

    std::cerr << "Inserting circle id=" << cid << " at (" << bg::get<0>(center) << ", " << bg::get<1>(center) << "), r=" << rad << std::endl;

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
        std::cerr << "Created first TourNode id=" << cid << std::endl;
        return;
    }

    // Utility function: reduce the energy of a node and trigger a reinsert if it hits zero
    auto processNeighbor = [&](TourNode*& neighborPtr) {
        neighborPtr->energy--;
        if (neighborPtr->energy == 0) {
            std::cerr << "Neighbor id=" << neighborPtr->id << " has zero energy, removing and reinserting its circles." << std::endl;
            auto linked = neighborPtr->circles;
            deleteTourNode(head, neighborPtr, rtreeP, rtreeS, map);
            std::sort(linked.begin(), linked.end(), [&](int a, int b) { return allNodes[a].r < allNodes[b].r; });
            for (int oc : linked) {
                insertCircle(oc, allNodes[oc], head, rtreeP, rtreeS, map, allNodes);
            }
        }
    };

    // 1) Try to link to an existing point via R-tree
    std::cerr << "Trying to link to an existing point via R-tree..." << std::endl;
    std::vector<PointValue> qres;
    rtreeP.query(bgi::nearest(center, 1), std::back_inserter(qres));
    if (!qres.empty()) {
        auto [pt, pid] = qres[0];
        if (bg::distance(pt, center) <= rad) {
            TourNode* node = map[pid];
            assert(node != nullptr);
            node->circles.push_back(cid);
            node->energy += 3; 
            node->insertions++;
            map[cid] = node;

            std::cerr << "Linked circle id=" << cid << " to existing node id=" << pid << std::endl;

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
    std::cerr << "Finding best edge to insert via segment R-tree..." << std::endl;
    constexpr int K = 16;
    std::vector<SegmentValue> sres;
    rtreeS.query(bgi::nearest(center, K), std::back_inserter(sres));
    double bestAdd = std::numeric_limits<double>::infinity();
    TourNode* bestU = nullptr;
    Point bestPt;

    for (auto& sv : sres) {
        int uid = sv.second;
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

    std::cerr << "Best edge for insertion: between id=" << bestU->id << " and id=" << bestU->next->id
              << " at point (" << bg::get<0>(bestPt) << ", " << bg::get<1>(bestPt) << "), addCost=" << bestAdd << std::endl;

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

    std::cerr << "Inserted new TourNode id=" << cid << " between id=" << L->id << " and id=" << R->id << std::endl;

    // Process the neighbors
    processNeighbor(node->prev);
    node = map[cid]; // The delete-reinsert cascade might have invalidated this pointer
    processNeighbor(node->next);
    
    node = map[cid]; // Re-fetch the node after potential deletion
    if (node->insertions >= 2 && (node->insertions & (node->insertions - 1)) == 0) {
        optimizeTourPoint(node, head, rtreeP, rtreeS, map, allNodes);
    }
}

// Given a merge tree (vector<TreeNode>), reconstruct the tour by "unmerging" nodes.
// At each step, remove the highest-weight internal node from the tour and insert its children.
// Returns the sequence of tour points.
std::vector<Point> unmerge(std::vector<TreeNode> treeNodes) {
    int N = treeNodes.size();

    std::cerr << "Starting unmerge with " << N << " tree nodes." << std::endl;

    // Priority queue to process nodes in descending order of weight (merge cost)
    std::priority_queue<std::pair<double, int>> pq;
    pq.push({treeNodes[N - 1].weight, N - 1}); // Start from the root

    TourNode* head = nullptr; // Head of the cyclic tour
    bgi::rtree<PointValue, bgi::rstar<16>> rtreeP; // R-tree for tour points
    bgi::rtree<SegmentValue, bgi::rstar<16>> rtreeS; // R-tree for tour segments
    std::vector<TourNode*> map(N, nullptr); // Map from node id to TourNode*

    // Insert the root node's circle into the tour
    std::cerr << "Inserting root node id=" << (N - 1) << std::endl;
    insertCircle(N - 1, treeNodes[N - 1], head, rtreeP, rtreeS, map, treeNodes);

    // Process nodes in order of decreasing weight
    int nodesProcessed = 0;
    while (!pq.empty()) {
        auto [w, idx] = pq.top();
        pq.pop();
        nodesProcessed++;
        if (treeNodes[idx].left < 0) continue; // Skip leaves

        std::cerr << "Unmerging node id=" << idx << " (weight=" << w << "), left=" << treeNodes[idx].left << ", right=" << treeNodes[idx].right << std::endl;

        // Remove the current node's circle from the tour
        TourNode* tnode = map[idx];
        if (tnode) {
            auto& vc = tnode->circles;
            // Remove this node's id from the circle list
            vc.erase(std::remove(vc.begin(), vc.end(), idx), vc.end());
            map[idx] = nullptr;

            if (vc.empty()) {
                std::cerr << "Deleting TourNode id=" << idx << " during unmerge." << std::endl;
                deleteTourNode(head, tnode, rtreeP, rtreeS, map);
            } else if (idx == tnode->id) {
                // We removed the main id; need to update id and R-trees
                int oldId = tnode->id;
                int newId = vc.front();

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

                std::cerr << "Updated TourNode id from " << oldId << " to " << newId << " and updated R-trees." << std::endl;
            }
        }

        // Insert the left and right children into the tour
        for (int c : {treeNodes[idx].left, treeNodes[idx].right}) {
            std::cerr << "Inserting child node id=" << c << std::endl;
            insertCircle(c, treeNodes[c], head, rtreeP, rtreeS, map, treeNodes);
            // If the child is not a leaf, add it to the queue for further unmerging
            if (treeNodes[c].left >= 0) {
                pq.push({treeNodes[c].weight, c});
            }
        }

        // Once in a while, optimize all nodes in the tour
        if ((nodesProcessed >= 2 && (nodesProcessed & (nodesProcessed - 1)) == 0) || pq.empty()) {
            // Optimize all points in the tour
            for (int nid = N - 1; nid >= 0; nid--) {
                if (map[nid] == nullptr) continue; // Skip if not in the tour
                TourNode* t = map[nid];
                if (t) {
                    t->insertions++;
                    if (t->insertions >= 2 && (t->insertions & (t->insertions - 1)) == 0) {
                        std::cerr << "Optimizing TourNode id=" << t->id << " at (" << bg::get<0>(t->pos) << ", " << bg::get<1>(t->pos) << ")" << std::endl;
                        optimizeTourPoint(t, head, rtreeP, rtreeS, map, treeNodes);
                    }
                }
            }
            
            // Interleave optimization with energy reduction
            // For every node in the pq, decrease energy and handle zero-energy nodes
            for (int nid = N - 1; nid >= 0; nid--) {
                std::cerr << "Processing node id=" << nid << " for energy reduction." << std::endl;
                if (map[nid] == nullptr) continue; // Skip if not in the tour
                TourNode* t = map[nid];
                if (t) {
                    t->energy--;
                    if (t->energy == 0) {
                        auto linked = t->circles;
                        deleteTourNode(head, t, rtreeP, rtreeS, map);
                        std::sort(linked.begin(), linked.end(), [&](int a, int b) { return treeNodes[a].r < treeNodes[b].r; });
                        for (int oc : linked) {
                            insertCircle(oc, treeNodes[oc], head, rtreeP, rtreeS, map, treeNodes);
                        }
                    }
                }
            }

            // Optimize all points in the tour
            for (int nid = N - 1; nid >= 0; nid--) {
                if (map[nid] == nullptr) continue; // Skip if not in the tour
                TourNode* t = map[nid];
                if (t) {
                    t->insertions++;
                    if (t->insertions >= 2 && (t->insertions & (t->insertions - 1)) == 0) {
                        std::cerr << "Optimizing TourNode id=" << t->id << " at (" << bg::get<0>(t->pos) << ", " << bg::get<1>(t->pos) << ")" << std::endl;
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

    std::cerr << "Tour reconstruction complete. Tour size: " << tour.size() << std::endl;

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

// Verify that the tour is valid: each circle should contain a point in the tour.
bool verifyTour(const std::vector<Point>& tour, const std::vector<Circle>& circles) {
    // Construct an R*-tree from the tour points
    bgi::rtree<Point, bgi::rstar<16>> rtree;
    for (const auto& p : tour) {
        rtree.insert(p);
    }

    // Verify each circle contains at least one point from the tour
    for (const auto& c : circles) {
        std::vector<Point> candidates;
        rtree.query(bgi::within(circleBox(c)), std::back_inserter(candidates));
        if (candidates.empty()) {
            std::cerr << "Circle with center (" << bg::get<0>(c.center) << ", " << bg::get<1>(c.center) 
                      << ") and radius " << c.r << " does not contain any tour points." << std::endl;
            return false;
        }
    }

    return true;
}

double totalTourDistance(const std::vector<Point>& tour) {
    if (tour.empty()) return 0.0;

    double totalDist = 0.0;
    for (size_t i = 0; i < tour.size(); ++i) {
        const Point& p1 = tour[i];
        const Point& p2 = tour[(i + 1) % tour.size()]; // Wrap around to the start
        totalDist += bg::distance(p1, p2);
    }
    return totalDist;
}


std::vector<Point> CETSP(const std::vector<Circle>& circles) {
    auto circlesCopy = circles; // Make a copy to avoid modifying the input

    // Initialize the best tour and distance
    double bestDist = std::numeric_limits<double>::infinity();
    std::vector<Point> bestTour;

    totInserts = 0; // Reset global insert counter

    // Step 1: Remove covering circles to simplify the problem
    removeCoveringCircles(circlesCopy);

    // If there's only one circle left, directly return its center as the tour
    if (circlesCopy.size() == 1) {
        std::cerr << "Only one circle remains after removing coverings. Returning its center." << std::endl;
        return {circlesCopy[0].center};
    }

    // Repeat for some number of iterations
    constexpr int REPEAT_ITERS = 10;
    for (int iter = 0; iter < REPEAT_ITERS; ++iter) {
        // Step 2: Build the merge tree from the remaining circles
        std::vector<TreeNode> mergeTree = buildMergeTree(circlesCopy);

        // Step 3: Generate the tour from the merge tree
        std::vector<Point> tour = unmerge(mergeTree);

        // Step 4: Verify the generated tour to ensure correctness and update the best tour if valid
        if (verifyTour(tour, circles)) {
            double curDist = totalTourDistance(tour);
            if (curDist < bestDist) {
                bestDist = curDist;
                bestTour = tour;
            }
        }
    }

    if (bestTour.empty()) {
        std::cout << "Tour verification failed." << std::endl;
    }
    else {
        std::cout << "Tour verification passed. Number of elements: " << circles.size() << ", Total inserts: " << totInserts << std::endl;
    }

    // Return the valid tour
    return bestTour;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::cerr.rdbuf(nullptr); // Uncomment to disable debug 
    // std::ofstream debugFile("./debug.txt");
    // std::cerr.rdbuf(debugFile.rdbuf());

    // Prepare input and output directories
    namespace fs = std::filesystem;
    const fs::path inputDir("./data/cmu_processed");
    const fs::path outputDir("./data/cmu_output");
    if (!fs::exists(outputDir)) {
        fs::create_directory(outputDir);
    }
    std::cerr << "Input directory: " << inputDir << std::endl;
    std::cerr << "Output directory: " << outputDir << std::endl;

    // List of specific files to read (empty means read everything)
    std::vector<std::string> filesToRead = {}; // Add filenames here if needed

    // Iterate through all files in the input directory
    for (const auto& entry : fs::directory_iterator(inputDir)) {
        if (!filesToRead.empty()) {
            // Check if the current file is in the list
            if (std::find(filesToRead.begin(), filesToRead.end(), entry.path().filename().string()) == filesToRead.end()) {
                continue;
            }
        }
        if (!entry.is_regular_file()) continue;

        const fs::path inputFile = entry.path();
        const fs::path outputFile = outputDir / inputFile.filename();

        std::cout << "Processing file: " << inputFile << std::endl;

        auto startTime = std::chrono::high_resolution_clock::now();

        // Read circles from the input file
        std::vector<Circle> circles;
        std::ifstream inFile(inputFile);
        if (!inFile) {
            std::cout << "Failed to open input file: " << inputFile << std::endl;
            continue;
        }

        double x, y, r;
        int nodeId = 0;
        while (inFile >> x >> y >> r) {
            circles.push_back({Point(x, y), r, nodeId++, -1, 0.0, {}});
        }
        inFile.close();

        // Generate the tour using the CETSP algorithm
        auto tour = CETSP(circles);
        if (tour.empty()) {
            std::cout << "Failed to generate a valid tour for file: " << inputFile << std::endl;
            continue;
        }

        // Write the tour points to the output file
        std::ofstream outFile(outputFile);
        if (!outFile) {
            std::cout << "Failed to open output file: " << outputFile << std::endl;
            continue;
        }

        outFile << std::fixed << std::setprecision(4);
        for (const auto& p : tour) {
            outFile << bg::get<0>(p) << " " << bg::get<1>(p) << "\n";
        }
        outFile.close();

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        std::cout << "Finished processing file: " << inputFile << " in " << duration << " ms" << std::endl;
    }

    return 0;
}

