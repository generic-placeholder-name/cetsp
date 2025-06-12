// Slow code
// Does not include several optimizations that I have done later
// Keep around since it helps understand the algorithm

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
typedef bg::model::box<Point> Box;
typedef std::pair<Box, int> RTreeValue;

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
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(overlapDepth, h);
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

    bgi::rtree<RTreeValue, bgi::rstar<16>> rtree_pr;
    std::vector<bool> removed(n, false);
    for (int id : order) {
        if (removed[id]) continue;
        std::vector<RTreeValue> candidates;
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
    bgi::rtree<RTreeValue, bgi::rstar<16>> rtree_kd;
    for (auto& c : circles) {
        rtree_kd.insert(std::make_pair(circleBox(c), c.nodeId));
    }

    // 2) Dynamic NN maintenance: each circle stores its nn and gap, and we keep a global min-heap
    std::set<std::pair<double, int>> heap;
    std::vector<std::pair<double, int>> handles(2 * n, {-1, -1});
    for (auto& c : circles) {
        // find nearest neighbor via two-NN query
        std::vector<RTreeValue> res;
        rtree_kd.query(bgi::nearest(c.center, 2), std::back_inserter(res));
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
        std::vector<RTreeValue> resn;
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

            std::vector<RTreeValue> resc;
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
    std::cout << "Combine complete. Total tree nodes=" << treeNodes.size() << std::endl;

    return treeNodes;
}

// A TourNode represents one point in the (cyclic) tour.
// - pos = its coordinates.
// - circles = list of all TreeNode‐indices currently "linked" here.
// - energy = integer energy.  If it ever hits 0, we delete this TourNode,
//   re‐insert all its circles (in ascending‐radius order), and remove the node.
// - prev/next pointers form a circular doubly‐linked list.
struct TourNode {
    Point pos;
    std::vector<int> circles;
    int energy;
    TourNode* prev;
    TourNode* next;
    explicit TourNode(const Point& p)
        : pos(p), energy(0), prev(nullptr), next(nullptr) {}
};

// Given a circle (center, radius) and an edge [p1,p2], returns the point
// on or within the circle that minimizes insertion‐cost between p1 and p2.
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
    double phi = 0.5 * (alpha + beta);

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

// Delete a TourNode from the cyclic list.  Also clears its mapping from
// `circleToTourNode` so that any circle index previously linked there is now nullptr.
void deleteTourNode(TourNode*& head,
                    TourNode* node,
                    std::vector<TourNode*>& circleToTourNode) {
    // Remove circle‐to‐node mappings
    for (int cid : node->circles) {
        circleToTourNode[cid] = nullptr;
    }
    // Splice out of doubly linked list
    if (node->prev) node->prev->next = node->next;
    if (node->next) node->next->prev = node->prev;
    if (head == node) {
        // If this was the head, advance head (unless it’s the only node)
        head = (node->next != node) ? node->next : nullptr;
    }
    delete node;
}

int totInserts = 0; // Global counter for debugging
// Insert the circle `cid` (with data in `tn`) into the current tour starting at `head`.
// Updates energies, and if any neighbor hits energy=0, recursively removes it and
// reinserts its linked circles in ascending‐radius order.
// `circleToTourNode[cid]` will be set to wherever the circle ultimately lives.
void insertCircle(int cid,
                  const TreeNode& tn,
                  TourNode*& head,
                  std::vector<TourNode*>& circleToTourNode,
                  const std::vector<TreeNode>& allTreeNodes) {
    const Point& cpos = tn.center;
    double cr = tn.r;
    totInserts++;

    std::cerr << "Inserting circle " << cid << " with center (" 
              << bg::get<0>(cpos) << ", " << bg::get<1>(cpos) 
              << ") and radius " << cr << std::endl;

    // Case 1: empty tour → create first node
    if (!head) {
        std::cerr << "Tour is empty. Creating first node." << std::endl;
        TourNode* first = new TourNode(cpos);
        first->circles.push_back(cid);
        first->energy = 0;
        first->prev = first->next = first;
        head = first;
        circleToTourNode[cid] = first;
        std::cerr << "Inserted circle " << cid << " as the first node." << std::endl;
        return;
    }

    // 1A) Check if any existing TourNode is within distance <= cr
    TourNode* found = nullptr;
    for (TourNode* it = head;; it = it->next) {
        if (bg::distance(cpos, it->pos) <= cr) {
            found = it;
            break;
        }
        if (it->next == head) break;
    }
    if (found) {
        std::cerr << "Found existing TourNode within radius. Linking circle " << cid << " to it." << std::endl;
        // Link circle into that node
        found->circles.push_back(cid);
        found->energy += 3;
        // Decrement neighbors’ energies
        found->prev->energy -= 1;
        found->next->energy -= 1;
        circleToTourNode[cid] = found;

        // If any neighbor hits zero, remove and reinsert its circles
        std::vector<TourNode*> toRemove;
        if (found->prev->energy == 0) toRemove.push_back(found->prev);
        if (found->prev != found->next && found->next->energy == 0) toRemove.push_back(found->next);

        for (TourNode* rm : toRemove) {
            std::cerr << "Neighbor TourNode energy hit zero. Removing and reinserting its circles." << std::endl;
            // Extract their linked circles
            auto linked = rm->circles;
            rm->circles.clear();
            // Sort by radius ascending, then reinsert
            std::sort(
                linked.begin(),
                linked.end(),
                [&](int a, int b) {
                    return allTreeNodes[a].r < allTreeNodes[b].r;
                }
            );
            deleteTourNode(head, rm, circleToTourNode);
            for (int ocid : linked) {
                std::cerr << "Reinserting circle " << ocid << " from removed TourNode." << std::endl;
                insertCircle(ocid, allTreeNodes[ocid], head, circleToTourNode, allTreeNodes);
            }
        }
        return;
    }

    // 1B) Otherwise, find the best edge (u→v) to insert a new node
    double bestAdd = std::numeric_limits<double>::infinity();
    TourNode *bestA = nullptr, *bestB = nullptr;
    for (TourNode* u = head;; u = u->next) {
        TourNode* v = u->next;
        Point cand = findOptimalPoint(cpos, cr, u->pos, v->pos);
        double addCost = bg::distance(u->pos, cand)
                       + bg::distance(cand, v->pos)
                       - bg::distance(u->pos, v->pos);
        if (addCost < bestAdd) {
            bestAdd = addCost;
            bestA = u;
            bestB = v;
        }
        if (v == head) break;
    }
    std::cerr << "Best edge found for insertion between nodes with positions (" 
              << bg::get<0>(bestA->pos) << ", " << bg::get<1>(bestA->pos) << ") and ("
              << bg::get<0>(bestB->pos) << ", " << bg::get<1>(bestB->pos) << ")." << std::endl;

    // Compute the final insertion point
    Point insPt = findOptimalPoint(cpos, cr, bestA->pos, bestB->pos);
    std::cerr << "Computed insertion point: (" << bg::get<0>(insPt) << ", " << bg::get<1>(insPt) << ")." << std::endl;

    // Create a new TourNode there
    TourNode* node = new TourNode(insPt);
    node->circles.push_back(cid);
    node->energy = 0;
    // Link into the cyclic list between bestA and bestB
    node->prev = bestA;
    node->next = bestB;
    bestA->next = node;
    bestB->prev = node;
    circleToTourNode[cid] = node;
    std::cerr << "Inserted new TourNode for circle " << cid << " into the tour." << std::endl;

    // Energy updates
    node->energy += 3;
    node->prev->energy -= 1;
    node->next->energy -= 1;

    // If neighbor hits zero, remove+reinsert
    std::vector<TourNode*> toDel;
    if (node->prev->energy == 0) toDel.push_back(node->prev);
    if (node->prev != node->next && node->next->energy == 0) toDel.push_back(node->next);
    for (TourNode* rm : toDel) {
        std::cerr << "Neighbor TourNode energy hit zero. Removing and reinserting its circles." << std::endl;
        auto linked = rm->circles;
        rm->circles.clear();
        std::sort(
            linked.begin(),
            linked.end(),
            [&](int a, int b) {
                return allTreeNodes[a].r < allTreeNodes[b].r;
            }
        );
        deleteTourNode(head, rm, circleToTourNode);
        for (int ocid : linked) {
            std::cerr << "Reinserting circle " << ocid << " from removed TourNode." << std::endl;
            insertCircle(ocid, allTreeNodes[ocid], head, circleToTourNode, allTreeNodes);
        }
    }
}

// Given a Combine-tree (vector<TreeNode>), produce the final tour (vector<Point>).
// TreeNodes are indexed 0..N−1, with the root at index N−1.  Every node has its own
// (center, r), and left/right children (−1/−1 for leaves).  Merge‐weight is in weight.
// We "unmerge" by always extracting the highest‐weight internal node, removing its circle
// from the tour, then inserting its two children, etc.
std::vector<Point> unmerge(std::vector<TreeNode> treeNodes) {
    int N = (int)treeNodes.size();
    if (N == 0) return {};

    std::cerr << "Starting unmerge process with " << N << " tree nodes." << std::endl;

    // Build a max‐heap keyed by weight, storing (weight, index).
    std::priority_queue<std::pair<double, int>> pq;
    int root = N - 1;
    pq.push({treeNodes[root].weight, root});

    std::cerr << "Root node index: " << root << ", weight: " << treeNodes[root].weight << std::endl;

    // `head` is our cyclic tour’s starting pointer.
    TourNode* head = nullptr;
    // map circle‐index ⇒ TourNode* (nullptr if not in tour)
    std::vector<TourNode*> circleToTourNode(N, nullptr);

    // 1) Start by inserting the root circle into an initially‐empty tour.
    insertCircle(root, treeNodes[root], head, circleToTourNode, treeNodes);
    std::cerr << "Inserted root circle into the tour." << std::endl;

    // 2) Pop until the queue is empty.
    while (!pq.empty()) {
        auto [w, idx] = pq.top();
        pq.pop();

        std::cerr << "Processing node index: " << idx << ", weight: " << w << std::endl;

        // If this is a leaf, skip.  Only internal nodes split into children.
        if (treeNodes[idx].left == -1) {
            std::cerr << "Node " << idx << " is a leaf. Skipping." << std::endl;
            continue;
        }

        // Remove circle `idx` from whichever TourNode it currently lives in:
        TourNode* tnode = circleToTourNode[idx];
        if (tnode) {
            std::cerr << "Removing circle " << idx << " from its current TourNode." << std::endl;
            // Erase idx from tnode→circles
            auto& vec = tnode->circles;
            vec.erase(std::remove(vec.begin(), vec.end(), idx), vec.end());
            circleToTourNode[idx] = nullptr;
            // If that TourNode is now empty, delete it
            if (vec.empty()) {
                std::cerr << "TourNode is now empty. Deleting it." << std::endl;
                deleteTourNode(head, tnode, circleToTourNode);
            }
        }

        // Insert left and right children:
        int c1 = treeNodes[idx].left;
        int c2 = treeNodes[idx].right;
        std::cerr << "Inserting children of node " << idx << ": left=" << c1 << ", right=" << c2 << std::endl;
        insertCircle(c1, treeNodes[c1], head, circleToTourNode, treeNodes);
        insertCircle(c2, treeNodes[c2], head, circleToTourNode, treeNodes);

        // If those children are themselves internal, push them onto pq:
        if (treeNodes[c1].left != -1) {
            std::cerr << "Child " << c1 << " is internal. Pushing onto priority queue." << std::endl;
            pq.push({treeNodes[c1].weight, c1});
        }
        if (treeNodes[c2].left != -1) {
            std::cerr << "Child " << c2 << " is internal. Pushing onto priority queue." << std::endl;
            pq.push({treeNodes[c2].weight, c2});
        }

        // Debug-print the current tour
        std::cerr << "Current tour after processing node " << idx << ":" << std::endl;
        if (head) {
            TourNode* curr = head;
            do {
                std::cerr << "(" << bg::get<0>(curr->pos) << ", " << bg::get<1>(curr->pos) << ") ";
                curr = curr->next;
            } while (curr != head);
            std::cerr << std::endl;
        } else {
            std::cerr << "Tour is empty." << std::endl;
        }
    }

    // Finally, walk the cyclic tour once to collect points in order.
    std::vector<Point> tour;
    if (!head) {
        std::cerr << "Tour is empty." << std::endl;
        return tour;
    }
    std::cerr << "Collecting points from the cyclic tour." << std::endl;
    TourNode* curr = head;
    do {
        tour.push_back(curr->pos);
        curr = curr->next;
    } while (curr != head);

    std::cout << "Final tour collected with " << tour.size() << " points." << std::endl;

    // Delete all TourNodes to free memory
    while (head) {
        deleteTourNode(head, head, circleToTourNode);
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
            if (totalTourDistance(tour) < bestDist) {
                bestDist = totalTourDistance(tour);
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

    std::cerr.rdbuf(nullptr); // Uncomment to disable debug output

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

