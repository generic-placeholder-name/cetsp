#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <queue>
#include <cmath>
#include <list>
#include <random>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/intrusive/list.hpp>
#include <chrono> 

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bi = boost::intrusive;
namespace fs = std::filesystem;

using Point = bg::model::point<double, 2, bg::cs::cartesian>;
using RTree = bgi::rtree<Point, bgi::rstar<16>>;
using IdxRTree = bgi::rtree<std::pair<Point, size_t>, bgi::rstar<16>>;
using Box = bg::model::box<Point>;
using BoxIdxRTree = bgi::rtree<std::pair<Box, size_t>, bgi::rstar<16>>;
using Segment = bg::model::segment<Point>;
using SegmentRTree = bgi::rtree<Segment, bgi::rstar<16>>;
using SegmentIdxRTree = bgi::rtree<std::tuple<Segment, size_t, size_t>, bgi::rstar<16>>;

// Function to calculate the midpoint between two points
Point midpointBetween(const Point& a, const Point& b) {
    return Point((bg::get<0>(a) + bg::get<0>(b)) / 2.0,
                 (bg::get<1>(a) + bg::get<1>(b)) / 2.0);
}

Point findOptimalPoint(const Point& center, double radius, const Point& p1, const Point& p2) {
    constexpr double RAD_EPSILON = 1e-12;

    // If radius is zero, return the center
    if (radius < RAD_EPSILON) {
        return center;
    }

    // Calculate the vector from p1 to p2
    double dx = bg::get<0>(p2) - bg::get<0>(p1);
    double dy = bg::get<1>(p2) - bg::get<1>(p1);

    // Handle the case where p1 and p2 are the same
    if (std::abs(dx) < 1e-9 && std::abs(dy) < 1e-9) {
        double dist = bg::distance(center, p1);
        if (dist <= radius) {
            return p1; // p1 is already optimal
        } else {
            // Closest point on the circle to p1
            double scale = radius / dist - RAD_EPSILON;
            return Point(bg::get<0>(center) + scale * (bg::get<0>(p1) - bg::get<0>(center)),
                         bg::get<1>(center) + scale * (bg::get<1>(p1) - bg::get<1>(center)));
        }
    }

    // Normalize the vector
    double magnitude = std::sqrt(dx * dx + dy * dy);
    dx /= magnitude;
    dy /= magnitude;

    // Project the center onto the infinite line defined by p1 and p2
    double t = ((bg::get<0>(center) - bg::get<0>(p1)) * dx +
                (bg::get<1>(center) - bg::get<1>(p1)) * dy);

    // Clamp the projection to the segment [p1, p2]
    double t_clamped = std::max(0.0, std::min(t, magnitude));
    Point clamped_projection(bg::get<0>(p1) + t_clamped * dx,
                              bg::get<1>(p1) + t_clamped * dy);

    // Check if the clamped projection is within the circle
    double distToProjection = bg::distance(center, clamped_projection);
    if (distToProjection <= radius) {
        return clamped_projection;
    } else {
        // Move the clamped projection to the circle's boundary
        double scale = radius / distToProjection - RAD_EPSILON;
        return Point(bg::get<0>(center) + scale * (bg::get<0>(clamped_projection) - bg::get<0>(center)),
                     bg::get<1>(center) + scale * (bg::get<1>(clamped_projection) - bg::get<1>(center)));
    }
}


struct TSPNode {
    Point basePoint, tourPoint;
    double distance;

    // Constructor with basePoint and distance, tourPoint defaults to basePoint
    TSPNode(const Point& basePoint, double distance = 0)
        : basePoint(basePoint),
          tourPoint(basePoint),
          distance(distance) {}

    // Constructor with basePoint, distance, and tourPoint explicitly specified
    TSPNode(const Point& basePoint, double distance, const Point& tourPoint)
        : basePoint(basePoint),
          tourPoint(tourPoint),
          distance(distance) {}

    // Function to perturb the tourPoint by a random point within distance * ratio
    void perturb(double ratio = 0.1) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dist(-1.0, 1.0);

        Point newPoint;
        do {
            double dx = dist(gen) * distance * ratio;
            double dy = dist(gen) * distance * ratio;
            newPoint = Point(bg::get<0>(tourPoint) + dx, bg::get<1>(tourPoint) + dy);
        } while (bg::distance(basePoint, newPoint) > distance);

        tourPoint = newPoint;
    }
};

// Function to get the previous iterator in a circular list
template <typename Iterator, typename List>
Iterator prevIt(Iterator it, List& list) {
    auto prv = (it == list.begin()) ? std::prev(list.end()) : std::prev(it);
    return prv;
}

// Function to get the next iterator in a circular list
template <typename Iterator, typename List>
Iterator nextIt(Iterator it, List& list) {
    auto nxt = std::next(it); 
    return nxt == list.end() ? list.begin() : nxt; 
}

double calculateTotalTourDistance(const std::vector<TSPNode>& tour) {
    if (tour.empty()) return 0.0;

    double totalDistance = 0.0;

    for (size_t i = 0; i < tour.size(); i++) {
        // Get next node
        const auto& nextNode = tour[i == tour.size() - 1 ? 0 : i + 1];

        // Calculate distance between consecutive nodes
        const auto& currentPoint = tour[i].tourPoint;
        const auto& nextPoint = nextNode.tourPoint;
        totalDistance += bg::distance(currentPoint, nextPoint);
    }

    return totalDistance;
}


// Struct to hold the distance between two points for the priority queue
struct DistanceEntry {
    double distance;
    size_t idxA, idxB;

    bool operator<(const DistanceEntry& other) const {
        return distance > other.distance;
    }
};

struct MergeNode {
    double distance;
    size_t idx;  // Index of the point in the `points` vector
    MergeNode* a;
    MergeNode* b;

    MergeNode(double dist, size_t index, MergeNode* childA = nullptr, MergeNode* childB = nullptr)
        : distance(dist), idx(index), a(childA), b(childB) {}

    ~MergeNode() {
        delete a;
        delete b;
    }

    bool isLeaf() const {
        return !a && !b;
    }
};

struct UnmergeNode {
    Point point; 
    std::unordered_set<size_t> indexes; 
    size_t capacity; 

    UnmergeNode(const Point& point, size_t idx) : point(point), capacity(1) {indexes.insert(idx);}
    UnmergeNode(const Point& point, std::unordered_set<size_t> indexes) : point(point), indexes(indexes), capacity(indexes.size()) {}
    UnmergeNode(const Point& point, std::vector<size_t> indexes) : point(point), indexes(indexes.begin(), indexes.end()), capacity(indexes.size()) {}
    UnmergeNode(const Point& point, std::unordered_set<size_t> indexes, size_t capacity) : point(point), indexes(indexes), capacity(capacity) {}
    UnmergeNode(const Point& point, std::vector<size_t> indexes, size_t capacity) : point(point), indexes(indexes.begin(), indexes.end()), capacity(capacity) {}
};

// TSP approximation
std::vector<size_t> approximateTSP(std::vector<Point> points) {
    IdxRTree tree;
    for (size_t i = 0; i < points.size(); i++) {
        tree.insert({points[i], i});
    }

    std::priority_queue<DistanceEntry> distanceQueue;
    std::vector<MergeNode*> nodes(points.size(), nullptr);
    std::vector<bool> pointExists(points.size(), true);

    // Initialize distances and priority queue
    std::vector<std::vector<size_t>> closestIdx(points.size()); // Points that have point i as their closest point

    for (size_t i = 0; i < points.size(); i++) {
        std::vector<std::pair<Point, size_t>> result;
        tree.query(bgi::nearest(points[i], 2), std::back_inserter(result)); // We query two points since the query also returns point i

        for (auto& point : result) {
            if (point.second == i) {
                continue; // R-Tree query returns the point itself, we do not process this
            }

            double dist = bg::distance(points[i], point.first);
            size_t nearestIdx = point.second;
            distanceQueue.push({dist, i, nearestIdx});
            closestIdx[nearestIdx].push_back(i);
        }

        // Create leaf nodes
        nodes[i] = new MergeNode(0.0, i);
    }

    // Iteratively pop pairs of closest points and replace with midpoints
    while (!distanceQueue.empty()) {
        auto [d, idxA, idxB] = distanceQueue.top();
        distanceQueue.pop();

        if (!pointExists[idxA] || !pointExists[idxB]) continue;

        Point midpoint = midpointBetween(points[idxA], points[idxB]);
        size_t newIdx = points.size();
        points.push_back(midpoint);
        pointExists.push_back(true);
        closestIdx.push_back({});

        tree.insert({midpoint, newIdx});
        tree.remove({points[idxA], idxA});
        tree.remove({points[idxB], idxB});
        pointExists[idxA] = pointExists[idxB] = false;

        // Create a new node for the midpoint
        nodes.push_back(new MergeNode(d, newIdx, nodes[idxA], nodes[idxB]));

        // Update nearest neighbors and distances for affected points
        std::vector<size_t> affectedPoints = closestIdx[idxA];
        affectedPoints.insert(affectedPoints.end(), closestIdx[idxB].begin(), closestIdx[idxB].end());
        affectedPoints.push_back(newIdx);

        closestIdx[idxA].clear();
        closestIdx[idxB].clear();

        for (size_t affectedIdx : affectedPoints) {
            if (!pointExists[affectedIdx]) continue;

            std::vector<std::pair<Point, size_t>> result;
            tree.query(bgi::nearest(points[affectedIdx], 2), std::back_inserter(result));

            for (auto& point : result) {
                if (point.second == affectedIdx) {
                    continue; // R-Tree query returns the point itself, we do not process this
                }

                double dist = bg::distance(points[affectedIdx], point.first);
                size_t nearestIdx = point.second;
                distanceQueue.push({dist, affectedIdx, nearestIdx});
                closestIdx[nearestIdx].push_back(affectedIdx);
            }
        }
    }

    // Store points in a circular linked list
    std::list<size_t> tour;
    std::vector<std::list<size_t>::iterator> tourNodes(points.size());
    tourNodes[points.size() - 1] = tour.insert(tour.end(), points.size() - 1);

    // Decoding phase
    auto nodeComparator = [&nodes](size_t lhs, size_t rhs) {
        return nodes[lhs]->distance < nodes[rhs]->distance; // Max-heap based on distance
    };

    std::priority_queue<size_t, std::vector<size_t>, decltype(nodeComparator)> maxHeap(nodeComparator);
    maxHeap.push(points.size() - 1);

    SegmentIdxRTree segTree;
    auto makeSeg = [&](size_t i, size_t j) {
        if (i > j) {
            std::swap(i, j);
        }
        return std::make_tuple(Segment(points[i], points[j]), i, j);
    };

    while (!maxHeap.empty()) {
        auto maxIdx = maxHeap.top();
        maxHeap.pop();

        auto maxIdxIt = tourNodes[maxIdx];
        MergeNode* maxNode = nodes[maxIdx];

        if (!maxNode->isLeaf()) {
            // Retrieve indices of the left and right children
            size_t leftIdx = maxNode->a->idx, rightIdx = maxNode->b->idx;
            bool leftIsLeaf = maxNode->a->isLeaf(), rightIsLeaf = maxNode->b->isLeaf();

            // Erase the current node
            if (tour.size() > 1) {
                auto prevItr = prevIt(maxIdxIt, tour);
                auto nextItr = nextIt(maxIdxIt, tour);
                segTree.remove(makeSeg(maxIdx, *prevItr));
                if (tour.size() > 2) {
                    segTree.remove(makeSeg(maxIdx, *nextItr));
                    segTree.insert(makeSeg(*prevItr, *nextItr));
                }
            }
            tour.erase(maxIdxIt);

            if (tour.empty()) {
                tourNodes[leftIdx] = tour.insert(tour.end(), leftIdx);
                tourNodes[rightIdx] = tour.insert(tour.end(), rightIdx);
                segTree.insert(makeSeg(leftIdx, rightIdx));
            } else if (tour.size() == 1) {
                auto itr = tour.begin();
                tourNodes[leftIdx] = tour.insert(tour.end(), leftIdx);
                tourNodes[rightIdx] = tour.insert(tour.end(), rightIdx);
                segTree.insert(makeSeg(leftIdx, rightIdx));
                segTree.insert(makeSeg(leftIdx, *itr));
                segTree.insert(makeSeg(*itr, rightIdx));
            } else {
                auto getBestSeg = [&](const Point& point) {
                    static constexpr size_t K = 26;
                    std::vector<std::tuple<Segment, size_t, size_t>> candidates;
                    segTree.query(bgi::nearest(point, K), std::back_inserter(candidates));
                    double minExtraDist = std::numeric_limits<double>::max();
                    std::tuple<Segment, size_t, size_t> minCandidate;
                    for (const auto& candidate : candidates) {
                        const auto& [p1, p2] = std::get<0>(candidate);
                        double extraDist = bg::distance(p1, point) + bg::distance(p2, point) - bg::distance(p1, p2);
                        if (extraDist < minExtraDist) {
                            minExtraDist = extraDist;
                            minCandidate = candidate;
                        }
                    }
                    return minCandidate;
                };

                auto insertIdx = [&](size_t idx) {
                    auto bestSeg = getBestSeg(points[idx]);
                    auto [seg, prevIdx, nextIdx] = bestSeg;
                    auto [prevPoint, nextPoint] = seg;
                    if (nextIt(tourNodes[prevIdx], tour) != tourNodes[nextIdx]) {
                        std::swap(prevIdx, nextIdx);
                        std::swap(prevPoint, nextPoint);
                        assert(nextIt(tourNodes[prevIdx], tour) == tourNodes[nextIdx]);
                    } 
                    tourNodes[idx] = tour.insert(tourNodes[nextIdx], idx);
                    segTree.remove(bestSeg);
                    segTree.insert(makeSeg(prevIdx, idx));
                    segTree.insert(makeSeg(idx, nextIdx));
                };
                
                insertIdx(leftIdx);
                insertIdx(rightIdx);
            }

            // Push children to the heap if they are not leaves
            if (!leftIsLeaf) {
                maxHeap.push(leftIdx);
            }
            if (!rightIsLeaf) {
                maxHeap.push(rightIdx);
            }
        }
    }

    // Convert back to vector for convenience
    std::vector<size_t> ans;
    ans.reserve(points.size());
    for (const auto& node: tour) {
        ans.push_back(node);
    }

    return ans;
}

// Function to query the closest circle to a point
template <size_t K = 16>
std::pair<Box, size_t> getClosestCircle(const BoxIdxRTree& tree, const Point& point, size_t idx = std::numeric_limits<size_t>::max()) {
    std::vector<std::pair<Box, size_t>> nearestBoxes;
    tree.query(bgi::nearest(point, K), std::back_inserter(nearestBoxes));

    std::pair<Box, size_t> closestBox;
    double minDistance = std::numeric_limits<double>::max();

    for (const auto& [box, i] : nearestBoxes) {
        if (idx == i) {
            continue;
        }

        Point center((bg::get<0>(box.min_corner()) + bg::get<0>(box.max_corner())) / 2.0,
                     (bg::get<1>(box.min_corner()) + bg::get<1>(box.max_corner())) / 2.0);
        double radius = (bg::get<0>(box.max_corner()) - bg::get<0>(box.min_corner())) / 2.0;
        double distance = bg::distance(point, center) - radius;

        if (distance < minDistance) {
            minDistance = distance;
            closestBox = {box, i};
        }
    }


    return closestBox;
}

template<size_t K = 16>
Segment findBestSegment(const SegmentRTree& tree, const Point& point) {
    std::vector<Segment> candidates;
    tree.query(bgi::nearest(point, K), std::back_inserter(candidates));

    Segment minCandidate;
    double minExtraDist = std::numeric_limits<double>::max();

    for (const auto& seg : candidates) {
        const auto& [p1, p2] = seg;
        double extraDist = bg::distance(p1, point) + bg::distance(point, p2) - bg::distance(p1, p2);
        if (extraDist < minExtraDist) {
            minExtraDist = extraDist;
            minCandidate = seg;
        }
    }

    return minCandidate;
}

void rotatePoints(std::vector<Point>& points, double angle) {
    double cosTheta = std::cos(angle);
    double sinTheta = std::sin(angle);

    for (auto& point : points) {
        double newX = bg::get<0>(point) * cosTheta - bg::get<1>(point) * sinTheta;
        double newY = bg::get<0>(point) * sinTheta + bg::get<1>(point) * cosTheta;
        point.set<0>(newX);
        point.set<1>(newY);
    }
}

double getRandomAngle() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(0, std::numbers::pi); // Range: [0, pi]
    return dis(gen);
}

// CETSP approximation function
std::vector<TSPNode> makePairCenterCETSP(std::vector<Point> points, std::vector<double> distances) {
    auto start = std::chrono::high_resolution_clock::now();

    double rot = getRandomAngle();
    rotatePoints(points, rot);

    BoxIdxRTree tree;

    // Initialize R-tree with bounding boxes
    for (size_t i = 0; i < points.size(); i++) {
        Point minCorner(bg::get<0>(points[i]) - distances[i], bg::get<1>(points[i]) - distances[i]);
        Point maxCorner(bg::get<0>(points[i]) + distances[i], bg::get<1>(points[i]) + distances[i]);
        Box box(minCorner, maxCorner);
        tree.insert({box, i});
    }

    std::priority_queue<DistanceEntry> distanceQueue;
    std::vector<MergeNode*> nodes(points.size(), nullptr);
    std::vector<bool> pointExists(points.size(), true);

    // Initialize distances and priority queue
    std::vector<std::vector<size_t>> closestIdx(points.size());

    for (size_t i = 0; i < points.size(); i++) {

        auto [closestBox, nearestIdx] = getClosestCircle(tree, points[i], i);
        Point center((bg::get<0>(closestBox.min_corner()) + bg::get<0>(closestBox.max_corner())) / 2.0,
                     (bg::get<1>(closestBox.min_corner()) + bg::get<1>(closestBox.max_corner())) / 2.0);
        double radius = (bg::get<0>(closestBox.max_corner()) - bg::get<0>(closestBox.min_corner())) / 2.0;
        double dist = bg::distance(points[i], center) - radius;

        // std::cout << "Init " << i << " " << nearestIdx << std::endl;

        distanceQueue.push({dist, i, nearestIdx});
        closestIdx[nearestIdx].push_back(i);

        // Create leaf nodes
        nodes[i] = new MergeNode(0.0, i);
    }

    // Iteratively pop pairs of closest points and replace with midpoints
    while (!distanceQueue.empty()) {
        auto [d, idxA, idxB] = distanceQueue.top();
        distanceQueue.pop();

        if (!pointExists[idxA] || !pointExists[idxB]) continue;

        // std::cout << "Merging " << idxA << " " << idxB << std::endl;

        Point midpoint = midpointBetween(points[idxA], points[idxB]);
        double midRad = (distances[idxA] + distances[idxB]) / 2.0;
        size_t newIdx = points.size();
        points.push_back(midpoint);
        distances.push_back(midRad);
        pointExists.push_back(true);
        closestIdx.push_back({});

        Point minCorner(bg::get<0>(midpoint) - midRad, bg::get<1>(midpoint) - midRad);
        Point maxCorner(bg::get<0>(midpoint) + midRad, bg::get<1>(midpoint) + midRad);
        Box newBox(minCorner, maxCorner);

        tree.insert({newBox, newIdx});
        tree.remove({Box(Point(bg::get<0>(points[idxA]) - distances[idxA], bg::get<1>(points[idxA]) - distances[idxA]),
                         Point(bg::get<0>(points[idxA]) + distances[idxA], bg::get<1>(points[idxA]) + distances[idxA])),
                     idxA});
        tree.remove({Box(Point(bg::get<0>(points[idxB]) - distances[idxB], bg::get<1>(points[idxB]) - distances[idxB]),
                         Point(bg::get<0>(points[idxB]) + distances[idxB], bg::get<1>(points[idxB]) + distances[idxB])),
                     idxB});
        pointExists[idxA] = pointExists[idxB] = false;

        nodes.push_back(new MergeNode(d, newIdx, nodes[idxA], nodes[idxB]));

        std::vector<size_t> affectedPoints = closestIdx[idxA];
        affectedPoints.insert(affectedPoints.end(), closestIdx[idxB].begin(), closestIdx[idxB].end());
        affectedPoints.push_back(newIdx);

        closestIdx[idxA].clear();
        closestIdx[idxB].clear();

        for (size_t affectedIdx : affectedPoints) {
            if (!pointExists[affectedIdx]) continue;

            auto [closestBox, nearestIdx] = getClosestCircle(tree, points[affectedIdx], affectedIdx);
            Point center((bg::get<0>(closestBox.min_corner()) + bg::get<0>(closestBox.max_corner())) / 2.0,
                         (bg::get<1>(closestBox.min_corner()) + bg::get<1>(closestBox.max_corner())) / 2.0);
            double radius = (bg::get<0>(closestBox.max_corner()) - bg::get<0>(closestBox.min_corner())) / 2.0;
            double dist = bg::distance(points[affectedIdx], center) - radius;

            distanceQueue.push({dist, affectedIdx, nearestIdx});
            closestIdx[nearestIdx].push_back(affectedIdx);
        }
    }

    rotatePoints(points, -rot);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Merge done in " << duration.count() << " ms" << std::endl;

    // Store points in a circular linked list
    std::list<UnmergeNode> tour;
    std::vector<std::list<UnmergeNode>::iterator> tourNodes(points.size()); 

    // Create the starting node and add it to the tour
    tourNodes[points.size() - 1] = tour.insert(tour.end(), UnmergeNode(points.back(), points.size() - 1));

    // Decoding phase
    auto nodeComparator = [&nodes](size_t lhs, size_t rhs) {
        return nodes[lhs]->distance < nodes[rhs]->distance; // Max-heap based on distance
    };

    std::priority_queue<size_t, std::vector<size_t>, decltype(nodeComparator)> maxHeap(nodeComparator);
    maxHeap.push(points.size() - 1);

    RTree pointTree; 

    SegmentRTree segTree; 

    struct PointComparator {
        bool operator()(const Point& lhs, const Point& rhs) const {
            if (bg::get<0>(lhs) != bg::get<0>(rhs))
                return bg::get<0>(lhs) < bg::get<0>(rhs);
            return bg::get<1>(lhs) < bg::get<1>(rhs);
        }
    };
    std::map<Point, std::list<UnmergeNode>::iterator, PointComparator> pointMap;

    auto removeNode = [&](std::list<UnmergeNode>::iterator itr) -> void {
        // std::cout << "removing node " << bg::get<0>(itr->point) << ' ' << bg::get<1>(itr->point) << std::endl;
        auto prevNode = *prevIt(itr, tour);
        auto nextNode = *nextIt(itr, tour);
        pointTree.remove(itr->point);
        pointMap.erase(itr->point);
        segTree.remove(Segment(itr->point, prevNode.point));
        if (!bg::equals(prevNode.point, nextNode.point)) {
            segTree.remove(Segment(itr->point, nextNode.point));
            segTree.insert(Segment(prevNode.point, nextNode.point));
        }
        tour.erase(itr);
    };

    auto insertNodeUtil = [&](size_t idx, bool depleteCapacity, const auto& self) -> void {
        // std::cout << "inserting " << idx << " " << depleteCapacity << std::endl;

        Point point = points[idx];
        double radius = distances[idx];

        // If closest point is within radius, we simply add to this point
        Point closestPoint;
        pointTree.query(bgi::nearest(point, 1), &closestPoint);
        if (bg::distance(closestPoint, point) <= radius) {
            assert(pointMap.find(closestPoint) != pointMap.end());
            auto it = pointMap[closestPoint];
            it->indexes.insert(idx);
            it->capacity++;
            tourNodes[idx] = it;
            return;
        }

        // If there is only one point in the tree: just add a new one
        if (tour.size() == 1) {
            auto p = tour.begin()->point;
            auto optPoint = findOptimalPoint(point, radius, p, p);
            auto itr = tour.insert(tour.end(), UnmergeNode(optPoint, idx));
            tourNodes[idx] = itr;
            segTree.insert(Segment(p, optPoint));
            pointTree.insert(optPoint);
            pointMap[optPoint] = itr;
            return;
        }

        // Else, query the closest segment and add a new point in the middle
        Segment seg = findBestSegment(segTree, point);
        auto [p1, p2] = seg;
        Point optPoint = findOptimalPoint(point, radius, p1, p2);
        if (pointMap.find(p1) == pointMap.end() || pointMap.find(p2) == pointMap.end()) {
            std::cout << "nononono" << std::endl;
            std::cout << bg::get<0>(p1) << ' ' << bg::get<1>(p1) << std::endl;
            std::cout << bg::get<0>(p2) << ' ' << bg::get<1>(p2) << std::endl;
            std::cout << "MAP: " << std::endl;
            for (const auto& [key, _]: pointMap) {
                std::cout << bg::get<0>(key) << ' ' << bg::get<1>(key) << std::endl;
            }
        }
        assert(pointMap.find(p1) != pointMap.end());
        assert(pointMap.find(p2) != pointMap.end());
        auto itrPrev = pointMap[p1], itrNext = pointMap[p2];
        if (nextIt(itrPrev, tour) != itrNext) {
            std::swap(itrPrev, itrNext);
            assert(nextIt(itrPrev, tour) == itrNext);
        }
        auto itr = tour.insert(itrNext, UnmergeNode(optPoint, idx));
        tourNodes[idx] = itr;
        segTree.remove(seg);
        segTree.insert(Segment(p1, optPoint));
        segTree.insert(Segment(p2, optPoint));
        pointTree.insert(optPoint);
        pointMap[optPoint] = itr;

        // Change the neighboring nodes
        if (depleteCapacity) {
            itrPrev->capacity--;
            if (itrPrev->capacity == 0) {
                std::vector<size_t> reinsertOrder(itrPrev->indexes.begin(), itrPrev->indexes.end());
                std::sort(reinsertOrder.begin(), reinsertOrder.end(), [&](size_t i, size_t j) {return distances[i] < distances[j];});
                removeNode(itrPrev);
                for (size_t idx: reinsertOrder) {
                    self(idx, false, self);
                }
            }

            itrNext->capacity--;
            if (itrNext->capacity == 0) {
                std::vector<size_t> reinsertOrder(itrNext->indexes.begin(), itrNext->indexes.end());
                std::sort(reinsertOrder.begin(), reinsertOrder.end(), [&](size_t i, size_t j) {return distances[i] < distances[j];});
                removeNode(itrNext);
                for (size_t idx: reinsertOrder) {
                    self(idx, false, self);
                }
            }
        }
    };

    auto insertNode = [&](size_t idx, bool depleteCapacity = true) {
        insertNodeUtil(idx, depleteCapacity, insertNodeUtil);
    };

    auto optimizeTour = [&]() {
        // std::cout << "opt " << tour.size() << std::endl;

        for (std::list<UnmergeNode>::iterator itr = tour.begin(), _nextItr; itr != tour.end(); itr = _nextItr) {
            _nextItr = std::next(itr);

            if (itr->indexes.size() > 1) {
                continue;
            }

            auto prevItr = prevIt(itr, tour);
            auto nextItr = nextIt(itr, tour);
            
            size_t idx = *(itr->indexes.begin());
            Point optPoint = findOptimalPoint(points[idx], distances[idx], prevItr->point, nextItr->point);

            static constexpr double EQ_TOLERANCE = 1e-12;
            if (bg::distance(optPoint, itr->point) > EQ_TOLERANCE) {
                if (pointMap.find(optPoint) != pointMap.end()) {
                    removeNode(itr);
                    tourNodes[idx] = pointMap[optPoint];
                    tourNodes[idx]->indexes.insert(idx);
                    tourNodes[idx]->capacity++;
                } 
                else {
                    segTree.remove(Segment(prevItr->point, itr->point));
                    segTree.remove(Segment(nextItr->point, itr->point));
                    segTree.insert(Segment(prevItr->point, optPoint));
                    segTree.insert(Segment(nextItr->point, optPoint));
                    pointTree.remove(itr->point);
                    pointTree.insert(optPoint);
                    pointMap.erase(itr->point);
                    pointMap[optPoint] = itr;
                    itr->point = optPoint;
                }
            }
        }
    };

    size_t turnCount = 0;
    auto isOptValue = [&](size_t n) -> bool {
        return n > 2 && (n & (n - 1)) == 0;
    };

    while (!maxHeap.empty()) {
        turnCount++;

        size_t maxIdx = maxHeap.top();
        maxHeap.pop();

        MergeNode* maxNode = nodes[maxIdx];

        if (!maxNode->isLeaf()) {
            // Retrieve indices of the left and right children
            size_t leftIdx = maxNode->a->idx, rightIdx = maxNode->b->idx;
            bool leftIsLeaf = maxNode->a->isLeaf(), rightIsLeaf = maxNode->b->isLeaf();

            // std::cout << "Processing " << maxIdx << ' ' << leftIdx << ' ' << rightIdx << ' ' << tour.size() << std::endl;

            // Handle the case with a single node separately  
            if (tour.size() == 1) {
                Point l = findOptimalPoint(points[leftIdx], distances[leftIdx], points[rightIdx], points[rightIdx]);
                Point r = findOptimalPoint(points[rightIdx], distances[rightIdx], l, l);
                pointTree.insert(l);
                pointTree.insert(r);
                segTree.insert(Segment(l, r));
                tour.clear();
                auto leftItr = tour.insert(tour.end(), UnmergeNode(l, leftIdx));
                auto rightItr = tour.insert(tour.end(), UnmergeNode(r, rightIdx));
                tourNodes[leftIdx] = leftItr;
                tourNodes[rightIdx] = rightItr;
                pointMap[l] = leftItr;
                pointMap[r] = rightItr;
            } else {
                // Remove the current node
                auto curNodeIt = tourNodes[maxIdx];
                auto& curNode = *curNodeIt;
                if (curNode.indexes.size() > 1) { // If this point is shared with other base points
                    curNode.indexes.erase(maxIdx); // Simply erase this index from the indexes
                } else {
                    removeNode(curNodeIt);
                }

                insertNode(leftIdx);
                insertNode(rightIdx);
            }

            // Push children to the heap if they are not leaves
            if (!leftIsLeaf) {
                maxHeap.push(leftIdx);
            }
            if (!rightIsLeaf) {
                maxHeap.push(rightIdx);
            }
        }

        // Optimize the provisional TSP every so often
        if (isOptValue(turnCount) && tour.size() > 2) {
            optimizeTour();
        }

        /*std::cout << "pointMap" << std::endl;
        for (auto& [point, itr]: pointMap) {
            std::cout << bg::get<0>(point) << " " << bg::get<1>(point) << std::endl;
        }

        std::cout << "pointTree" << std::endl;
        for (auto& point: pointTree) {
            std::cout << bg::get<0>(point) << " " << bg::get<1>(point) << std::endl;
        }*/

        if (pointMap.size() != pointTree.size()) {
            std::cout << "MAP: " << std::endl;
            for (const auto& [key, _]: pointMap) {
                std::cout << bg::get<0>(key) << ' ' << bg::get<1>(key) << std::endl;
            }
            std::cout << "TREE: " << std::endl;
            for (const auto& key: pointTree) {
                std::cout << bg::get<0>(key) << ' ' << bg::get<1>(key) << std::endl;
            }
        }

        assert(pointMap.size() == pointTree.size());
    }

    // Final opt
    optimizeTour();

    // Convert to vector format
    std::vector<TSPNode> tsp;
    for (const auto& node : tour) {
        for (const auto& idx : node.indexes) {
            tsp.emplace_back(points[idx], distances[idx], node.point);
        }
    }

    /*std::vector<Point> p;
    for (const auto& node: tsp) {
        p.push_back(node.tourPoint);
    }

    auto reorder = approximateTSP(p);
    std::vector<TSPNode> newTSP;
    for (const auto& idx : reorder) {
        newTSP.push_back(tsp[idx]);
    }
    tsp = std::move(newTSP);*/

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "All done in " << duration.count() << " ms" << std::endl;

    return tsp;
}

std::vector<TSPNode> approximateCETSP(std::vector<Point> points, std::vector<double> distances) {
    static constexpr size_t RETRIES = 10;
    double best = std::numeric_limits<double>::max();
    std::vector<TSPNode> bestTSP;
    for (size_t i = 0; i < RETRIES; i++) {
        std::vector<TSPNode> tsp = makePairCenterCETSP(points, distances);
        auto dist = calculateTotalTourDistance(tsp);
        if (dist < best) {
            best = dist; 
            bestTSP = tsp;
        }
    }
    return bestTSP;
}

// Function to read data from a file
bool readData(const fs::path& filename, std::vector<Point>& points, std::vector<double>& distances) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    points.clear();
    distances.clear();

    double x, y, distance;
    while (infile >> x >> y >> distance) {
        points.emplace_back(x, y);
        distances.push_back(distance);
    }

    infile.close();
    return true;
}

// Function to write results to a file
void writeResults(const fs::path& filename, const std::vector<TSPNode>& tour) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error writing to file: " << filename << std::endl;
        return;
    }

    outfile << "Decoded TSP tour:\n";

    // Iterate through the TSPNodeList
    for (const auto& node: tour) {
        const auto& b = node.basePoint;
        const auto& p = node.tourPoint;
        outfile << bg::get<0>(b) << " " << bg::get<1>(b) << " " << node.distance << " " << bg::get<0>(p) << " " << bg::get<1>(p) << "\n";
    }

    // Write the total distance
    outfile << "Total tour distance: " << calculateTotalTourDistance(tour) << std::endl;

    outfile.close();
}

int main() {
    const fs::path inputDir("./data/cmu_processed");
    const fs::path outputDir("./data/cmu_output");

    std::cout << "BEGIN!" << std::endl;

    double totalDist = 0.0;

    try {
        // Ensure the output directory exists
        if (!fs::exists(outputDir)) {
            fs::create_directories(outputDir);
        }

        // Iterate through all files in the input directory
        for (const auto& entry : fs::directory_iterator(inputDir)) {
            if (!entry.is_regular_file()) continue;

            fs::path inputFile = entry.path();
            fs::path outputFile = outputDir / entry.path().filename();

            std::cout << "Processing: " << inputFile << " -> " << outputFile << std::endl;

            std::vector<Point> points;
            std::vector<double> distances;

            // Read data from the input file
            if (!readData(inputFile, points, distances)) {
                std::cerr << "Failed to read data from: " << inputFile << std::endl;
                continue;
            }

            // Run approximate CETSP
            auto tspResult = approximateCETSP(points, distances);

            for (const auto& node: tspResult) {
                if (bg::distance(node.basePoint, node.tourPoint) > node.distance) {
                    std::cout << "WTF\n";
                    const auto& b = node.basePoint;
                    const auto& p = node.tourPoint;
                    std::cout << bg::get<0>(b) << " " << bg::get<1>(b) << " " << node.distance << " " << bg::get<0>(p) << " " << bg::get<1>(p) << std::endl;
                    assert(0 && "Algorithm pushed point outside range, exiting.");
                }
            }

            totalDist += calculateTotalTourDistance(tspResult);

            // Write results to the output file
            writeResults(outputFile, tspResult);

            std::cout << "Processed: " << inputFile << " -> " << outputFile << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    std::cout << "Processing complete!" << std::endl;
    std::cout << "Total distance: " << totalDist << std::endl;

}
