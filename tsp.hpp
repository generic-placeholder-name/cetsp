#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <list>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using Point = bg::model::point<double, 2, bg::cs::cartesian>;
// Use R*-tree instead of the standard quadratic R-tree
using RTree = bgi::rtree<std::pair<Point, size_t>, bgi::rstar<16>>;

// Struct to hold the distance between two points for the priority queue
struct DistanceEntry {
    double distance;
    size_t idxA, idxB;

    bool operator<(const DistanceEntry& other) const {
        return distance > other.distance;
    }
};

struct Node {
    double distance;
    Point point;
    Node* a;
    Node* b;

    Node(double dist, const Point& pt, Node* childA = nullptr, Node* childB = nullptr)
        : distance(dist), point(pt), a(childA), b(childB) {}
    
    bool isLeaf() const {
        return !a && !b;
    }
};

// Function to calculate the midpoint between two points
Point midpointBetween(const Point& a, const Point& b) {
    return Point((bg::get<0>(a) + bg::get<0>(b)) / 2.0,
                 (bg::get<1>(a) + bg::get<1>(b)) / 2.0);
}

// Main TSP approximation function
void approximate_tsp(std::vector<Point> points) {
    RTree tree;
    for (size_t i = 0; i < points.size(); i++) {
        tree.insert({points[i], i});
    }

    std::priority_queue<DistanceEntry> distanceQueue;
    std::vector<Node*> nodes(points.size(), nullptr);
    std::vector<bool> pointExists(points.size(), true);

    // Initialize distances and priority queue
    std::vector<std::vector<size_t>> closestIdx(points.size()); // Points that have point i as their closest point

    for (size_t i = 0; i < points.size(); i++) {
        std::vector<std::pair<Point, size_t>> result;
        tree.query(bgi::nearest(points[i], 2), std::back_inserter(result));

        if (result.size() >= 2) {
            auto nearestPoint = result[1];
            double dist = bg::distance(points[i], nearestPoint.first);
            size_t nearestIdx = nearestPoint.second;
            distanceQueue.push({dist, i, nearestIdx});
            closestIdx[nearestIdx].push_back(i);
        }

        // Create leaf nodes
        nodes[i] = new Node(0.0, points[i]);
    }

    // Iteratively pop pairs of closest points and replace with midpoints
    while (!distanceQueue.empty()) {
        auto [d, idxA, idxB] = distanceQueue.top();
        distanceQueue.pop();

        if (!pointExists[idxA] || !pointExists[idxB]) continue;

        auto& a = points[idxA], b = points[idxB];
        Point midpoint = midpointBetween(a, b);
        size_t newIdx = points.size();
        points.push_back(midpoint);
        pointExists.push_back(true);

        tree.insert({midpoint, newIdx});
        tree.remove({a, idxA});
        tree.remove({b, idxB});
        pointExists[idxA] = pointExists[idxB] = false;

        // Create a new node for the midpoint
        nodes.push_back(new Node(d, midpoint, nodes[idxA], nodes[idxB]));

        // Update nearest neighbors and distances for affected points
        // Points that had idxA or idxB as their closest neighbor need updates
        std::vector<size_t> affectedPoints = closestIdx[idxA];
        affectedPoints.insert(affectedPoints.end(), closestIdx[idxB].begin(), closestIdx[idxB].end());

        closestIdx[idxA].clear();
        closestIdx[idxB].clear();

        for (size_t affectedIdx : affectedPoints) {
            if (!pointExists[affectedIdx]) continue;

            std::vector<std::pair<Point, size_t>> result;
            tree.query(bgi::nearest(points[affectedIdx], 2), std::back_inserter(result));

            if (result.size() >= 2) {
                auto nearestPoint = result[1];
                double dist = bg::distance(points[affectedIdx], nearestPoint.first);
                size_t nearestIdx = nearestPoint.second;
                distanceQueue.push({dist, affectedIdx, nearestIdx});
                closestIdx[nearestIdx].push_back(affectedIdx);
            }
        }
    }

    // Store points in a circular linked list
    std::list<Node*> tour;
    for (size_t i = 0; i < points.size(); i++) {
        if (pointExists[i]) {
            tour.push_back(nodes[i]);
        }
    }

    assert(tour.size() == 1);

    // Decoding phase
    std::priority_queue<Node*> maxHeap;
    maxHeap.push(tour.back());

    while (!maxHeap.empty()) {
        Node* maxNode = maxHeap.top();
        maxHeap.pop();

        auto it = std::find(tour.begin(), tour.end(), maxNode);
        if (it != tour.end()) {
            tour.erase(it);
            tour.push_back(maxNode->a);
            tour.push_back(maxNode->b);

            if (!maxNode->a->isLeaf()) {
                maxHeap.push(maxNode->a);
            }
            if (!maxNode->b->isLeaf()) {
                maxHeap.push(maxNode->b);
            }
        }
    }

    // Print the tour
    std::cout << "Decoded TSP tour:" << std::endl;
    for (Node* node : tour) {
        const auto& p = node->point;
        std::cout << "(" << bg::get<0>(p) << ", " << bg::get<1>(p) << ")" << std::endl;
    }
}

int main() {
    // Example points
    std::vector<Point> points = {
        Point(0, 0), Point(1, 1), Point(2, 2), Point(3, 3), Point(4, 4)
    };

    approximate_tsp(points);

    std::cout << "TSP approximation complete." << std::endl;
    return 0;
}
