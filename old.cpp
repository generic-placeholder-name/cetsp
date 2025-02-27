#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <vector>
#include <queue>
#include <cmath>
#include <list>
#include <random>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/intrusive/list.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bi = boost::intrusive;
namespace fs = std::filesystem;

using Point = bg::model::point<double, 2, bg::cs::cartesian>;
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
    size_t idx;  // Index of the point in the `points` vector
    Node* a;
    Node* b;

    Node(double dist, size_t index, Node* childA = nullptr, Node* childB = nullptr)
        : distance(dist), idx(index), a(childA), b(childB) {}

    ~Node() {
        delete a;
        delete b;
    }

    bool isLeaf() const {
        return !a && !b;
    }
};

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

    // Hook for Boost intrusive list
    // bi::list_member_hook<> memberHook;

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

struct IndexNode : public bi::list_base_hook<> {
    size_t index;
    IndexNode(const size_t& index) : index(index) {}
};

// using TSPNodeList = bi::list<TSPNode, bi::member_hook<TSPNode, bi::list_member_hook<>, &TSPNode::memberHook>>;
using IndexList = bi::list<IndexNode>;

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

// TSP approximation
std::vector<size_t> approximateTSP(std::vector<Point> points) {
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
        nodes[i] = new Node(0.0, i);
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
        nodes.push_back(new Node(d, newIdx, nodes[idxA], nodes[idxB]));

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
    IndexList tour;

    // Create the starting node and add it to the tour
    IndexNode* tourStart = new IndexNode(points.size() - 1);
    tour.push_back(*tourStart);

    // Decoding phase
    auto nodeComparator = [&nodes](IndexNode* lhs, IndexNode* rhs) {
        return nodes[lhs->index]->distance < nodes[rhs->index]->distance; // Max-heap based on distance
    };

    std::priority_queue<IndexNode*, std::vector<IndexNode*>, decltype(nodeComparator)> maxHeap(nodeComparator);
    maxHeap.push(tourStart);

    while (!maxHeap.empty()) {
        IndexNode* curNode = maxHeap.top();
        maxHeap.pop();

        size_t maxIdx = curNode->index;
        Node* maxNode = nodes[maxIdx];

        if (!maxNode->isLeaf()) {
            // Retrieve indices of the left and right children
            size_t leftIdx = maxNode->a->idx, rightIdx = maxNode->b->idx;
            bool leftIsLeaf = maxNode->a->isLeaf(), rightIsLeaf = maxNode->b->isLeaf();

            // Find neighboring nodes in the tour
            auto curNodeIt = IndexList::s_iterator_to(*curNode);
            auto prevNodeIt = prevIt(curNodeIt, tour);
            auto nextNodeIt = nextIt(curNodeIt, tour);

            IndexNode& prevNode = *prevNodeIt;
            IndexNode& nextNode = *nextNodeIt;

            // Calculate distances for insertion orders
            double distPrevLeftRightNext =
                bg::distance(points[prevNode.index], points[leftIdx]) +
                bg::distance(points[rightIdx], points[nextNode.index]);

            double distPrevRightLeftNext =
                bg::distance(points[prevNode.index], points[rightIdx]) +
                bg::distance(points[leftIdx], points[nextNode.index]);

            // Initialize new nodes corresponding to the children of the current node
            IndexNode* leftNode = new IndexNode(leftIdx);
            IndexNode* rightNode = new IndexNode(rightIdx);

            // Insert nodes based on the smaller distance
            if (distPrevLeftRightNext < distPrevRightLeftNext) {
                tour.insert(std::next(curNodeIt), *leftNode);
                tour.insert(std::next(IndexList::s_iterator_to(*leftNode)), *rightNode);
            } else {
                tour.insert(std::next(curNodeIt), *rightNode);
                tour.insert(std::next(IndexList::s_iterator_to(*rightNode)), *leftNode);
            }

            // Erase the current node
            tour.erase(curNodeIt);
            delete curNode;

            // Push children to the heap if they are not leaves
            if (!leftIsLeaf) {
                maxHeap.push(leftNode);
            }
            if (!rightIsLeaf) {
                maxHeap.push(rightNode);
            }
        }
    }

    // Convert back to vector for convenience
    std::vector<size_t> ans;
    ans.reserve(points.size());
    for (const auto& node: tour) {
        ans.push_back(node.index);
    }

    // Clean up
    for (auto it = tour.begin(); it != tour.end();) {
        auto& node = *it;
        ++it;  
        tour.erase(tour.iterator_to(node));
        delete &node;
    }

    return ans;
}

std::vector<TSPNode> approximateCETSP(std::vector<Point> points, std::vector<double> distances) {
    std::vector<size_t> indexes = approximateTSP(points);

    std::vector<TSPNode> tsp;
    for (const auto& index: indexes) {
        tsp.emplace_back(points[index], distances[index]);
    }

    auto adjust = [&tsp](size_t outerIters, size_t innerIters) -> void {
        for (int it = 0; it < outerIters; it++) {
            for (size_t i = 0; i < tsp.size(); i++) {
                if (tsp[i].distance == 0) {
                    continue;
                }

                size_t prevIdx = (i == 0 ? tsp.size() - 1 : i - 1);
                size_t nextIdx = (i == tsp.size() - 1 ? 0 : i + 1);
                double extraDist = bg::distance(tsp[prevIdx].tourPoint, tsp[i].tourPoint) 
                                 + bg::distance(tsp[i].tourPoint, tsp[nextIdx].tourPoint)
                                 - bg::distance(tsp[prevIdx].tourPoint, tsp[nextIdx].tourPoint);
                double rat = std::min(extraDist / tsp[i].distance, 1.0);
                tsp[i].perturb(rat / (it + 8));
            }
            for (int iters = 0; iters < innerIters; iters++) {
                for (size_t i = 0; i < tsp.size(); i++) {
                    size_t prevIdx = (i == 0 ? tsp.size() - 1 : i - 1);
                    size_t nextIdx = (i == tsp.size() - 1 ? 0 : i + 1);
                    tsp[i].tourPoint = findOptimalPoint(tsp[i].basePoint, tsp[i].distance, tsp[prevIdx].tourPoint, tsp[nextIdx].tourPoint);
                }
            }
        }
    };

    auto trim = [&tsp]() -> void {
        const size_t n = tsp.size();

        std::vector<size_t> anchor(n);
        std::iota(anchor.begin(), anchor.end(), 0);

        // Indexes of nodes
        std::vector<size_t> indexes(n);
        std::iota(indexes.begin(), indexes.end(), 0);

        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::shuffle(indexes.begin(), indexes.end(), gen);

        RTree tree;
        for (size_t index: indexes) {
            if (tree.empty()) {
                tree.insert({tsp[index].tourPoint, index});
                continue;
            }

            std::pair<Point, size_t> point;
            tree.query(bgi::nearest(tsp[index].tourPoint, 1), &point);
            double dist = bg::distance(tsp[index].basePoint, point.first);
            if (dist > tsp[index].distance) {
                tree.insert({tsp[index].tourPoint, index});
            } else {
                anchor[index] = point.second;
            }
        }

        std::vector<std::vector<size_t>> nodesInIndex(n);
        for (size_t i : indexes) {
            nodesInIndex[anchor[i]].push_back(i);
        }

        std::vector<Point> anchorPoints;
        std::vector<size_t> anchorIndexes;

        for(size_t i = 0; i < n; i++) {
            if (!nodesInIndex[i].empty()) {
                anchorPoints.push_back(tsp[i].tourPoint);
                anchorIndexes.push_back(i);
            }
        }

        std::vector<size_t> reorderedIndexes = approximateTSP(anchorPoints);

        std::vector<TSPNode> newTSP;
        for (size_t x: reorderedIndexes) {
            size_t i = anchorIndexes[x];
            for (const auto& idx: nodesInIndex[i]) {
                newTSP.push_back(tsp[idx]);
                newTSP.back().tourPoint = tsp[i].tourPoint;
            }
        }

        for (size_t i = 0; i < n; i++) {
            size_t prevIdx = (i == 0 ? n - 1 : i - 1);
            size_t nextIdx = (i == n - 1 ? 0 : i + 1);
            newTSP[i].tourPoint = findOptimalPoint(newTSP[i].basePoint, newTSP[i].distance, newTSP[prevIdx].tourPoint, newTSP[nextIdx].tourPoint);
        }

        tsp = std::move(newTSP);
    };

    auto rearrangeTSP = [&tsp]() -> void {
        std::vector<Point> points;
        for (auto& node : tsp) {
            points.push_back(node.tourPoint);
        }

        auto indexes = approximateTSP(points);

        std::vector<TSPNode> newTSP;

        // Rebuild the list based on the new order 
        for (const auto& index : indexes) {
            newTSP.push_back(tsp[index]);
        }

        tsp = std::move(newTSP); 
    };

    std::cout << "Current TSP dist: " << calculateTotalTourDistance(tsp) << std::endl;
    adjust(2, 2);
    std::cout << "Current TSP dist: " << calculateTotalTourDistance(tsp) << std::endl;
    trim();
    std::cout << "Current TSP dist: " << calculateTotalTourDistance(tsp) << std::endl;
    adjust(3, 8);
    std::cout << "Current TSP dist: " << calculateTotalTourDistance(tsp) << std::endl;
    trim();
    std::cout << "Current TSP dist: " << calculateTotalTourDistance(tsp) << std::endl;
    // adjust(3, 16);
    // std::cout << "Current TSP dist: " << calculateTotalTourDistance(tsp) << std::endl;

    return tsp;
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
