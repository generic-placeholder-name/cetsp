#pragma once

#include <cstddef>
#include <unordered_set>
#include <absl/container/flat_hash_set.h>
#include <boost/geometry.hpp>

// Choose Abseil or std hash set
#if USE_ABSEIL_HASH_SET
template<typename T>
using HashSet = absl::flat_hash_set<T>;
#else
template<typename T>
using HashSet = std::unordered_set<T>;
#endif

namespace bg = boost::geometry;

// Forward declare Point (so we don’t need full geometry.hpp here)
typedef bg::model::point<double, 2, bg::cs::cartesian> Point;

// ==================== Structures ====================

/*
* TreeNode represents a node in the merge tree.
* Since merge trees are stored as std::vectors of TreeNode,
* we use size_t indices to refer to child nodes.
*/
struct TreeNode { 
    size_t left, right; 
    double weight; 
    Point center;
    double r;

    TreeNode(size_t l, size_t r, double w, Point c = Point(), double ra = 0.0) 
        : left(l), right(r), weight(w), center(c), r(ra) {}
};

// Circle representation
struct Circle {
    Point center; // center of the circle
    double r; // radius
    size_t nodeId; // index of the corresponding TreeNode in the merge tree
    size_t nn; // nearest neighbor circle index
    double gap; // gap distance to nearest neighbor
    HashSet<size_t> rev; // reverse neighbors (circles that have this as nn)
};

// Node in the cyclic tour
struct TourNode {
    Point pos;
    HashSet<size_t> circles;
    size_t energy; // we pop and reinsert when this is 2^n
    size_t insertions; // we optimize point location when this is 2^n
    TourNode* prev;
    TourNode* next;
    size_t id; // circle index

    TourNode(const Point& p, size_t id_ = static_cast<size_t>(-1))
        : pos(p), energy(0), insertions(0), prev(nullptr), next(nullptr), id(id_) {
        circles.insert(id_);
    }
};