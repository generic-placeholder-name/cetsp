#pragma once

#include <cetsp/types.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <limits>

using TreeNodeId = std::size_t;

// Internal representation of the binary merge history. Nodes are stored in a
// vector and children always precede their parent.
struct TreeNode {
    TreeNodeId left;
    TreeNodeId right;
    double weight;
    Point center;
    double r;

    [[nodiscard]] static TreeNode leaf(const Point& center, double radius) {
        return TreeNode(noChild, noChild, 0.0, center, radius);
    }

    [[nodiscard]] static TreeNode branch(
        TreeNodeId left,
        TreeNodeId right,
        double weight,
        const Point& center,
        double radius) {
        assert(left != noChild && right != noChild);
        return TreeNode(left, right, weight, center, radius);
    }

    [[nodiscard]] bool isLeaf() const noexcept {
        assert((left == noChild) == (right == noChild) &&
               "a merge-tree node must have either zero or two children");
        return left == noChild;
    }

    [[nodiscard]] std::array<TreeNodeId, 2> children() const noexcept {
        assert(!isLeaf() && "a merge-tree leaf has no children");
        return {left, right};
    }

private:
    static constexpr TreeNodeId noChild =
        std::numeric_limits<TreeNodeId>::max();

    TreeNode(
        TreeNodeId left,
        TreeNodeId right,
        double weight,
        const Point& center,
        double radius)
        : left(left), right(right), weight(weight), center(center), r(radius) {}
};
