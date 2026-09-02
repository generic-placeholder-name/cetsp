#include "merge_tree.hpp"

#include <boost/geometry/core/access.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace bg = boost::geometry;

MergeTree::MergeTree(std::vector<Node> nodes)
    : nodes_(std::move(nodes)) {
    validate();
}

bool MergeTree::empty() const noexcept {
    return nodes_.empty();
}

std::size_t MergeTree::size() const noexcept {
    return nodes_.size();
}

std::size_t MergeTree::leafCount() const noexcept {
    return nodes_.empty() ? 0 : (nodes_.size() + 1) / 2;
}

std::optional<MergeTree::NodeId> MergeTree::root() const noexcept {
    if (nodes_.empty()) {
        return std::nullopt;
    }
    return nodes_.size() - 1;
}

const Circle& MergeTree::neighborhood(NodeId id) const {
    return node(id).neighborhood;
}

bool MergeTree::isLeaf(NodeId id) const {
    return !node(id).branch.has_value();
}

std::array<MergeTree::NodeId, 2> MergeTree::children(NodeId id) const {
    const Node& value = node(id);
    if (!value.branch) {
        throw std::logic_error("a merge-tree leaf has no children");
    }
    return value.branch->children;
}

double MergeTree::mergeGap(NodeId id) const {
    const Node& value = node(id);
    if (!value.branch) {
        throw std::logic_error("a merge-tree leaf has no merge gap");
    }
    return value.branch->mergeGap;
}

const MergeTree::Node& MergeTree::node(NodeId id) const {
    if (id >= nodes_.size()) {
        throw std::out_of_range("merge-tree node ID is out of range");
    }
    return nodes_[id];
}

void MergeTree::validate() const {
    if (nodes_.empty()) {
        return;
    }
    if (nodes_.size() % 2 == 0) {
        throw std::invalid_argument(
            "a nonempty full binary merge tree must have an odd node count");
    }

    const std::size_t leaves = leafCount();
    std::vector<std::size_t> parentCounts(nodes_.size(), 0);
    for (NodeId id = 0; id < nodes_.size(); ++id) {
        const Node& value = nodes_[id];
        const Circle& circle = value.neighborhood;
        if (!std::isfinite(bg::get<0>(circle.center)) ||
            !std::isfinite(bg::get<1>(circle.center)) ||
            !std::isfinite(circle.r) || circle.r < 0.0) {
            throw std::invalid_argument(
                "merge-tree neighborhoods must be finite with non-negative radii");
        }

        const bool expectedLeaf = id < leaves;
        if (expectedLeaf != !value.branch.has_value()) {
            throw std::invalid_argument(
                "merge-tree leaves must precede all branch nodes");
        }
        if (!value.branch) {
            continue;
        }
        if (!std::isfinite(value.branch->mergeGap)) {
            throw std::invalid_argument("merge-tree gaps must be finite");
        }

        const auto [left, right] = value.branch->children;
        if (left == right || left >= id || right >= id) {
            throw std::invalid_argument(
                "merge-tree branch children must be distinct earlier nodes");
        }
        ++parentCounts[left];
        ++parentCounts[right];
    }

    const NodeId rootId = nodes_.size() - 1;
    for (NodeId id = 0; id < nodes_.size(); ++id) {
        const std::size_t expectedParents = id == rootId ? 0 : 1;
        if (parentCounts[id] != expectedParents) {
            throw std::invalid_argument(
                "every non-root merge-tree node must have exactly one parent");
        }
    }
}
