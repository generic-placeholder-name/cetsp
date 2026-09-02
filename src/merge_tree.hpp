#pragma once

#include <cetsp/types.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace cetsp_detail {
class MergeTreeBuilder;
}

// Immutable full-binary merge history. Every non-root node has exactly one
// parent, every branch has two children, and node IDs remain stable for the
// lifetime of the tree.
class MergeTree {
public:
    using NodeId = std::size_t;

    MergeTree() = default;

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::size_t leafCount() const noexcept;
    [[nodiscard]] std::optional<NodeId> root() const noexcept;

    [[nodiscard]] const Circle& neighborhood(NodeId node) const;
    [[nodiscard]] bool isLeaf(NodeId node) const;
    [[nodiscard]] std::array<NodeId, 2> children(NodeId node) const;
    [[nodiscard]] double mergeGap(NodeId node) const;

private:
    struct Branch {
        std::array<NodeId, 2> children;
        double mergeGap;
    };

    struct Node {
        Circle neighborhood;
        std::optional<Branch> branch;
    };

    explicit MergeTree(std::vector<Node> nodes);

    [[nodiscard]] const Node& node(NodeId id) const;
    void validate() const;

    std::vector<Node> nodes_;

    friend class cetsp_detail::MergeTreeBuilder;
};
