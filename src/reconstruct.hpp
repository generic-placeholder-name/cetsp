#pragma once

#include <cetsp/types.hpp>

#include "merge_tree.hpp"

#include <vector>

// ==================== Tour reconstruction ====================

// Given a merge tree (vector<TreeNode>), reconstruct the tour by "unmerging" nodes.
// At each step, remove the highest-gap internal node from the tour and insert its children.
// Returns the sequence of tour points.
std::vector<Point> reconstructTour(const std::vector<TreeNode>& treeNodes);
