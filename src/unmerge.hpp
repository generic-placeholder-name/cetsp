#pragma once

#include <vector>
#include "structures.hpp"
#include "geometry.hpp"

// ==================== Tour reconstruction ====================

// Given a merge tree (vector<TreeNode>), reconstruct the tour by "unmerging" nodes.
// At each step, remove the highest-weight internal node from the tour and insert its children.
// Returns the sequence of tour points.
std::vector<Point> unmerge(std::vector<TreeNode> treeNodes);
