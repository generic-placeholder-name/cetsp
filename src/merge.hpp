#pragma once

#include <cetsp/types.hpp>

#include "merge_tree.hpp"

#include <random>
#include <vector>

// ==================== Merge phase ====================

// Remove circles that are completely covered by others.
// Mutates the vector but does not attach solver state to the remaining inputs.
void removeCoveringCircles(std::vector<Circle>& circles);

// Build a merge tree from the given circles.
// Returns a vector of TreeNodes representing the merge process.
std::vector<TreeNode> buildMergeTree(
    const std::vector<Circle>& circles,
    std::mt19937_64& randomEngine);

