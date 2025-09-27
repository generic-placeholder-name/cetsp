#pragma once

#include <vector>
#include "structures.hpp"
#include "geometry.hpp"

// ==================== Merge phase ====================

// Remove circles that are completely covered by others.
// Mutates the vector of circles and reassigns nodeIds.
void removeCoveringCircles(std::vector<Circle>& circles);

// Build a merge tree from the given circles.
// Returns a vector of TreeNodes representing the merge process.
std::vector<TreeNode> buildMergeTree(std::vector<Circle> circles);

