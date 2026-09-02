#pragma once

#include <cetsp/types.hpp>

#include "merge_tree.hpp"

#include <random>
#include <vector>

// ==================== Merge phase ====================

// Remove redundant circles that completely contain another circle. Any point
// visiting the contained circle necessarily visits the containing circle too.
// Mutates the vector but does not attach solver state to the remaining inputs.
void removeCoveringCircles(std::vector<Circle>& circles);

// Build an immutable merge history from the given circles.
MergeTree buildMergeTree(
    const std::vector<Circle>& circles,
    std::mt19937_64& randomEngine);

