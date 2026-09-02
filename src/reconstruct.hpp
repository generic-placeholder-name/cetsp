#pragma once

#include <cetsp/types.hpp>

#include "merge_tree.hpp"

#include <vector>

// ==================== Tour reconstruction ====================

// Reconstruct a tour by unmerging the highest-gap branches and inserting their
// children.
std::vector<Point> reconstructTour(const MergeTree& mergeTree);
