#include "CETSP.hpp"
#include "merge.hpp"
#include "unmerge.hpp"
#include "debug.hpp"

#include <vector>
#include <limits>
#include <iostream>
#include <iterator>

// Verify that the tour is valid: each circle should contain a point in the tour.
bool verifyTour(const std::vector<Point>& tour, const std::vector<Circle>& circles) {
    // Construct an R*-tree from the tour points
    bgi::rtree<Point, bgi::rstar<16>> rtree;
    for (const auto& p : tour) {
        rtree.insert(p);
    }

    // Verify each circle contains at least one point from the tour
    for (const auto& c : circles) {
        std::vector<Point> candidates;
        rtree.query(bgi::within(circleBox(c)), std::back_inserter(candidates));
        if (candidates.empty()) {
            DBG("Circle with center (" << bg::get<0>(c.center) << ", " << bg::get<1>(c.center) 
                << ") and radius " << c.r << " does not contain any tour points.");
            return false;
        }
    }

    return true;
}

double totalTourDistance(const std::vector<Point>& tour) {
    if (tour.empty()) return 0.0;

    double totalDist = 0.0;
    for (size_t i = 0; i < tour.size(); ++i) {
        const Point& p1 = tour[i];
        const Point& p2 = tour[(i + 1) % tour.size()]; // Wrap around to the start
        totalDist += bg::distance(p1, p2);
    }
    return totalDist;
}

std::vector<Point> CETSP(const std::vector<Circle>& circles, int numRepeats) {
    auto circlesCopy = circles; // Make a copy to avoid modifying the input

    // Initialize the best tour and distance
    double bestDist = std::numeric_limits<double>::infinity();
    std::vector<Point> bestTour;

    // Step 1: Remove covering circles to simplify the problem
    removeCoveringCircles(circlesCopy);

    // If there's only one circle left, directly return its center as the tour
    if (circlesCopy.size() == 1) {
        DBG("Only one circle remains after removing coverings. Returning its center.");
        return {circlesCopy[0].center};
    }

    // Repeat for some number of iterations
    for (int iter = 0; iter < numRepeats; ++iter) {
        // Step 2: Build the merge tree from the remaining circles
        std::vector<TreeNode> mergeTree = buildMergeTree(circlesCopy);

        // Step 3: Generate the tour from the merge tree
        std::vector<Point> tour = unmerge(mergeTree);

        // Step 4: Verify the generated tour to ensure correctness and update the best tour if valid
        if (verifyTour(tour, circles)) {
            double curDist = totalTourDistance(tour);
            if (curDist < bestDist) {
                bestDist = curDist;
                bestTour = tour;
            }
        }
    }

    if (bestTour.empty()) {
        DBG("Tour verification failed.");
    }
    else {
        DBG("Tour verification passed. Number of elements: " << circles.size());
    }

    // Return the valid tour
    return bestTour;
}