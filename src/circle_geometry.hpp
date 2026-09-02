#pragma once

#include <cetsp/types.hpp>

#include <boost/geometry/geometries/box.hpp>

#include <random>

// Circle geometry shared across solver translation units.
using Box = boost::geometry::model::box<Point>;

// Create a bounding box around a circle. Each bound is expanded outward by one
// representable double so points on a computed boundary remain queryable.
Box circleBox(const Circle& c, double padding = 0.0);

// Compute "gap" distance (center distance minus the sum of radii).
double gapDist(const Circle& first, const Circle& second);

// Choose a low-cost insertion point inside a circle. Intersecting edges get a
// zero-cost point; otherwise one safeguarded Newton step refines a boundary
// point selected from deterministic geometric seeds.
[[nodiscard]] Point chooseInsertionPoint(
    const Circle& neighborhood,
    const Point& edgeStart,
    const Point& edgeEnd);

// Compute a combined "representative" circle from two circles
Circle makeCombinedCircle(
    const Circle& first,
    const Circle& second,
    std::mt19937_64& randomEngine);
