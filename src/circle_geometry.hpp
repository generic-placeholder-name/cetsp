#pragma once

#include <cetsp/types.hpp>

#include <boost/geometry/geometries/box.hpp>

#include <random>
#include <utility>

// Circle geometry shared across solver translation units.
using Box = boost::geometry::model::box<Point>;

// Create a bounding box around a circle. Each bound is expanded outward by one
// representable double so points on a computed boundary remain queryable.
Box circleBox(const Circle& c, double padding = 0.0);

// Compute "gap" distance (center-dist minus sum of radii)
double gapDist(const Point& a, const Point& b, double ra, double rb);

// Choose a low-cost insertion point inside a circle. Intersecting edges get a
// zero-cost point; otherwise one safeguarded Newton step refines a boundary
// point selected from deterministic geometric seeds.
[[nodiscard]] Point chooseInsertionPoint(
    const Point& center,
    double radius,
    const Point& edgeStart,
    const Point& edgeEnd);

// Compute a combined "representative" circle from two circles
std::pair<Point, double> makeCombinedCircle(
    const Point& p1, double r1,
    const Point& p2, double r2,
    std::mt19937_64& randomEngine
);
