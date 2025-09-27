#pragma once

#include <utility>
#include <random>
#include <cmath>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "structures.hpp"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Aliases
typedef bg::model::segment<Point> Segment;
typedef bg::model::box<Point> Box;

typedef std::pair<Point, size_t> PointValue;
typedef std::pair<Box, size_t> BoxValue;
typedef std::pair<Segment, size_t> SegmentValue;

// ==================== Geometry helpers ====================

// Create bounding box around a circle
Box circleBox(const Circle& c);

// Compute "gap" distance (center-dist minus sum of radii)
double gapDist(const Point& a, const Point& b, double ra, double rb);

// Compute a combined "representative" circle from two circles
std::pair<Point, double> makeCombinedCircle(
    const Point& p1, double r1,
    const Point& p2, double r2
);

