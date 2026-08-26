#pragma once

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>

// Cartesian 2D point used by the public API.
using Point = boost::geometry::model::point<
    double,
    2,
    boost::geometry::cs::cartesian>;

// Geometric CETSP neighborhood. Merge and tour state belongs to the solver,
// not to caller-owned input data.
struct Circle {
    Point center{};
    double r{};
};
