#include "circle_geometry.hpp"

#include <boost/geometry.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace bg = boost::geometry;

Box circleBox(const Circle& c, double padding) {
    const double radius = c.r + padding;
    const double infinity = std::numeric_limits<double>::infinity();
    const double centerX = bg::get<0>(c.center);
    const double centerY = bg::get<1>(c.center);

    return Box(
        Point(std::nextafter(centerX - radius, -infinity),
              std::nextafter(centerY - radius, -infinity)),
        Point(std::nextafter(centerX + radius, infinity),
              std::nextafter(centerY + radius, infinity))
    );
}

double gapDist(const Point& a, const Point& b, double ra, double rb) {
    return bg::distance(a, b) - (ra + rb);
}

std::pair<Point, double> makeCombinedCircle(
    const Point& p1,
    double r1,
    const Point& p2,
    double r2,
    std::mt19937_64& randomEngine) {
    // Basically, a circle approximating the "lens" shape formed by the intersection of two circles.
    // This choice is largely arbitrary but seems to work better than just taking the midpoint.
    double x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
    double x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
    double dx = x2 - x1, dy = y2 - y1;
    double d = std::hypot(dx, dy);

    // If one circle is completely inside the other
    if (d + std::min(r1, r2) <= std::max(r1, r2)) {
        return (r1 < r2) ? std::make_pair(p1, r1) : std::make_pair(p2, r2);
    }

    // Normalize direction vector
    double ux = dx / d;
    double uy = dy / d;

    // Points on the edges of the circles along the line connecting centers
    double xEdge1 = x1 + ux * r1;
    double yEdge1 = y1 + uy * r1;
    double xEdge2 = x2 - ux * r2;
    double yEdge2 = y2 - uy * r2;

    // Midpoint of these edge points is the new center
    Point center((xEdge1 + xEdge2) / 2.0, (yEdge1 + yEdge2) / 2.0);

    // Handle non-overlapping case (set radius = 0)
    if (d >= r1 + r2) {
        return std::make_pair(center, 0.0);
    }

    // Circles intersect — compute potential radius range
    double overlapDepth = (r1 + r2 - d) / 2.0;

    // Compute distance from center of first circle to chord midpoint
    double a = (r1*r1 - r2*r2 + d*d) / (2*d);
    double h = std::sqrt(std::max(0.0, r1*r1 - a*a));  // half chord length

    // Pick between the axial overlap depth and the half-chord length. Depending
    // on the circles' relative sizes, either geometric quantity can be larger.
    const auto [minimumRadius, maximumRadius] = std::minmax(overlapDepth, h);
    std::uniform_real_distribution<> dis(minimumRadius, maximumRadius);
    double radius = dis(randomEngine);

    return std::make_pair(center, radius);
}
