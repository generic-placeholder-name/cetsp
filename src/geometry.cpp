#include "geometry.hpp"

Box circleBox(const Circle& c) {
    double adjr = c.r + 1e-12; // Add a small epsilon to avoid precision issues
    return Box(
        Point(bg::get<0>(c.center) - adjr, bg::get<1>(c.center) - adjr),
        Point(bg::get<0>(c.center) + adjr, bg::get<1>(c.center) + adjr)
    );
}

double gapDist(const Point& a, const Point& b, double ra, double rb) {
    return bg::distance(a, b) - (ra + rb);
}

std::pair<Point, double> makeCombinedCircle(const Point& p1, double r1, const Point& p2, double r2) {
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

    // Radius range: [overlapDepth, h]
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(overlapDepth, h);
    double radius = dis(gen);

    return std::make_pair(center, radius);
}