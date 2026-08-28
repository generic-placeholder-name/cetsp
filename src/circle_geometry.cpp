#include "circle_geometry.hpp"

#include <boost/geometry.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace bg = boost::geometry;

namespace {

constexpr double fullTurn = 2.0 * std::numbers::pi;

Point closestPointOnSegment(
    const Point& point,
    const Point& segmentStart,
    const Point& segmentEnd) {
    const double startX = bg::get<0>(segmentStart);
    const double startY = bg::get<1>(segmentStart);
    const double dx = bg::get<0>(segmentEnd) - startX;
    const double dy = bg::get<1>(segmentEnd) - startY;
    const double lengthSquared = dx * dx + dy * dy;
    if (lengthSquared == 0.0) {
        return segmentStart;
    }

    const double offsetX = bg::get<0>(point) - startX;
    const double offsetY = bg::get<1>(point) - startY;
    const double fraction = std::clamp(
        (offsetX * dx + offsetY * dy) / lengthSquared,
        0.0,
        1.0);
    return Point{startX + fraction * dx, startY + fraction * dy};
}

double bearingFrom(const Point& origin, const Point& target) {
    return std::atan2(
        bg::get<1>(target) - bg::get<1>(origin),
        bg::get<0>(target) - bg::get<0>(origin));
}

struct BoundaryPathObjective {
    const Point& center;
    double radius;
    const Point& edgeStart;
    const Point& edgeEnd;

    [[nodiscard]] Point pointAt(double angle) const {
        return Point{
            bg::get<0>(center) + radius * std::cos(angle),
            bg::get<1>(center) + radius * std::sin(angle)};
    }

    [[nodiscard]] double costAt(double angle) const {
        const Point point = pointAt(angle);
        return bg::distance(edgeStart, point) +
               bg::distance(point, edgeEnd);
    }

    struct AngularDerivatives {
        double first;
        double second;
    };

    [[nodiscard]] AngularDerivatives derivativesAt(double angle) const {
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        const Point point = pointAt(angle);
        const double pointX = bg::get<0>(point);
        const double pointY = bg::get<1>(point);
        const double velocityX = -radius * sine;
        const double velocityY = radius * cosine;
        const double accelerationX = -radius * cosine;
        const double accelerationY = -radius * sine;
        const double velocitySquared = radius * radius;

        AngularDerivatives result{0.0, 0.0};
        const std::array<const Point*, 2> endpoints{
            &edgeStart,
            &edgeEnd,
        };
        for (const Point* endpoint : endpoints) {
            const double offsetX = pointX - bg::get<0>(*endpoint);
            const double offsetY = pointY - bg::get<1>(*endpoint);
            const double distance = std::hypot(offsetX, offsetY);
            const double directionalRate =
                offsetX * velocityX + offsetY * velocityY;
            result.first += directionalRate / distance;
            result.second +=
                (velocitySquared +
                 offsetX * accelerationX +
                 offsetY * accelerationY) /
                    distance -
                directionalRate * directionalRate /
                    (distance * distance * distance);
        }
        return result;
    }
};

struct BoundarySample {
    double angle;
    double cost;
};

Point chooseBoundaryPoint(
    const BoundaryPathObjective& objective,
    const Point& closestEdgePoint) {
    const double firstBearing =
        bearingFrom(objective.center, objective.edgeStart);
    const double secondBearing =
        bearingFrom(objective.center, objective.edgeEnd);
    const double bisectorX =
        std::cos(firstBearing) + std::cos(secondBearing);
    const double bisectorY =
        std::sin(firstBearing) + std::sin(secondBearing);
    const double closestBearing =
        bearingFrom(objective.center, closestEdgePoint);
    const double bisectorBearing =
        std::hypot(bisectorX, bisectorY) >
                32.0 * std::numeric_limits<double>::epsilon()
            ? std::atan2(bisectorY, bisectorX)
            : closestBearing;

    const std::array<double, 2> seedAngles{
        closestBearing,
        bisectorBearing,
    };
    // These cover the two useful geometric views of the edge: its closest
    // approach to the circle and the average direction of its endpoints.
    BoundarySample best{
        seedAngles.front(),
        objective.costAt(seedAngles.front()),
    };
    for (std::size_t index = 1; index < seedAngles.size(); ++index) {
        const double angle = seedAngles[index];
        const double cost = objective.costAt(angle);
        if (cost < best.cost) {
            best = {angle, cost};
        }
    }

    const BoundaryPathObjective::AngularDerivatives derivatives =
        objective.derivativesAt(best.angle);
    const double curvatureTolerance =
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(best.cost));
    if (std::isfinite(derivatives.first) &&
        std::isfinite(derivatives.second) &&
        derivatives.second > curvatureTolerance) {
        // One refinement only. The modulo preserves the point represented by
        // a large Newton step while keeping trigonometric range reduction tame.
        const double rawStep = -derivatives.first / derivatives.second;
        if (std::isfinite(rawStep)) {
            const double refinedAngle =
                best.angle + std::remainder(rawStep, fullTurn);
            const double refinedCost = objective.costAt(refinedAngle);
            if (refinedCost < best.cost) {
                best = {refinedAngle, refinedCost};
            }
        }
    }

    return objective.pointAt(best.angle);
}

} // namespace

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

Point chooseInsertionPoint(
    const Point& center,
    double radius,
    const Point& edgeStart,
    const Point& edgeEnd) {
    if (radius < 0.0) {
        throw std::invalid_argument(
            "insertion-point circle radius must be non-negative");
    }
    if (radius == 0.0) {
        return center;
    }

    const Point closestEdgePoint =
        closestPointOnSegment(center, edgeStart, edgeEnd);
    if (bg::distance(center, closestEdgePoint) <= radius) {
        return closestEdgePoint;
    }

    const BoundaryPathObjective objective{
        center,
        radius,
        edgeStart,
        edgeEnd,
    };
    return chooseBoundaryPoint(objective, closestEdgePoint);
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
