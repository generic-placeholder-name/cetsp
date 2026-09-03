#include "circle_geometry.hpp"

#include <boost/geometry.hpp>

#include <algorithm>
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

struct AngularDerivatives {
    double first;
    double second;
};

struct BoundaryCandidate {
    double angle;
    double cost;
};

struct BoundaryMotion {
    double pointX;
    double pointY;
    double velocityX;
    double velocityY;
    double accelerationX;
    double accelerationY;
    double velocitySquared;
};

AngularDerivatives endpointDistanceDerivatives(
    const BoundaryMotion& motion,
    const Point& endpoint) {
    const double offsetX =
        motion.pointX - bg::get<0>(endpoint);
    const double offsetY =
        motion.pointY - bg::get<1>(endpoint);
    const double distance = std::hypot(offsetX, offsetY);
    const double directionalRate =
        offsetX * motion.velocityX + offsetY * motion.velocityY;
    return {
        directionalRate / distance,
        (motion.velocitySquared +
         offsetX * motion.accelerationX +
         offsetY * motion.accelerationY) /
                distance -
            directionalRate * directionalRate /
                (distance * distance * distance),
    };
}

struct BoundaryObjective {
    const Circle& neighborhood;
    const Point& edgeStart;
    const Point& edgeEnd;

    [[nodiscard]] Point pointAt(double angle) const {
        return Point{
            bg::get<0>(neighborhood.center) +
                neighborhood.r * std::cos(angle),
            bg::get<1>(neighborhood.center) +
                neighborhood.r * std::sin(angle)};
    }

    [[nodiscard]] BoundaryCandidate evaluate(double angle) const {
        const Point point = pointAt(angle);
        return {
            angle,
            bg::distance(edgeStart, point) +
                bg::distance(point, edgeEnd),
        };
    }

    [[nodiscard]] AngularDerivatives derivativesAt(double angle) const {
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        const Point point = pointAt(angle);
        const double radius = neighborhood.r;
        const BoundaryMotion motion{
            bg::get<0>(point),
            bg::get<1>(point),
            -radius * sine,
            radius * cosine,
            -radius * cosine,
            -radius * sine,
            radius * radius,
        };

        const AngularDerivatives start =
            endpointDistanceDerivatives(motion, edgeStart);
        const AngularDerivatives end =
            endpointDistanceDerivatives(motion, edgeEnd);
        return {
            start.first + end.first,
            start.second + end.second,
        };
    }
};

BoundaryCandidate chooseInitialBoundaryCandidate(
    const BoundaryObjective& objective,
    const Point& closestApproach) {
    const double firstBearing =
        bearingFrom(objective.neighborhood.center, objective.edgeStart);
    const double secondBearing =
        bearingFrom(objective.neighborhood.center, objective.edgeEnd);
    const double bisectorX =
        std::cos(firstBearing) + std::cos(secondBearing);
    const double bisectorY =
        std::sin(firstBearing) + std::sin(secondBearing);
    const double closestBearing =
        bearingFrom(objective.neighborhood.center, closestApproach);
    const double bisectorBearing =
        std::hypot(bisectorX, bisectorY) >
                32.0 * std::numeric_limits<double>::epsilon()
            ? std::atan2(bisectorY, bisectorX)
            : closestBearing;

    // Compare the edge's closest approach with the circular mean of the two
    // endpoint directions.
    const BoundaryCandidate closest = objective.evaluate(closestBearing);
    const BoundaryCandidate bisector = objective.evaluate(bisectorBearing);
    return bisector.cost < closest.cost ? bisector : closest;
}

BoundaryCandidate refineBoundaryCandidateOnce(
    const BoundaryObjective& objective,
    const BoundaryCandidate& candidate) {
    const AngularDerivatives derivatives =
        objective.derivativesAt(candidate.angle);
    const double curvatureTolerance =
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(candidate.cost));
    if (!std::isfinite(derivatives.first) ||
        !std::isfinite(derivatives.second) ||
        derivatives.second <= curvatureTolerance) {
        return candidate;
    }

    // One refinement only. The modulo preserves the point represented by a
    // large Newton step while keeping trigonometric range reduction tame.
    const double rawStep = -derivatives.first / derivatives.second;
    if (!std::isfinite(rawStep)) {
        return candidate;
    }

    const BoundaryCandidate refined = objective.evaluate(
        candidate.angle + std::remainder(rawStep, fullTurn));
    return refined.cost < candidate.cost ? refined : candidate;
}

Point chooseBoundaryInsertionPoint(
    const BoundaryObjective& objective,
    const Point& closestApproach) {
    const BoundaryCandidate initial =
        chooseInitialBoundaryCandidate(objective, closestApproach);
    const BoundaryCandidate refined =
        refineBoundaryCandidateOnce(objective, initial);
    return objective.pointAt(refined.angle);
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

double gapDist(const Circle& first, const Circle& second) {
    return bg::distance(first.center, second.center) - (first.r + second.r);
}

Point chooseInsertionPoint(
    const Circle& neighborhood,
    const Point& edgeStart,
    const Point& edgeEnd) {
    const Point& center = neighborhood.center;
    const double radius = neighborhood.r;
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

    const BoundaryObjective objective{
        neighborhood,
        edgeStart,
        edgeEnd,
    };
    return chooseBoundaryInsertionPoint(objective, closestEdgePoint);
}

Circle makeCombinedCircle(
    const Circle& first,
    const Circle& second,
    std::mt19937_64& randomEngine) {
    // Unpack the circles for ease of working with their components.
    const Point& p1 = first.center;
    const double r1 = first.r;
    const Point& p2 = second.center;
    const double r2 = second.r;

    // Basically, a circle approximating the "lens" shape formed by the intersection of two circles.
    // This choice is largely arbitrary but seems to work better than just taking the midpoint.
    double x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
    double x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
    double dx = x2 - x1, dy = y2 - y1;
    double d = std::hypot(dx, dy);

    // If one circle is completely inside the other
    if (d + std::min(r1, r2) <= std::max(r1, r2)) {
        return (r1 < r2) ? first : second;
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
        return Circle{center, 0.0};
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

    return Circle{center, radius};
}
