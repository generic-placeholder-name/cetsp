#include <cetsp/cetsp.hpp>

#include <boost/geometry/core/access.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

int failures = 0;

void expect(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void expectNear(double actual, double expected, double tolerance, std::string_view message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr << "FAIL: " << message << " (actual=" << actual
                  << ", expected=" << expected << ")\n";
        ++failures;
    }
}

Circle makeCircle(double x, double y, double radius) {
    return {Point{x, y}, radius};
}

bool exactlyEqualTours(
    const std::vector<Point>& first,
    const std::vector<Point>& second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (std::size_t index = 0; index < first.size(); ++index) {
        if (boost::geometry::get<0>(first[index]) !=
                boost::geometry::get<0>(second[index]) ||
            boost::geometry::get<1>(first[index]) !=
                boost::geometry::get<1>(second[index])) {
            return false;
        }
    }
    return true;
}

void testTourDistance() {
    expectNear(totalTourDistance({}), 0.0, 0.0, "an empty tour has zero length");
    expectNear(totalTourDistance({Point{2.0, 3.0}}), 0.0, 0.0,
               "a singleton tour has zero length");

    const std::vector<Point> triangle{
        Point{0.0, 0.0}, Point{3.0, 0.0}, Point{3.0, 4.0}};
    expectNear(totalTourDistance(triangle), 12.0, 1e-12,
               "tour distance includes the closing edge");
}

void testExactTourVerification() {
    const std::vector<Circle> circles{makeCircle(0.0, 0.0, 1.0)};
    const double infinity = std::numeric_limits<double>::infinity();

    expect(verifyTour({Point{1.0, 0.0}}, circles),
           "a point on a circle boundary covers that neighborhood");
    expect(verifyTour({Point{std::nextafter(1.0, infinity), 0.0}}, circles),
           "broad and narrow verification agree just beyond a rounded boundary");
    expect(!verifyTour({Point{0.9, 0.9}}, circles),
           "a point inside the bounding box but outside the disk is rejected");
    expect(!verifyTour({}, circles), "a nonempty instance is not covered by an empty tour");
    expect(verifyTour({}, {}), "the empty CETSP instance is vacuously covered");
}

void testBoundaryInputs() {
    expect(CETSP({}).empty(), "an empty instance returns an empty tour");

    const std::vector<Circle> singleton{makeCircle(2.0, -4.0, 3.0)};
    const auto singletonTour = CETSP(singleton, 1);
    expect(singletonTour.size() == 1, "one circle produces one tour point");
    expect(verifyTour(singletonTour, singleton), "the singleton tour covers its circle");

    bool threw = false;
    try {
        static_cast<void>(CETSP(singleton, 0));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "a non-positive repetition count is rejected");

    threw = false;
    try {
        static_cast<void>(CETSP({}, 0));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "invalid options are rejected even for an empty instance");

    threw = false;
    try {
        static_cast<void>(CETSP({makeCircle(0.0, 0.0, -1.0)}, 1));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "a negative circle radius is rejected");
}

void testSeededRunsAreValid() {
    const std::vector<Circle> circles{
        makeCircle(-2.0, -1.0, 0.75),
        makeCircle(2.0, -1.0, 0.5),
        makeCircle(2.0, 2.0, 0.9),
        makeCircle(-2.0, 2.0, 0.6),
        makeCircle(0.0, 0.5, 0.4),
    };
    CetspOptions options;
    options.numRepeats = 8;
    options.seed = 0x9e3779b97f4a7c15ULL;
    options.maxThreads = 1;

    const auto first = solveCetsp(circles, options);
    options.maxThreads = 4;
    const auto second = solveCetsp(circles, options);

    expect(!first.empty() && verifyTour(first, circles),
           "an explicitly seeded run produces a valid tour");
    expect(!second.empty() && verifyTour(second, circles),
            "reusing a seed remains valid for every set backend");
    expect(exactlyEqualTours(first, second),
           "worker count does not change a seeded run within one process");
}

void testSmallCetspInstances() {
    const std::vector<Circle> disjoint{
        makeCircle(-5.0, 0.0, 1.0),
        makeCircle(5.0, 0.0, 1.0),
    };
    const auto tour = CETSP(disjoint, 3);
    expect(!tour.empty(), "a two-circle instance produces a tour");
    expect(verifyTour(tour, disjoint), "the two-circle tour covers both neighborhoods");

    const std::vector<Circle> nested{
        makeCircle(0.0, 0.0, 5.0),
        makeCircle(1.0, 0.0, 1.0),
    };
    const auto nestedTour = CETSP(nested, 1);
    expect(nestedTour.size() == 1,
           "nested CETSP neighborhoods need only one tour point");
    expect(verifyTour(nestedTour, nested),
           "the retained smaller neighborhood also covers the larger one");

    const std::vector<Circle> fourPointInstance{
        makeCircle(-4.0, -3.0, 0.0),
        makeCircle(4.0, -3.0, 0.0),
        makeCircle(4.0, 3.0, 0.0),
        makeCircle(-4.0, 3.0, 0.0),
    };
    const auto fourPointTour = CETSP(fourPointInstance, 4);
    expect(fourPointTour.size() == fourPointInstance.size(),
           "four disjoint point neighborhoods remain distinct tour nodes");
    expect(verifyTour(fourPointTour, fourPointInstance),
           "the multi-node reconstruction path preserves coverage");
}

} // namespace

int main() {
    testTourDistance();
    testExactTourVerification();
    testBoundaryInputs();
    testSmallCetspInstances();
    testSeededRunsAreValid();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }

    std::cout << "All public CETSP tests passed\n";
    return 0;
}
