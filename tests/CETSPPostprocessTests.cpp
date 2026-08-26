#include <cetsp/cetsp.hpp>
#include <cetsp/postprocess/pipeline.hpp>

#include <boost/geometry/core/access.hpp>

#include <cstddef>
#include <iostream>
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

bool sameTour(
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

void testDisabledPipelineIsIdentity() {
    const std::vector<Circle> circles{
        {Point{0.0, 0.0}, 1.0},
        {Point{4.0, 0.0}, 1.0},
    };
    const std::vector<Point> tour{
        Point{0.0, 0.0}, Point{4.0, 0.0}};
    cetsp::postprocess::PostprocessOptions options;
    const auto result = cetsp::postprocess::polishTour(
        tour, circles, options);
    expect(sameTour(result, tour),
           "a pipeline with no enabled backend preserves the tour");
}

void testInvalidRoundCountIsRejected() {
    cetsp::postprocess::PostprocessOptions options;
    options.rounds = 0;
    bool threw = false;
    try {
        static_cast<void>(cetsp::postprocess::polishTour({}, {}, options));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "zero postprocessing rounds are rejected");
}

void testSocpDropsAnUnneededVisitWithoutLaunchingSolver() {
    const std::vector<Circle> circles{{Point{0.0, 0.0}, 1.0}};
    const std::vector<Point> tour{
        Point{0.0, 0.0}, Point{10.0, 10.0}};
    cetsp::postprocess::SocpOptions options;
    const auto result = cetsp::postprocess::improvePointsWithSocp(
        tour, circles, options);
    expect(result.size() == 1 && verifyTour(result, circles),
           "SOCP preprocessing removes visits assigned no circle");
}

} // namespace

int main() {
    testDisabledPipelineIsIdentity();
    testInvalidRoundCountIsRejected();
    testSocpDropsAnUnneededVisitWithoutLaunchingSolver();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }
    std::cout << "All CETSP postprocess tests passed\n";
    return 0;
}
