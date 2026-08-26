#include "lkh_parser.hpp"

#include <cstddef>
#include <iostream>
#include <sstream>
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

void expectRejected(
    std::string_view contents,
    std::size_t expectedNodeCount,
    std::string_view message) {
    std::istringstream input{std::string(contents)};
    bool rejected = false;
    try {
        static_cast<void>(
            cetsp::postprocess::detail::parseLkhTour(
                input, expectedNodeCount));
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    expect(rejected, message);
}

void testValidTour() {
    std::istringstream input{
        "NAME : sample\n"
        "TYPE : TOUR\n"
        "DIMENSION : 3\n"
        "TOUR_SECTION\n"
        "3 1\n"
        "2\n"
        "-1\n"
        "EOF\n"};

    const std::vector<std::size_t> expected{2, 0, 1};
    expect(cetsp::postprocess::detail::parseLkhTour(
               input, expected.size()) == expected,
           "a valid LKH permutation is converted to zero-based IDs");

    std::istringstream emptyTour{"TOUR_SECTION -1 EOF"};
    expect(cetsp::postprocess::detail::parseLkhTour(emptyTour, 0).empty(),
           "an explicitly empty LKH tour is accepted for an empty instance");
}

void testInvalidIds() {
    expectRejected("TOUR_SECTION 0 -1", 1,
                   "zero cannot underflow into an unsigned node ID");
    expectRejected("TOUR_SECTION -2 -1", 1,
                   "negative values other than the terminator are rejected");
    expectRejected("TOUR_SECTION 1 2 4 -1", 3,
                   "node IDs above the expected dimension are rejected");
    expectRejected("TOUR_SECTION 1junk -1", 1,
                   "partially numeric node IDs are rejected");
    expectRejected("TOUR_SECTION 9223372036854775808 -1", 1,
                   "node IDs outside the signed parser range are rejected");
}

void testInvalidTourStructure() {
    expectRejected("NAME sample EOF", 1,
                   "a missing TOUR_SECTION is rejected");
    expectRejected("TOUR_SECTION 1 EOF", 1,
                   "a missing signed terminator is rejected");
    expectRejected("TOUR_SECTION 1 1 2 -1", 3,
                   "duplicate node IDs are rejected");
    expectRejected("TOUR_SECTION 1 2 -1", 3,
                   "an incomplete permutation is rejected");
}

} // namespace

int main() {
    testValidTour();
    testInvalidIds();
    testInvalidTourStructure();

    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }

    std::cout << "All LKH parser tests passed\n";
    return 0;
}
