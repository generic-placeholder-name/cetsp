#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include <cetsp/types.hpp>

bool verifyTour(const std::vector<Point>& tour, const std::vector<Circle>& circles);

double totalTourDistance(const std::vector<Point>& tour);

struct CetspOptions {
    int numRepeats = 10;
    std::optional<std::uint64_t> seed;
    // Zero selects hardware concurrency; positive values are capped by numRepeats.
    std::size_t maxThreads = 0;
};

// Throws std::invalid_argument for invalid circle geometry or a non-positive
// repetition count. An empty instance with a positive count returns an empty tour.
std::vector<Point> solveCetsp(
    const std::vector<Circle>& circles,
    const CetspOptions& options);

// Compatibility entry point. Without an explicit seed, each call uses fresh entropy.
std::vector<Point> CETSP(const std::vector<Circle>& circles, int numRepeats = 10);
