#pragma once

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace cetsp::postprocess::detail {

inline std::vector<std::size_t> parseLkhTour(
    std::istream& input,
    std::size_t expectedNodeCount) {
    std::vector<std::size_t> tour;
    tour.reserve(expectedNodeCount);
    std::vector<bool> seen(expectedNodeCount, false);

    std::string token;
    bool readingTour = false;
    bool terminated = false;
    while (input >> token) {
        if (!readingTour) {
            if (token == "TOUR_SECTION") {
                readingTour = true;
            }
            continue;
        }
        if (token == "EOF") {
            break;
        }

        std::int64_t externalId = 0;
        const char* const begin = token.data();
        const char* const end = begin + token.size();
        const auto [parsedEnd, error] =
            std::from_chars(begin, end, externalId);
        if (error != std::errc{} || parsedEnd != end) {
            throw std::runtime_error(
                "invalid LKH tour node ID: " + token);
        }
        if (externalId == -1) {
            terminated = true;
            break;
        }
        if (externalId <= 0 ||
            static_cast<std::uint64_t>(externalId) > expectedNodeCount) {
            throw std::runtime_error(
                "LKH tour node ID is outside the expected range: " + token);
        }

        const std::size_t nodeId =
            static_cast<std::size_t>(externalId - 1);
        if (seen[nodeId]) {
            throw std::runtime_error(
                "LKH tour contains a duplicate node ID: " + token);
        }
        seen[nodeId] = true;
        tour.push_back(nodeId);
    }

    if (input.bad()) {
        throw std::runtime_error("failed while reading LKH tour output");
    }
    if (!readingTour) {
        throw std::runtime_error("LKH output has no TOUR_SECTION");
    }
    if (!terminated) {
        throw std::runtime_error("LKH TOUR_SECTION has no -1 terminator");
    }
    if (tour.size() != expectedNodeCount) {
        throw std::runtime_error(
            "LKH tour does not contain the expected number of nodes");
    }
    return tour;
}

} // namespace cetsp::postprocess::detail
