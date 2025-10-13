#pragma once

#include <vector>

#include "structures.hpp"
#include "geometry.hpp"

// Function signatures
bool verifyTour(const std::vector<Point>& tour, const std::vector<Circle>& circles);

double totalTourDistance(const std::vector<Point>& tour);

std::vector<Point> CETSP(const std::vector<Circle>& circles, int numRepeats = 10);
