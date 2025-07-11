#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <cmath>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>

namespace bg = boost::geometry;
using Point = bg::model::point<double, 2, bg::cs::cartesian>;

#include "LKH.hpp"

double totalTourDistance(const std::vector<Point>& tour) {
    if (tour.empty()) return 0.0;

    double totalDist = 0.0;
    for (size_t i = 0; i < tour.size(); ++i) {
        const Point& p1 = tour[i];
        const Point& p2 = tour[(i + 1) % tour.size()]; // Wrap around to the start
        totalDist += bg::distance(p1, p2);
    }
    return totalDist;
}

std::vector<Point> read_points_from_file(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y;
        if (!(iss >> x >> y)) continue;
        points.emplace_back(x, y);
    }
    return points;
}

int main() {
    const std::string filename = "../data/cmu_output/rotatingDiamonds4.txt";

    try {
        // Get absolute path for input file
        std::filesystem::path input_path = std::filesystem::absolute(filename);
        auto points = read_points_from_file(input_path.string());
        if (points.empty()) {
            std::cerr << "No points read from file.\n";
            return 1;
        }

        double origDist = totalTourDistance(points);
        std::cout << "Original tour distance: " << origDist << "\n";

        // Pass the absolute path to LKH options if needed
        auto tour_indices = solve_tsp_with_lkh(points);

        std::vector<Point> new_tour;
        for (size_t idx : tour_indices) {
            if (idx >= points.size()) throw std::runtime_error("Invalid tour index");
            new_tour.push_back(points[idx]);
        }
        double newDist = totalTourDistance(new_tour);
        std::cout << "LKH tour distance: " << newDist << "\n";
        std::cout << "Delta: " << origDist - newDist << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
