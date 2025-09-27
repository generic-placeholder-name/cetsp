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

// Normalize points so that average pairwise distance ≈ scale_target
void normalize_points(std::vector<Point>& points, double scale_target = 1048576) {
    double total = 0.0;
    size_t count = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            total += bg::distance(points[i], points[j]);
            ++count;
        }
    }
    if (count == 0) return;
    double avg = total / count;
    double scale = scale_target / avg;
    for (auto& pt : points) {
        bg::set<0>(pt, bg::get<0>(pt) * scale);
        bg::set<1>(pt, bg::get<1>(pt) * scale);
    }
}

// Write .tsp file in TSPLIB EUC_2D format
void write_tsp_file(const std::string& path, const std::vector<Point>& points) {
    std::ofstream ofs(path);
    ofs << "NAME : TSP\n";
    ofs << "TYPE : TSP\n";
    ofs << "DIMENSION : " << points.size() << "\n";
    ofs << "EDGE_WEIGHT_TYPE : EUC_2D\n";
    ofs << "NODE_COORD_SECTION\n";
    for (size_t i = 0; i < points.size(); ++i) {
        ofs << i + 1 << " " << bg::get<0>(points[i]) << " " << bg::get<1>(points[i]) << "\n";
    }
    ofs << "EOF\n";
}

// Write LKH parameter file
void write_par_file(const std::string& path, const std::string& tsp_path, const std::string& tour_path) {
    std::ofstream ofs(path);
    ofs << "PROBLEM_FILE = " << tsp_path << "\n";
    ofs << "OUTPUT_TOUR_FILE = " << tour_path << "\n";
    ofs << "RUNS = 1\n";
}

// Run LKH on the given parameter file
void run_lkh(const std::string& lkh_bin, const std::string& par_path) {
    std::string command = lkh_bin + " " + par_path;
    std::cerr << "[LKH] Command: " << command << "\n";
    int result = std::system(command.c_str());
    if (result != 0) {
        throw std::runtime_error("LKH failed to run.");
    }
}

// Parse LKH tour output
std::vector<size_t> parse_tour_file(const std::string& path) {
    std::ifstream ifs(path);
    std::vector<size_t> tour;
    std::string line;
    bool reading = false;
    while (std::getline(ifs, line)) {
        if (line == "TOUR_SECTION") {
            reading = true;
            continue;
        }
        if (!reading) continue;
        if (line == "-1" || line == "EOF") break;
        tour.push_back(std::stoul(line) - 1);  // convert from 1-based to 0-based
    }
    return tour;
}

// This will be a build option parameter, but for now we use a default path 
const std::string default_lkh_path = R"(C:\Users\minhk\Documents\Code\LKH\x64\Release\LKH-3.exe)";

// Main TSP solver function
std::vector<size_t> solve_tsp_with_lkh(std::vector<Point> points, const std::string& lkh_bin = default_lkh_path) {
    if (points.size() < 2)
        return {};

    std::cerr << "[LKH] Normalizing points...\n";
    normalize_points(points);

    std::string tmp_dir = std::filesystem::temp_directory_path().string();
    std::string tsp_path = tmp_dir + "lkh_problem.tsp";
    std::string par_path = tmp_dir + "lkh_params.par";
    std::string tour_path = tmp_dir + "lkh_output.tour";

    std::cerr << "[LKH] Writing TSP file: " << tsp_path << "\n";
    write_tsp_file(tsp_path, points);
    std::cerr << "[LKH] Writing parameter file: " << par_path << "\n";
    write_par_file(par_path, tsp_path, tour_path);
    std::cerr << "[LKH] Running LKH: " << lkh_bin << " " << par_path << "\n";
    run_lkh(lkh_bin, par_path);
    std::cerr << "[LKH] Parsing tour file: " << tour_path << "\n";

    return parse_tour_file(tour_path);
}
