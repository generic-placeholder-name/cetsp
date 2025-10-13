#include "CETSP.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <filesystem>
#include <optional>

// Simple settings structure
struct Settings {
    std::filesystem::path inputDir = "./data/cmu_processed";
    std::filesystem::path outputDir = "./data/cmu_output";
    std::vector<std::string> filesToRead; // empty = process all
    std::optional<std::filesystem::path> logFile; // optional log file
    int numRepeats = 10; // Number of repetitions for CETSP
};

// Read settings from a file (very simple: key value per line)
Settings loadSettings(const std::string& settingsFile) {
    Settings s;
    std::ifstream in(settingsFile);
    if (!in) {
        std::cerr << "No settings file found. Using defaults." << std::endl;
        return s;
    }

    std::string key;
    while (in >> key) {
        if (key == "inputDir") {
            in >> s.inputDir;
        } else if (key == "outputDir") {
            in >> s.outputDir;
        } else if (key == "file") {
            std::string fname;
            in >> fname;
            s.filesToRead.push_back(fname);
        } else if (key == "logFile") {
            std::string path;
            in >> path;
            s.logFile = path;
        } else if (key == "numRepeats") {
            in >> s.numRepeats;
        }
    }
    return s;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Load settings
    Settings settings = loadSettings("settings.txt");
    std::cout << "Settings:\n";
    std::cout << "  Input directory: " << settings.inputDir << "\n";
    std::cout << "  Output directory: " << settings.outputDir << "\n";
    if (!settings.filesToRead.empty()) {
        std::cout << "  Files to read:\n";
        for (const auto& f : settings.filesToRead) {
            std::cout << "    " << f << "\n";
        }
    } else {
        std::cout << "  Processing all files in input directory.\n";
    }
    std::cout << "  Log file: " << (settings.logFile ? settings.logFile->string() : "stdout") << "\n";
    std::cout << "  Number of repetitions: " << settings.numRepeats << "\n";

    // Set up logging
    std::ofstream logFile;
    std::ostream& logStream = [&]() -> std::ostream& {
        if (settings.logFile) {
            logFile.open(*settings.logFile); // overwrite mode
            if (logFile) {
                return logFile;
            } else {
                std::cerr << "Failed to open log file. Falling back to stdout.\n";
            }
        }
        return std::cout;
    }();

    namespace fs = std::filesystem;
    if (!fs::exists(settings.outputDir)) {
        fs::create_directory(settings.outputDir);
    }
    logStream << "Input directory: " << settings.inputDir << "\n";
    logStream << "Output directory: " << settings.outputDir << "\n";

    // Iterate through all files in the input directory
    for (const auto& entry : fs::directory_iterator(settings.inputDir)) {
        if (!settings.filesToRead.empty()) {
            if (std::find(settings.filesToRead.begin(),
                          settings.filesToRead.end(),
                          entry.path().filename().string()) == settings.filesToRead.end()) {
                continue;
            }
        }
        if (!entry.is_regular_file()) continue;

        const fs::path inputFile = entry.path();
        const fs::path outputFile = settings.outputDir / inputFile.filename();

        logStream << "Processing file: " << inputFile << "\n";

        auto startTime = std::chrono::high_resolution_clock::now();

        // Read circles from the input file
        std::vector<Circle> circles;
        std::ifstream inFile(inputFile);
        if (!inFile) {
            logStream << "Failed to open input file: " << inputFile << "\n";
            continue;
        }

        double x, y, r;
        size_t nodeId = 0;
        while (inFile >> x >> y >> r) {
            circles.push_back({Point(x, y), r, nodeId++, static_cast<size_t>(-1), 0.0, {}});
        }
        inFile.close();

        // Generate the tour using the CETSP algorithm
        auto tour = CETSP(circles, settings.numRepeats);
        if (tour.empty()) {
            logStream << "Failed to generate a valid tour for file: " << inputFile << "\n";
            continue;
        }

        // Write the tour points to the output file
        std::ofstream outFile(outputFile);
        if (!outFile) {
            logStream << "Failed to open output file: " << outputFile << "\n";
            continue;
        }

        double dist = totalTourDistance(tour);
        outFile << std::fixed << std::setprecision(4);
        outFile << dist << "\n";
        for (const auto& p : tour) {
            outFile << bg::get<0>(p) << " " << bg::get<1>(p) << "\n";
        }
        outFile.close();

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        logStream << "Finished processing file: " << inputFile
                  << " in " << duration << " ms, distance=" << dist << "\n";
    }

    return 0;
}
