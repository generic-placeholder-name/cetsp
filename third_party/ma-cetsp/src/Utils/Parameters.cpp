/**
 * Parameters.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "Utils/cmdline.h"
#include "Utils/Parameters.hpp"

Parameters::Parameters(int argc, char **argv) {
    cmdline::parser parser;
    parser.add<int>("random", 'g', "random number for LKH", false, 0);
    parser.add<int>("instance", 'i', "instance index", false, INSTANCE_INDEX);
    parser.add<int>("seed", 's', "seed for random generator", false, SEED);
    parser.add<std::string>("init", '\0', "initialization", false, INITIALIZATION);
    parser.add<std::string>("select", '\0', "selection", false, SELECTION);
    parser.add<std::string>("crossover", '\0', "crossover", false, CROSSOVER);
    parser.add<std::string>("improvement", '\0', "improvement", false, IMPROVEMENT);
    parser.add<std::string>("greed", '\0', "greed", false, GREEDY_ALGO);
    parser.add<std::string>("distance", '\0', "distance", false, DISTANCE);
    parser.add<std::string>("instance_file", '\0', "plain instance file", false, "");
    parser.add<std::string>("seed_file", '\0', "initial population file", false, "");
    parser.add<std::string>("result_file", '\0', "result file", false, "");
    parser.add<std::string>("lkh_executable", '\0', "LKH executable", false, "");
    parser.add<std::string>("lkh_temp_root", '\0', "LKH temporary directory", false, "");
    parser.add<std::string>("snapshot_dir", '\0', "population snapshot directory", false, "");
    // parameters
    parser.add<int>("pop_size", 'p', "population size", false, POPULATION_SIZE);
    parser.add<int>("iteration", 'r', "iteration", false, ITERATION);
    parser.add<int>("patience", '\0', "generations without improvement", false, 0);
    parser.add<int>("initial_patience", '\0', "initial generations without improvement", false, 0);
    parser.add<double>("max_time", 't', "max running time", false, MAX_TIME);
    parser.add<double>("fit_beta", 'b', "coefficient for fitness function", false, FIT_BETA);
    parser.add<int>("dist_th", 'd', "distance threshold", false, DISTANCE_THRESHOLD);
    parser.add<int>("neighbor_size", 'n', "neighbor size", false, NEIGHBOR_SIZE);
    parser.add<int>("solver_threads", '\0', "Gurobi threads per solve", false, 0);
    parser.add<int>("snapshot_interval", '\0', "population snapshot interval", false, 0);

    parser.parse_check(argc, argv);

    random_num = parser.get<int>("random");
    instance_index = parser.get<int>("instance");
    seed = parser.get<int>("seed");
    init = parser.get<std::string>("init");
    select = parser.get<std::string>("select");
    crossover = parser.get<std::string>("crossover");
    improvement = parser.get<std::string>("improvement");
    greed = parser.get<std::string>("greed");
    distance = parser.get<std::string>("distance");
    instance_file = parser.get<std::string>("instance_file");
    seed_file = parser.get<std::string>("seed_file");
    result_file = parser.get<std::string>("result_file");
    lkh_executable = parser.get<std::string>("lkh_executable");
    lkh_temp_root = parser.get<std::string>("lkh_temp_root");
    snapshot_dir = parser.get<std::string>("snapshot_dir");
    population_size = parser.get<int>("pop_size");
    iteration = parser.get<int>("iteration");
    patience = parser.get<int>("patience");
    if (patience == 0) {
        patience = iteration / 10;
    }
    initial_patience = parser.get<int>("initial_patience");
    if (initial_patience == 0) {
        initial_patience = patience;
    }
    max_time = parser.get<double>("max_time");
    fit_beta = parser.get<double>("fit_beta");
    dist_th = parser.get<int>("dist_th");
    neighbor_size = parser.get<int>("neighbor_size");
    solver_threads = parser.get<int>("solver_threads");
    snapshot_interval = parser.get<int>("snapshot_interval");
    if (patience <= 0 || initial_patience <= 0) {
        throw std::invalid_argument("patience values must be positive");
    }
    if (snapshot_interval < 0) {
        throw std::invalid_argument("snapshot_interval must be nonnegative");
    }
    if (solver_threads < 0) {
        throw std::invalid_argument("solver_threads must be nonnegative");
    }
    timestamp = std::to_string(std::time(nullptr));
}

void Parameters::print() const {
    std::cout << "[PARAMS] instance_index: " << instance_index
              << " random_num: " << random_num
              << " seed: " << seed
              << " population_size: " << population_size
              << " iteration: " << iteration
              << " patience: " << patience
              << " initial_patience: " << initial_patience
              << " max_time: " << max_time
              << " fit_beta: " << fit_beta
              << " dist_th: " << dist_th
              << " neighbor_size: " << neighbor_size
              << " solver_threads: " << solver_threads
              << " instance_file: " << instance_file
              << " seed_file: " << seed_file
              << " snapshot_dir: " << snapshot_dir
              << " snapshot_interval: " << snapshot_interval
              << " timestamp: " << timestamp
              << std::endl;
}
