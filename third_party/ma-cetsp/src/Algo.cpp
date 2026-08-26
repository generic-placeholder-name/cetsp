/**
 * Algo.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "Algo.hpp"

Algo::Algo(Parameters *params): population(params), data(params) {
    this->params = params;
    this->random = new Random(params->seed);
    this->iteration = params->iteration;
    this->patience_threshold = params->patience;
    this->initial_patience_threshold = params->initial_patience;
    this->instance_index = params->instance_index;
    this->timestamp = params->timestamp;
}

Algo::~Algo(){
    delete random;
}

void Algo::run() {
    // read data;
    Centers centers = data.getData();
    // set context
    population.setContext(centers, random, timestamp);

    auto start_run = std::chrono::high_resolution_clock::now();

    int iter = 0;
    int best_iter = 0;
    int patience = 0;
    bool improved = false;
    bool escaped_initial_incumbent = false;
    std::chrono::duration<double> best_running_time{};
    // init population
    List *best_solution = population.initPopulation();
    const auto end_initialization =
        std::chrono::high_resolution_clock::now();
    const double initialization_time =
        std::chrono::duration<double>(end_initialization - start_run).count();
    best_running_time = end_initialization - start_run;
    double best_solution_value = best_solution->getValue();
    if (!params->snapshot_dir.empty()) {
        population.writeSnapshot(params->snapshot_dir, 0, "initial");
    }

    // iteration
    int last_iter = 0;
    while (iter++ < iteration) {
        last_iter = iter;
        std::cout << std::endl << "Iteration " << iter << " : " << std::endl;

        auto start_iter = std::chrono::high_resolution_clock::now();
        // next population
        best_solution = population.nextPopulation(patience);
        improved = best_solution->getValue() - best_solution_value < -DELTA ? true : false;
        best_solution_value = best_solution->getValue();

        auto end_iter = std::chrono::high_resolution_clock::now();

        if (improved) {
            best_iter = iter;
            patience = 0;
            escaped_initial_incumbent = true;
            best_running_time = end_iter - start_run;
            if (LOG) data.write(best_solution, iter, std::to_string(best_running_time.count()));
        } else {
            ++patience;
        }

        if (!params->snapshot_dir.empty()) {
            if (improved) {
                population.writeSnapshot(params->snapshot_dir, iter, "best");
            } else if (params->snapshot_interval > 0 &&
                       iter % params->snapshot_interval == 0) {
                population.writeSnapshot(params->snapshot_dir, iter, "periodic");
            }
        }

        std::cout << "[LOG] iter: " << iter
                  << " best_iter: " << best_iter
                  << " best_value: " << best_solution->getValue()
                  << " iter_time: " << std::chrono::duration<double> (end_iter - start_iter).count()
                  << " total_time: " << std::chrono::duration<double> (end_iter - start_run).count()
                  << std::endl;

        const int active_patience_threshold = escaped_initial_incumbent
            ? patience_threshold
            : initial_patience_threshold;
        if (patience >= active_patience_threshold) {
            std::cout << std::endl << "[STOP] best solution hasn't been improved since " << active_patience_threshold << " iterations" << std::endl;
            break;
        }

        if (std::chrono::duration<double> (end_iter - start_run).count() > params->max_time) {
            std::cout << std::endl << "[STOP] running time exceeds " << params->max_time << " seconds" << std::endl;
            break;
        }
    }

    auto end_run = std::chrono::high_resolution_clock::now();

    if (!params->snapshot_dir.empty()) {
        population.writeSnapshot(params->snapshot_dir, last_iter, "final");
    }

    if (!LOG)  data.write(best_solution, best_iter, std::to_string(best_running_time.count()));

    std::cout << std::endl
            << "[SUMMARY] instance: " << data.getName()
            << " best_value: " << best_solution->getValue()
            << " best_time: " << std::to_string(best_running_time.count())
            << " initialization_time: " << initialization_time
            << " total_time: " << std::chrono::duration<double> (end_run - start_run).count()
            << " result_file: " << data.getResultFilename()
            << std::endl;
}
