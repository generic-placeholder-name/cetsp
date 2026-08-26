/**
 * Parameters.hpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#ifndef CETSP_PARAMETERS_HPP
#define CETSP_PARAMETERS_HPP


#include <iostream>
#include <random>
#include <ctime>
#include "cmdline.h"
#include "Defs.hpp"

class Parameters {
public:
    int seed = 0;
    int random_num;
    std::string timestamp;
    std::string init;
    std::string select;
    std::string crossover;
    std::string improvement;
    std::string greed;
    std::string distance;
    std::string instance_file;
    std::string seed_file;
    std::string result_file;
    std::string lkh_executable;
    std::string lkh_temp_root;
    std::string snapshot_dir;
    int instance_index;
    int population_size;
    int iteration;
    int patience;
    int initial_patience;
    double max_time;
    double fit_beta;
    int dist_th;
    int neighbor_size;
    int solver_threads;
    int snapshot_interval;
    Parameters(int argc, char **argv);
    Parameters() = default;
    ~Parameters() = default;
    void print() const;
};


#endif //CETSP_PARAMETERS_HPP
