/**
 * Population.cpp
 * created on : March 03 2022
 * author : Z.LEI
 **/

#include "Genetic/Population.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

Population::Population(Parameters *params): neighbor(params->neighbor_size), ls(params) {
    this->initialization = params->init;
    this->selection = params->select;
    this->crossover_type = params->crossover;
    this->population_size = params->population_size;
    this->fit_beta = params->fit_beta;
    this->dist_th = params->dist_th;
    this->seed_file = params->seed_file;
    this->best_solution = nullptr;
    this->crossover = nullptr;
}

Population::~Population() {
    for (int i = 0; i < population.size(); ++i) {
        delete population[i];
    }
    delete best_solution;
    delete crossover;
}

void Population::setContext(Centers &centers, Random *random, std::string timestamp) {
    this->random = random;
    this->centers = centers;
    neighbor.setContext(centers.size());
    ls.setContext(random, centers, &neighbor, timestamp);
    kmeans.setContext(random, centers);
    crossover = CrossoverFactory::createCrossover(crossover_type);
    crossover->setContext(random);
}

List* Population::randomSolution() {
    std::vector<int> ids(centers.size() - 1);
    std::iota(ids.begin(), ids.end(), 1);
    random->permutation(ids);
    List *solution = new List();
    Node* head = new Node(0, centers[0][0], centers[0][1]);
    solution->add(head);
    for (int i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        double theta = random->randomDoubleDistr(0, 2 * PI);
        double r = random->randomDoubleDistr(0, 1);
        double x = centers[id][0] + r * centers[id][2] * cos(theta);
        double y = centers[id][1] + r * centers[id][2] * sin(theta);
        if (pow(x - centers[id][0], 2) + pow(y - centers[id][1], 2) > pow(centers[id][2], 2)){
            std::cout << "ERROR : init point out of circle" << std::endl;
        }
        Node* node = new Node(id, x, y);
        solution->add(node);
    }
    return solution;
}

List* Population::kmeansSolution() {
    std::vector<std::vector<int>> groups = kmeans.getGroups();
    // random.permutation(groups);
    List *solution = new List();
    for (int i = 0; i < groups.size(); ++i) {
        random->permutation(groups[i]);
        for (int j = 0; j < groups[i].size(); ++j) {
            int id = groups[i][j];
            double x = centers[id][0], y = centers[id][1];
            if (id != 0) {
                double theta = random->randomDoubleDistr(0, 2 * PI);
                double r = random->randomDoubleDistr(0, 1);
                x = centers[id][0] + r * centers[id][2] * cos(theta);
                y = centers[id][1] + r * centers[id][2] * sin(theta);
                if (pow(x - centers[id][0], 2) + pow(y - centers[id][1], 2) > pow(centers[id][2], 2)){
                    std::cout << "ERROR : init point out of circle" << std::endl;
                }
            }
            Node* node = new Node(id, x, y);
            solution->add(node);
        }
    }
    Node* p = solution->head();
    while (p->id != 0) {
        p = p->next;
    }
    solution->setHead(p);
    return solution;
}

List *Population::initSolution() {
    List* solution = nullptr;
    if (initialization == "RANDOM") {
        solution = randomSolution();
    } else if (initialization == "KMEANS") {
        solution = kmeansSolution();
    } else {
        std::cerr << "[ERROR] invalid initialization method" << std::endl;
    }
    solution = ls.initSolOpt(solution);
    return solution;
}

std::vector<List*> Population::readSeedSolutions() {
    std::vector<List*> solutions;
    if (seed_file.empty()) {
        return solutions;
    }

    std::ifstream input(seed_file);
    if (!input) {
        throw std::runtime_error("could not open seed file: " + seed_file);
    }
    std::string magic;
    int version = 0;
    std::size_t node_count = 0;
    std::size_t solution_count = 0;
    if (!(input >> magic >> version >> node_count >> solution_count) ||
        magic != "MA_CETSP_SEEDS" || version != 1 ||
        node_count != centers.size()) {
        throw std::runtime_error("invalid seed file header: " + seed_file);
    }

    solutions.reserve(solution_count);
    for (std::size_t solution_index = 0;
         solution_index < solution_count; ++solution_index) {
        std::string marker;
        if (!(input >> marker) || marker != "TOUR") {
            throw std::runtime_error("invalid seed tour marker");
        }
        auto* solution = new List();
        std::vector<bool> seen(node_count, false);
        for (std::size_t node_index = 0; node_index < node_count; ++node_index) {
            int id = -1;
            double x = 0.0;
            double y = 0.0;
            if (!(input >> id >> x >> y) || id < 0 ||
                static_cast<std::size_t>(id) >= node_count || seen[id]) {
                delete solution;
                throw std::runtime_error("invalid node in seed tour");
            }
            seen[id] = true;
            solution->add(new Node(id, x, y));
        }
        Node* head = solution->head();
        while (head->id != 0) {
            head = head->next;
        }
        solution->setHead(head);
        solution->evaluate();
        solutions.push_back(solution);
    }
    return solutions;
}

List *Population::initPopulation() {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<List*> seeded = readSeedSolutions();
    if (initialization == "KMEANS") {
        auto start_kmeans = std::chrono::high_resolution_clock::now();
        kmeans.run();
        auto end_kmeans = std::chrono::high_resolution_clock::now();
        std::cout << "kmeans time : "<< std::chrono::duration <double> (end_kmeans - start_kmeans).count() << " s"  << std::endl;
    }

    List* best = nullptr;
    auto add_initial = [&](List* solution) {
        solution = ls.initSolOpt(solution);
        if (!best || solution->getValue() < best->getValue()) {
            delete best;
            best = new List(*solution);
        }
        if (!insertSolution(solution)) {
            delete solution;
            return false;
        }
        return true;
    };

    std::size_t seeded_accepted = 0;
    for (List* solution : seeded) {
        if (population.size() >= static_cast<std::size_t>(population_size)) {
            delete solution;
            continue;
        }
        if (add_initial(solution)) {
            ++seeded_accepted;
        }
    }

    int attempts = 0;
    std::size_t generated_accepted = 0;
    const int max_attempts = std::max(10, population_size * 10);
    while (population.size() < static_cast<std::size_t>(population_size) &&
           attempts++ < max_attempts) {
        List* solution = nullptr;
        if (initialization == "RANDOM") {
            solution = randomSolution();
        } else if (initialization == "KMEANS") {
            solution = kmeansSolution();
        } else {
            throw std::runtime_error("invalid initialization method");
        }
        if (add_initial(solution)) {
            ++generated_accepted;
        }
    }

    std::cout << "[INIT] seed_tours: " << seeded.size()
              << " seeded_accepted: " << seeded_accepted
              << " generated_accepted: " << generated_accepted
              << " population: " << population.size() << std::endl;

    if (!best || population.size() < 2) {
        delete best;
        throw std::runtime_error("could not construct an initial population");
    }
    best_solution = best;

    std::sort(population.begin(), population.end(), [] (List* s1, List* s2) {
        return s1->getValue() < s2->getValue();
    });

    auto end = std::chrono::high_resolution_clock::now();
    if (LOG) {
        std::cout << initialization << " init time : "<< std::chrono::duration <double> (end - start).count() << " s"  << std::endl;
    }
    neighbor.updateNeighbors();
    return best_solution;
}

std::pair<List*, List*> Population::chooseParent() {
    int size = population.size();
    int i = -1, j = -1;
    if (selection == "RANDOM") {
        while (i == j) {
            i = random->randomInt(size);
            j = random->randomInt(size);
        }
    } else if (selection == "ROULETTE") {
        int sum_portion = (1 + size) * size / 2;
        int lucky;
        while (i == j) {
            lucky = random->randomInt(sum_portion) + 1;
            i = size - int(sqrt(2 * lucky));
            lucky = random->randomInt(sum_portion) + 1;
            j = size - int(sqrt(2 * lucky));
        }
        if (i > j) {
            i = i ^ j;
            j = i ^ j;
            i = i ^ j;
        }
    }
    if (LOG) {
        std::cout << "parents indices : " << i << " " << j << std::endl;
    }
    return {population[i], population[j]};
}

List *Population::nextPopulation(int patience) {
    std::pair<List*, List*> parents;
    double dist1 = 0, dist2 = 0;
    List* offspring = nullptr;
    int try_times = 5;
    while (dist1 == 0 || dist2 == 0) {
        if (try_times-- <= 0) {
            randomSwap(offspring);
        } else {
            parents = chooseParent();
            offspring = crossover->run(parents.first, parents.second);
        }
        dist1 = Distance::run(offspring, parents.first);
        dist2 = Distance::run(offspring, parents.second);
        if (LOG) {
            std::cout << "dists between offspring and parents :" << dist1 << " " << dist2 << std::endl;
        }
    }

    // mutation
    if (random->randomInt(1000) < patience) {
        randomSwap(offspring);
    }

    offspring = ls.VND(offspring);
    insertSolution(offspring);

    populationManagement();

    if (LOG) {
        std::cout << "standard population size : " << population_size << ", current population size : " << population.size() << std::endl;
        std::cout << "population distances : " << std::endl;
        for (auto &p : population) {
            std::cout << p->getDistance() << " ";
        }
        std::cout << std::endl;
        std::cout << "population costs : " << std::endl;
        for (auto &p : population) {
            std::cout << p->getValue() << " ";
        }
        std::cout << std::endl;
    }

    return best_solution;
}

bool Population::insertSolution(List *s) {
    double distance_threshold = dist_th;
    double min_dist = INT_MAX;
    std::vector<double> distances(population.size(), 0);
    for (int i = 0; i < population.size(); ++i) {
        double dist = Distance::run(s, population[i]);
        min_dist = std::min(min_dist, dist);
        distances[i] = std::min(population[i]->getDistance(), dist);
    }
    s->setDistance(min_dist);
    if ((min_dist > 0 && best_solution && s->getValue() < best_solution->getValue()) || min_dist > distance_threshold) {
        for (int i = 0; i < population.size(); ++i) {
            population[i]->setDistance(distances[i]);
        }
        population.emplace_back(s);
        neighbor.updateCentroids(s);
        return true;
    } else {
        return false;
    }
}

void Population::updateDistances() {
    std::vector<double> distances(population.size(), INT_MAX);
    for (int i = 0; i < population.size(); ++i) {
        for (int j = i + 1; j < population.size(); ++j) {
            double dist = Distance::run(population[i], population[j]);
            distances[i] = std::min(distances[i], dist);
            distances[j] = std::min(distances[j], dist);
        }
    }
    for (int i = 0; i < population.size(); ++i) {
        population[i]->setDistance(distances[i]);
    }
}

void Population::populationManagement() {
    std::unordered_map<List*, std::vector<double>> rank;

    // value rank
    std::sort(population.begin(), population.end(), [] (List* s1, List* s2) {
        return s1->getValue() < s2->getValue();
    });

    for (int i = 0; i < population.size(); ++i) {
        rank[population[i]].emplace_back(100.0 * i / (population.size() - 1));
    }

    // distance rank
    std::sort(population.begin(), population.end(), [] (List* s1, List* s2) {
        return s1->getDistance() > s2->getDistance();
    });

    for (int i = 0; i < population.size(); ++i) {
        rank[population[i]].emplace_back(100.0 * i / (population.size() - 1));
    }

    for (int i = 0; i < population.size(); ++i) {
        double alpha = 1, beta = fit_beta;
        double fitness = alpha * rank[population[i]][0] + beta * rank[population[i]][1];
        population[i]->setFitness(fitness);
    }

    std::sort(population.begin(), population.end(), [] (List* s1, List* s2) {
        return s1->getFitness() < s2->getFitness();
    });

    if (population.size() >= 1.5 * population_size) {
        population.resize(population_size);
        neighbor.updateNeighbors();
        updateDistances();
    }

    List* best = *std::min_element(population.begin(), population.end(), [] (List* s1, List* s2) {
        return s1->getValue() < s2->getValue();
    });

    if (best->getValue() < best_solution->getValue()) {
        delete best_solution;
        best_solution = new List(*best);
    }

}

void Population::randomSwap(List* s) {
    int steps = 5;
    while (steps-- > 0) {
        int i = random->randomInt(s->size()-1) + 1;
        int j = random->randomInt(s->size()-1) + 1;
        while (i == j && abs(i - j) == 1) {
            i = random->randomInt(s->size()-1) + 1;
            j = random->randomInt(s->size()-1) + 1;
        }
        Node* p1 = s->head();
        Node* p2 = s->head();
        while(i-- > 0) p1 = p1->next;
        while (j-- > 0) p2 = p2->next;
        Node::swap(p1, p2);
    }
}

void Population::writeSnapshot(
    const std::filesystem::path& directory,
    int generation,
    const std::string& reason) {
    if (directory.empty()) {
        return;
    }
    std::filesystem::create_directories(directory);

    std::ostringstream base_name;
    base_name << "generation-" << std::setw(6) << std::setfill('0')
              << generation << '-' << reason;
    const std::filesystem::path seed_path =
        directory / (base_name.str() + ".seeds");
    const std::filesystem::path member_path =
        directory / (base_name.str() + ".csv");

    std::vector<List*> ordered = population;
    std::sort(ordered.begin(), ordered.end(), [](List* left, List* right) {
        return left->getValue() < right->getValue();
    });

    std::ofstream seed_output(seed_path);
    if (!seed_output) {
        throw std::runtime_error("could not create population snapshot: " +
                                 seed_path.string());
    }
    seed_output << std::setprecision(std::numeric_limits<double>::max_digits10)
                << "MA_CETSP_SEEDS 1 " << centers.size() << ' '
                << ordered.size() << '\n';
    for (List* solution : ordered) {
        seed_output << "TOUR\n";
        Node* node = solution->head();
        for (int index = 0; index < solution->size(); ++index) {
            seed_output << node->id << ' ' << node->x << ' ' << node->y << '\n';
            node = node->next;
        }
    }

    std::vector<double> nearest(ordered.size(),
                                std::numeric_limits<double>::infinity());
    double pairwise_sum = 0.0;
    double pairwise_min = std::numeric_limits<double>::infinity();
    double pairwise_max = 0.0;
    std::size_t pairwise_count = 0;
    for (std::size_t first = 0; first < ordered.size(); ++first) {
        for (std::size_t second = first + 1; second < ordered.size(); ++second) {
            const double distance = Distance::run(ordered[first], ordered[second]);
            nearest[first] = std::min(nearest[first], distance);
            nearest[second] = std::min(nearest[second], distance);
            pairwise_sum += distance;
            pairwise_min = std::min(pairwise_min, distance);
            pairwise_max = std::max(pairwise_max, distance);
            ++pairwise_count;
        }
    }

    std::ofstream member_output(member_path);
    if (!member_output) {
        throw std::runtime_error("could not create population metadata: " +
                                 member_path.string());
    }
    member_output << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "rank,cost,nearest_edit_distance,fitness\n";

    std::vector<std::size_t> distance_order(ordered.size());
    std::iota(distance_order.begin(), distance_order.end(), 0);
    std::sort(distance_order.begin(), distance_order.end(),
              [&nearest](std::size_t left, std::size_t right) {
                  return nearest[left] > nearest[right];
              });
    std::vector<double> distance_rank(ordered.size(), 0.0);
    for (std::size_t rank = 0; rank < distance_order.size(); ++rank) {
        distance_rank[distance_order[rank]] = ordered.size() == 1
            ? 0.0
            : 100.0 * static_cast<double>(rank) /
                static_cast<double>(ordered.size() - 1);
    }
    for (std::size_t index = 0; index < ordered.size(); ++index) {
        const double value_rank = ordered.size() == 1
            ? 0.0
            : 100.0 * static_cast<double>(index) /
                static_cast<double>(ordered.size() - 1);
        const double fitness = value_rank + fit_beta * distance_rank[index];
        member_output << index + 1 << ',' << ordered[index]->getValue() << ','
                      << nearest[index] << ',' << fitness << '\n';
    }

    const std::filesystem::path manifest_path = directory / "manifest.csv";
    const bool reset_manifest = generation == 0 && reason == "initial";
    const bool write_header = reset_manifest ||
        !std::filesystem::exists(manifest_path) ||
        std::filesystem::file_size(manifest_path) == 0;
    const auto manifest_mode = reset_manifest
        ? std::ios::trunc
        : std::ios::app;
    std::ofstream manifest_output(manifest_path, manifest_mode);
    if (!manifest_output) {
        throw std::runtime_error("could not create population manifest: " +
                                 manifest_path.string());
    }
    if (write_header) {
        manifest_output
            << "generation,reason,population,best_cost,mean_cost,worst_cost,"
               "min_pairwise_edit_distance,mean_pairwise_edit_distance,"
               "max_pairwise_edit_distance,seed_file,member_file\n";
    }
    const double mean_cost = std::accumulate(
        ordered.begin(), ordered.end(), 0.0,
        [](double total, List* solution) {
            return total + solution->getValue();
        }) / static_cast<double>(ordered.size());
    const double mean_pairwise = pairwise_count == 0
        ? 0.0
        : pairwise_sum / static_cast<double>(pairwise_count);
    if (pairwise_count == 0) {
        pairwise_min = 0.0;
    }
    manifest_output
        << std::setprecision(std::numeric_limits<double>::max_digits10)
        << generation << ',' << reason << ',' << ordered.size() << ','
        << ordered.front()->getValue() << ',' << mean_cost << ','
        << ordered.back()->getValue() << ',' << pairwise_min << ','
        << mean_pairwise << ',' << pairwise_max << ','
        << seed_path.filename().string() << ','
        << member_path.filename().string() << '\n';
}
