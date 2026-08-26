/**
 * Solver.cpp
 * created on : Nov 30 2022
 * author : Z.LEI
 **/

#include "LocalSearch/Solver.hpp"

#include <cmath>
#include <limits>

Solver::Solver() {}
Solver::~Solver() {}

void Solver::setContext(Centers &centers) {
    this->centers = centers;
}

void Solver::setThreads(int value) {
    threads = value;
}

void Solver::solve(List* solution) {
    auto start = std::chrono::high_resolution_clock::now();
    // params
    const int n = centers.size();               // number of points
    std::vector<double> cx(n);                  // center x
    std::vector<double> cy(n);                  // center y
    std::vector<double> r(n);                   // radius

    Node *p = solution->head();
    for (int i = 0; i < n; ++i) {
        cx[i] = centers[p->id][0];
        cy[i] = centers[p->id][1];
        r[i] = centers[p->id][2];
        p = p->next;

    }

    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        // env.set("LogFile", "solve.log");
        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
        env.set(GRB_DoubleParam_BarQCPConvTol, 1e-10);
        if (threads > 0) {
            env.set(GRB_IntParam_Threads, threads);
        }
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        double* lbs = new double[n];
        double* ubs = new double[n];
        char* types = new char[n];
        for (int i = 0; i < n; ++i) {
            lbs[i] = -GRB_INFINITY;
            ubs[i] = GRB_INFINITY;
            types[i] = GRB_CONTINUOUS;
        }

        // Create variables
        GRBVar* dist = model.addVars(lbs, ubs, nullptr, types, nullptr, n);         // distance
        GRBVar* dx = model.addVars(lbs, ubs, nullptr, types, nullptr, n);           // x diff
        GRBVar* dy = model.addVars(lbs, ubs, nullptr, types, nullptr, n);           // y diff
        GRBVar* s = model.addVars(lbs, ubs, nullptr, types, nullptr, n);            // turing point x
        GRBVar* t = model.addVars(lbs, ubs, nullptr, types, nullptr, n);            // turing point y
        GRBVar* x = model.addVars(lbs, ubs, nullptr, types, nullptr, n);            // turing point x
        GRBVar* y = model.addVars(lbs, ubs, nullptr, types, nullptr, n);            // turing point y

        // Set objective
        GRBLinExpr obj = 0;
        for (int i = 0; i < n; ++i) {
            obj += dist[i];
        }
        model.setObjective(obj, GRB_MINIMIZE);

        for (int i = 0; i < n - 1; ++i) {
            model.addConstr(dx[i] == x[i] - x[i+1]);
            model.addConstr(dy[i] == y[i] - y[i+1]);
        }

        model.addConstr(dx[n-1] == x[n-1] - x[0]);
        model.addConstr(dy[n-1] == y[n-1] - y[0]);

        for (int i = 0; i < n; ++i) {
            model.addConstr(s[i] == cx[i] - x[i]);
            model.addConstr(t[i] == cy[i] - y[i]);
            model.addQConstr(dist[i] * dist[i] >= dx[i] * dx[i] + dy[i] * dy[i]);
            model.addQConstr(r[i] * r[i] >= s[i] * s[i] + t[i] * t[i]);
            model.addConstr(dist[i] >= 0);
        }

        // Optimize model
        model.optimize();
        const int status = model.get(GRB_IntAttr_Status);
        const int solution_count = model.get(GRB_IntAttr_SolCount);
        if (solution_count == 0) {
            if (LOG) {
                std::cout << "SOCP produced no solution; status = "
                          << status << std::endl;
            }
            return;
        }
        Node* p = solution->head();
        for (int i = 0 ; i < n; ++i) {
            p->x = x[i].get(GRB_DoubleAttr_X);
            p->y = y[i].get(GRB_DoubleAttr_X);
            const double offset_x = p->x - cx[i];
            const double offset_y = p->y - cy[i];
            const double center_distance = std::hypot(offset_x, offset_y);
            if (center_distance > r[i]) {
                const double interior_radius =
                    std::nextafter(r[i], 0.0);
                const double scale = center_distance == 0.0
                    ? 0.0
                    : interior_radius / center_distance;
                p->x = cx[i] + offset_x * scale;
                p->y = cy[i] + offset_y * scale;
            }
            p = p->next;
        }

        solution->evaluate();

        delete[] lbs;
        delete[] ubs;
        delete[] types;

    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    if (LOG) std::cout << "socp solution : " << solution->getValue() << " time : " << std::chrono::duration<double> (end - start).count() << " s" << std::endl;
}
