#include "gurobi_c++.h"
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>
#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Point;

// Solve: choose p_i = (x[i],y[i]) in each circle so that the
// sum of cyclic edge distances is minimized.
std::vector<Point> minimizeCycleLength(
    const std::vector<Point>& centers,
    const std::vector<double>& radii)
{
    int n = centers.size();
    if ((int)radii.size() != n)
        throw GRBException("centers/radii size mismatch", -1);

    // 1) Set up Gurobi environment & model
    GRBEnv   env(true);
    env.set("LogFile", "socp_boost.log");
    env.start();
    GRBModel model(env);

    // 2) Create variables: coordinates x[i], y[i]; slack t[i] ≥ 0
    std::vector<GRBVar> x(n), y(n), t(n);
    for (int i = 0; i < n; ++i) {
        x[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                            "x_" + std::to_string(i));
        y[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                            "y_" + std::to_string(i));
        t[i] = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS,
                            "t_" + std::to_string(i));
    }

    // 3) Circle constraints: (x_i - cx)^2 + (y_i - cy)^2 ≤ r_i^2
    for (int i = 0; i < n; ++i) {
        double cx = bg::get<0>(centers[i]);
        double cy = bg::get<1>(centers[i]);
        double  r = radii[i];

        GRBQuadExpr qc = (x[i] - cx) * (x[i] - cx)
                       + (y[i] - cy) * (y[i] - cy);
        model.addQConstr(qc <= r * r, "circle_" + std::to_string(i));
    }

    // 4) Edge‐length slacks: (x_i - x_j)^2 + (y_i - y_j)^2 ≤ t_i^2
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        GRBQuadExpr qc = (x[i] - x[j]) * (x[i] - x[j])
                       + (y[i] - y[j]) * (y[i] - y[j]);
        model.addQConstr(qc <= t[i] * t[i], "edge_" + std::to_string(i));
    }

    // 5) Objective: minimize ∑ t_i
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // 6) Optimize
    model.optimize();
    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        throw GRBException("No optimal solution", model.get(GRB_IntAttr_Status));
    }

    // 7) Extract solution into Boost Points
    std::vector<Point> sol(n);
    for (int i = 0; i < n; ++i) {
        double xi = x[i].get(GRB_DoubleAttr_X);
        double yi = y[i].get(GRB_DoubleAttr_X);
        sol[i] = Point(xi, yi);
    }
    return sol;
}

int main() {
    // Example centers and radii
    std::vector<Point> centers = {
        Point(0.0, 0.0),
        Point(1.0, 0.0),
        Point(0.5, 0.8)
    };
    std::vector<double> radii = {0.2, 0.3, 0.25};

    try {
        auto sol = minimizeCycleLength(centers, radii);
        std::cout << "Optimal points:\n";
        for (size_t i = 0; i < sol.size(); ++i) {
            std::cout << " p" << i << " = ("
                      << bg::get<0>(sol[i]) << ", "
                      << bg::get<1>(sol[i]) << ")\n";
        }
    } catch (GRBException& e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
        return 1;
    }
    return 0;
}
