//
// Created by Giacomo Merlo on 04/12/25.
//

#ifndef MONTECARLO_1_POLYTOPE_TPP
#define MONTECARLO_1_POLYTOPE_TPP

#include <cmath>
#include <algorithm> // For std::pow
#include "../geometry.hpp"
#include "integration_domain.hpp"
#include <vector>

using namespace geom;
using namespace std;

// Command to generate normals and offsets: qhull Qt Qx Fn < points.txt > hull.txt
// points.txt -- text file with first line: <num_points> <dim>, then coordinates of each point per line
template <size_t dim>
PolyTope<dim>::PolyTope(const std::vector<Point<dim>>&   vertices,
                        const std::vector<array<double, dim>>&  norms,
                        const std::vector<double>&      offs)
    : vec(vertices)
    , normals(norms)
    , offsets(offs)
{
    if (normals.size() != offsets.size()) {
        throw std::runtime_error(
            "PolyTope: normals and offsets must have the same size.");
    }

    if (vec.empty()) {
        throw std::runtime_error(
            "PolyTope: vertices list cannot be empty.");
    }
}

template<size_t dim>
Bounds<dim> PolyTope<dim>::getBounds() const{
    Bounds<dim> bounds;
    for (int i = 0; i < dim; ++i) {
        double max = vec[0][i];
        double min = vec[0][i];
        for (Point p: vec) {
            if (p[i] > max) max = p[i];
            if (p[i] < min) min = p[i];
        }
        bounds[i] = make_pair(min, max);
    }
    return bounds;
}

template<size_t dim>
double PolyTope<dim>::getBoxVolume() const {
    Bounds<dim> bou = this->getBounds();
    double vol = 1;
    for (int i = 0; i <dim; ++i) {
        vol *= (bou[i].second - bou[i].first);
    }
    return vol;
}

template<size_t dim>
bool PolyTope<dim>::isInside(const Point<dim> &point) const {
    const double tol = 1e-12; //tolleranza numerica
    for (size_t i = 0; i < normals.size(); ++i) {
        double s = 0.0;
        //Moltiplicando la normale per il punto costruisco l'iperpiano, di vettori perpendicolari
        //alla normale, passante per quel punto.
        for (size_t k = 0; k < dim; ++k)
            s += normals[i][k] * point[k];   // n_i · x

        if (s > offsets[i] + tol)
            return false;                //Se l'iperpiano viola la condizione
    }
    return true;
}




#endif //MONTECARLO_1_POLYTOPE_TPP