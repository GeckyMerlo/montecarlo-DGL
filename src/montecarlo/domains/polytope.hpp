//
// Created by Giacomo Merlo on 04/12/25.
//

#ifndef MONTECARLO_1_POLYTOPE_H
#define MONTECARLO_1_POLYTOPE_H

#include "integration_domain.hpp"
#include "../geometry.hpp"
#include <vector>

using namespace geom;
using namespace std;

template <size_t dim>
class PolyTope : public IntegrationDomain<dim> {
public:

    PolyTope(const std::vector<geom::Point<dim>>&   vertices,
                        const std::vector<array<double, dim>>&  norms,
                        const std::vector<double>&      offs);


    geom::Bounds<dim> getBounds() const override;


    double getBoxVolume() const override;


    bool isInside(const geom::Point<dim> &p) const override;

private:
    vector<Point<dim>> vec;
    vector<array<double, dim>> normals;
    vector<double> offsets;
};

#include "polytope.tpp"

#endif //MONTECARLO_1_POLYTOPE_H