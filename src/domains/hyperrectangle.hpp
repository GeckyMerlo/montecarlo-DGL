#ifndef MONTECARLO_1_HYPERRECTANGLE_HPP
#define MONTECARLO_1_HYPERRECTANGLE_HPP

#include "integration_domain.hpp"
#include <utility>
#include "../geometry.hpp"
using namespace std;
using namespace geom;

template <size_t dim>
class HyperRectangle : public IntegrationDomain<dim>{
public:
    HyperRectangle(array<double, dim> &dims);

    Bounds<dim> getBounds() const override;

    double getVolume() const override;

    bool isInside(const Point<dim> &point) const override;
private:
    array<double, dim> dimensions;
};

#include "hyperrectangle.tpp"

#endif //MONTECARLO_1_HYPERRECTANGLE_HPP