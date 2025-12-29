#ifndef HYPERSPHERE_HPP
#define HYPERSPHERE_HPP

#include "integration_domain.hpp"
#include <utility>
#include "../geometry.hpp"

using namespace std;
using namespace geom;

template <size_t dim>

class Hypersphere : public IntegrationDomain<dim>{
public:
    Hypersphere(double rad);

    Bounds<dim> getBounds() const override;

    double getVolume() const override;

    bool isInside(const Point<dim> &point) const override;
private:
    double radius;
};

#include "hypersphere.tpp"

#endif //HYPERSPHERE_HPP