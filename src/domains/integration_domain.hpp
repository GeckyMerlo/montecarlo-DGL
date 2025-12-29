#ifndef MONTECARLO_1_INTEGRATION_DOMAIN_HPP
#define MONTECARLO_1_INTEGRATION_DOMAIN_HPP

#include <array>
#include "../geometry.hpp"

using namespace std;
using namespace geom;

template <size_t dim>
class IntegrationDomain {
public:

    // returns the boundaries o f each dimension of the integration domain
    virtual Bounds<dim> getBounds() const = 0;

    // returns the volume of the integration domain
    virtual double getVolume() const = 0;

    // returns true if the given point is inside the integration domain
    virtual bool isInside(const Point<dim> &point) const = 0;

    virtual ~IntegrationDomain() = default;

};

#endif //MONTECARLO_1_INTEGRATION_DOMAIN_HPP