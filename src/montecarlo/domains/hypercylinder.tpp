#include <cmath>
#include <algorithm> // For std::pow
#include "../geometry.hpp"
#include "integration_domain.hpp"

using namespace geom;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <size_t dim>
HyperCylinder<dim>::HyperCylinder(double rad, double h)
    : radius(rad), height(h)
{
    static_assert(dim >= 2, "HyperCylinder requires at least 2 dimensions");
}

template <size_t dim>
auto HyperCylinder<dim>::getBounds() const -> Bounds<dim> {
    Bounds<dim> bounds;

    // Set bounds for the hypersphere base dimensions (0 to n-2)
    for (size_t i = 0; i < dim - 1; ++i) {
        bounds[i] = make_pair(-radius, radius);
    }

    // Set bounds for the height dimension (last dimension, n-1)
    bounds[dim - 1] = make_pair(0.0, height);

    return bounds;
}

template <size_t dim>
double HyperCylinder<dim>::getBoxVolume() const {
    // Volume of the bounding hypercube: (2*r)^(dim-1) * height
    return std::pow(2.0 * radius, static_cast<double>(dim - 1)) * height;
}

template <size_t dim>
bool HyperCylinder<dim>::isInside(const Point<dim> &point) const {
    // 1. Check the height constraint (last dimension)
    double h_val = point[dim - 1];
    if (h_val < 0.0 || h_val > height) {
        return false;
    }

    // 2. Check the radial constraint (hypersphere base)
    // Sum of squares of the first (n-1) coordinates
    double radial_dist_squared = 0.0;
    for (size_t i = 0; i < dim - 1; ++i) {
        double coord = point[i];
        radial_dist_squared += coord * coord;
    }

    return radial_dist_squared <= (radius * radius);
}