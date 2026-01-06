#ifndef MONTECARLO_1_HYPERCYLINDER_HPP
#define MONTECARLO_1_HYPERCYLINDER_HPP

#include "integration_domain.hpp"
#include "../geometry.hpp"
#include <cmath>

/**
 * @brief N-dimensional HyperCylinder implementation.
 * Defined as an (n-1)-dimensional hypersphere of radius 'r'
 * extended along the n-th dimension by height 'h'.
 */
template <size_t dim>
class HyperCylinder : public IntegrationDomain<dim> {
public:
    /**
     * @param rad Radius of the hypersphere base.
     * @param h   Height along the last dimension.
     */
    HyperCylinder(double rad, double h);

    /**
     * @brief Returns the bounding hypercube.
     * Bounds are [-r, r] for the first dim-1 dimensions,
     * and [0, h] for the last dimension.
     */
    geom::Bounds<dim> getBounds() const override;

    /**
     * @brief Returns the volume of the bounding box (Sampling Volume).
     * Used for Monte Carlo weight calculation: (2r)^(n-1) * h.
     */
    double getBoxVolume() const override;

    /**
     * @brief Checks if point is inside the hypercylinder.
     */
    bool isInside(const geom::Point<dim> &p) const override;

private:
    double radius;
    double height;
};

// Ensure the implementation code is saved as "hypercylinder.tpp"
#include "hypercylinder.tpp"

#endif //MONTECARLO_1_HYPERCYLINDER_HPP