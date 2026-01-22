/**
 * @file MCintegrator.tpp
 * @brief MontecarloIntegrator template implementation.
 * @details Contains inline implementations for classic uniform Monte Carlo
 * and importance-weighted integration.
 */

// MCintegrator.tpp
#include "integrator.hpp"
#include "../geometry.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace mc::integrators {

template <size_t dim>
MontecarloIntegrator<dim>::MontecarloIntegrator(const mc::domains::IntegrationDomain<dim>& d)
    : Integrator<dim>(d) {}

template <size_t dim>
double MontecarloIntegrator<dim>::OLDintegrate(const std::function<double(const mc::geom::Point<dim>&)>& f,
                                               int n_samples)
{
    std::vector<mc::geom::Point<dim>> points = this->initializeRandomizer(n_samples);

    double sum = 0.0;
    for (const auto& p : points) {
        if (this->domain.isInside(p)) sum += f(p);
    }

    const double volume = this->domain.getBoxVolume();
    return (sum / static_cast<double>(n_samples)) * volume;
}

template <size_t dim>
double MontecarloIntegrator<dim>::integrate(const std::function<double(const mc::geom::Point<dim>&)>& f,
                                            int n_samples,
                                            const mc::proposals::Proposal<dim>&,
                                            std::uint32_t seed)
{
    mc::estimators::MCMeanEstimator<dim> mean_estimator;
    mc::estimators::MeanEstimate<dim> mean_estimate =
        mean_estimator.estimate(this->domain, seed, static_cast<std::size_t>(n_samples), f);

    return mean_estimate.mean * this->domain.getBoxVolume();
}

} // namespace mc::integrators