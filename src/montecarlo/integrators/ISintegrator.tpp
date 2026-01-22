/**
 * @file ISintegrator.tpp
 * @brief ISMontecarloIntegrator template implementation.
 * @details Contains inline implementations for importance sampling integration.
 */
//
// Created by Giacomo Merlo on 15/01/26.
//
#include "integrator.hpp"
#include "../geometry.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <omp.h>

namespace mc::integrators {


template <size_t dim>
ISMontecarloIntegrator<dim>::ISMontecarloIntegrator(const mc::domains::IntegrationDomain<dim> &d)
    : Integrator<dim>(d) {}


template <size_t dim>
double ISMontecarloIntegrator<dim>::integrate(
    const std::function<double(const mc::geom::Point<dim>&)>& f,
    int n_samples,
    const mc::proposals::Proposal<dim>& proposal,
    std::uint32_t seed)
{
    mc::estimators::ISMeanEstimator<dim> mean_estimator;
    mc::estimators::ImportanceEstimate<dim> mean_estimate = mean_estimator.estimate(this->domain, seed, n_samples, proposal, f);
    return mean_estimate.mean;
}

} // namespace mc::integrators
