//
// Created by Giacomo Merlo on 15/01/26.
//
#include "integrator.hpp"
#include "../geometry.hpp"
#include "../RngManager.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <omp.h>

template <size_t dim>
ISMontecarloIntegrator<dim>::ISMontecarloIntegrator(const IntegrationDomain<dim> &d)
    : Integrator<dim>(d) {}


template <size_t dim>
double ISMontecarloIntegrator<dim>::integrate(
    const function<double(const Point<dim>&)>& f,
    int n_samples,
    const Proposal<dim>& proposal,
    uint32_t seed)
{
    ISMeanEstimator<dim> mean_estimator;
    ImportanceEstimate<dim> mean_estimate = mean_estimator.estimate(this->domain, seed, n_samples, proposal, f);
    return mean_estimate.mean * this->domain.getBoxVolume();
}