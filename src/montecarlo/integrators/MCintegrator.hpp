// MCintegrator.hpp
/**
 * @file MCintegrator.hpp
 * @brief Monte Carlo integration engine with multiple sampling strategies
 *
 * Implements classic Monte Carlo, importance sampling, and Metropolis-Hastings
 * integration methods for computing definite integrals over complex domains.
 */

#ifndef MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP
#define MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP

#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "../proposals/uniformProposal.hpp"
#include "../mcmc/metropolisHastingsSampler.hpp"
#include "../estimators/VolumeEstimatorMC.hpp"
#include "../estimators/ISMeanEstimator.hpp"
#include "../estimators/MCMeanEstimator.hpp"
#include "../geometry.hpp"
#include "integrator.hpp"

#include <cstdint>
#include <functional>

template <std::size_t dim>
class MontecarloIntegrator : public Integrator<dim> {
public:
    explicit MontecarloIntegrator(const IntegrationDomain<dim>& d);

    double OLDintegrate(const std::function<double(const Point<dim>&)>& f, int n_samples);

    double integrate(const std::function<double(const Point<dim>&)>& f,
                     int n_samples,
                     const Proposal<dim>& proposal,
                     std::uint32_t seed) override;
};

#include "MCintegrator.tpp"

#endif // MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP