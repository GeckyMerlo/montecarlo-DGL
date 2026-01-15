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
#include <functional>

/**
 * @namespace montecarlo
 * @brief Core library namespace containing integrators, domains, samplers, and estimators.
 */

/**
 * @brief Primary Monte Carlo integration class supporting multiple sampling strategies
 * @tparam dim Dimensionality of the integration domain; must match the domain type
 * 
 * Provides three integration methods:
 * 1. Classic uniform Monte Carlo (integrate)
 * 2. Importance sampling with custom proposal (integrate_importance)
 * 3. Metropolis-Hastings MCMC sampling (integrate_with_mh)
 * 
 * All methods estimate ∫_D f(x) dx for a given domain D and integrand f.
 */
template <std::size_t dim>
class MontecarloIntegrator : public Integrator<dim> {
public:
    /**
     * @brief Construct integrator for a specific domain
     * @param d Integration domain (hypersphere, rectangle, cylinder, polytope, etc.)
     */
    explicit MontecarloIntegrator(const IntegrationDomain<dim> &d);

    /**
     * @brief Classic Monte Carlo integration with uniform sampling
     * @param f Integrand function to integrate
     * @param n_samples Number of Monte Carlo samples
     * @return Estimated integral value ∫_D f(x) dx
     * 
     * Uses Hit-or-Miss sampling within the bounding box. Error scales as O(1/√n).
     */
    // Calcola l'integrale di una funzione 'f' usando Monte Carlo
    double OLDintegrate(const std::function<double(const Point<dim>&)>& f, int n_samples);


    double integrate(const function<double(const Point<dim>&)>& f, int n_samples, const Proposal<dim>& proposal, uint32_t seed) override;
};

#include "MCintegrator.tpp"

#endif // MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP
