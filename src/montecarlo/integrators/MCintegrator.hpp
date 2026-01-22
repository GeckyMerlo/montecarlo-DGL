/**
 * @file MCintegrator.hpp
 * @brief Classic Monte Carlo integration engine.
 * @details Implements uniform sampling-based Monte Carlo integration for computing
 * definite integrals over N-dimensional domains with optional importance sampling.
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

namespace mc::integrators {

/**
 * @class MontecarloIntegrator
 * @brief Uniform-sampling Monte Carlo integrator for N-dimensional domains.
 * @tparam dim Dimensionality of the integration domain.
 * 
 * @details Computes the integral:
 * ∫_Ω f(x) dx ≈ V_Ω · (1/N) ∑ᵢ f(xᵢ)
 * 
 * where xᵢ are samples drawn uniformly from the bounding box and f(xᵢ) is zero
 * if xᵢ ∉ Ω (hit-or-miss style).
 * 
 * Modern implementation uses MCMeanEstimator for variance-reduced sampling.
 */
template <std::size_t dim>
class MontecarloIntegrator : public Integrator<dim> {
public:
    /**
     * @brief Construct a uniform Monte Carlo integrator for a domain.
     * @param d Reference to the integration domain.
     */
    explicit MontecarloIntegrator(const mc::domains::IntegrationDomain<dim>& d);

    /**
     * @brief Legacy integration routine (deprecated).
     * @deprecated Use integrate() instead with proposal and seed.
     * @param f Integrand function.
     * @param n_samples Number of sample points.
     * @return Estimated integral value.
     */
    double OLDintegrate(const std::function<double(const mc::geom::Point<dim>&)>& f, int n_samples);

    /**
     * @brief Compute the integral using uniform Monte Carlo sampling.
     * @param f Integrand function: ℝⁿ → ℝ.
     * @param n_samples Number of sample points to evaluate.
     * @param proposal Proposal distribution (ignored for uniform MC; for API consistency).
     * @param seed Random seed for reproducibility.
     * @return Estimated value of ∫_Ω f(x) dx.
     * 
     * @details Uses MCMeanEstimator to compute mean of f over the domain,
     * then multiplies by domain's bounding box volume.
     */
    double integrate(const std::function<double(const mc::geom::Point<dim>&)>& f,
                     int n_samples,
                     const mc::proposals::Proposal<dim>& proposal,
                     std::uint32_t seed) override;
};

} // namespace mc::integrators

#include "MCintegrator.tpp"

#endif // MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP