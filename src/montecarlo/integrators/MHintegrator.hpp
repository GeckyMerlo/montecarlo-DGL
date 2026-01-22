/**
 * @file MHintegrator.hpp
 * @brief Metropolis-Hastings MCMC Monte Carlo integration engine.
 * @details Implements Markov Chain Monte Carlo integration for complex domains
 * and multimodal integrands, combining volume estimation with MCMC sampling.
 */

#ifndef MONTECARLO_DGL_MHINTEGRATOR_HPP
#define MONTECARLO_DGL_MHINTEGRATOR_HPP

#include "integrator.hpp"
#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "../mcmc/metropolisHastingsSampler.hpp"
#include "../estimators/VolumeEstimatorMC.hpp"
#include "../geometry.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>

namespace mc::integrators {

/**
 * @class MHMontecarloIntegrator
 * @brief Metropolis-Hastings MCMC integrator for complex domains.
 * @tparam dim Dimensionality of the integration domain.
 * 
 * @details Two-stage integration approach:
 * 1. Volume Estimation: Hit-or-Miss MC to estimate V_Ω
 * 2. Mean Estimation: MH sampler with target density p(x) to estimate E[f]
 * 
 * Integral: ∫_Ω f(x) dx ≈ V̂_Ω · (1/M) ∑ᵢ f(xᵢ) where xᵢ ~ π(x) ∝ p(x)
 * 
 * This approach is especially effective for:
 * - Complex domains with intricate boundaries
 * - Multimodal or highly non-uniform integrands
 * - Cases where uniform sampling has poor convergence
 * 
 * @note Configuration must be set via setConfig() before integrate().
 */
template <std::size_t dim>
class MHMontecarloIntegrator : public Integrator<dim> {
public:
    using Point = mc::geom::Point<dim>;
    using Func  = std::function<double(const Point&)>;

    /**
     * @brief Construct an MH-based integrator for a domain.
     * @param d Reference to the integration domain.
     */
    explicit MHMontecarloIntegrator(const mc::domains::IntegrationDomain<dim>& d);

    /**
     * @brief Compute the integral using MH sampling combined with volume estimation.
     * @param f Integrand function: ℝⁿ → ℝ.
     * @param n_samples Number of post-burn-in, post-thinning samples to keep.
     * @param proposal Unused (MH uses internal sampler).
     * @param seed Random seed for reproducibility.
     * @return Estimated integral ∫_Ω f(x) dx.
     * @throws std::runtime_error if setConfig() not called first.
     * @throws std::invalid_argument if n_samples ≤ 0.
     * @throws std::runtime_error if all samples yield non-finite f(x).
     */
    double integrate(const Func& f,
                     int n_samples,
                     const mc::proposals::Proposal<dim>& proposal,
                     std::uint32_t seed) override;

    /**
     * @brief Configure the MCMC sampler parameters.
     * @param burn_in_ Number of initial samples to discard (warmup).
     * @param thinning_ Keep every k-th sample (k=thinning_) for autocorrelation reduction.
     * @param n_samples_volume_ Number of samples for volume estimation.
     * @param deviation_ Proposal standard deviation (controls MH acceptance rate).
     * @param p_ Target density: π(x) ∝ p(x); should be proportional to integrand.
     * @param x0_ Initial point in the domain (burn-in starting point).
     * @throws std::invalid_argument if parameters are invalid (e.g., thinning=0, x0 ∉ Ω).
     * 
     * @note Must be called before integrate().
     */
    void setConfig(std::size_t burn_in_,
                   std::size_t thinning_,
                   std::size_t n_samples_volume_,
                   double deviation_,
                   Func p_,
                   Point x0_);

private:
    /** Configuration flag; true after setConfig() is called. */
    bool configured = false;

    /** Number of burn-in iterations (discarded). */
    std::size_t burn_in = 0;
    /** Thinning factor: keep every k-th sample. */
    std::size_t thinning = 1;
    /** Samples for volume estimation. */
    std::size_t n_samples_volume = 0;
    /** MH proposal deviation (Gaussian proposal standard deviation). */
    double deviation = 1.0;

    /** Target density proportional to integrand. */
    std::function<double(const Point&)> p;
    /** Initial MCMC state (must be in domain). */
    geom::Point<dim> x0;
};
 
} // namespace mc::integrators

#include "MHintegrator.tpp"

#endif // MONTECARLO_DGL_MHINTEGRATOR_HPP