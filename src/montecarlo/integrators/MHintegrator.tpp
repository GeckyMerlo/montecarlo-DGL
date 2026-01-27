/**
 * @file MHintegrator.tpp
 * @brief MHMontecarloIntegrator template implementation.
 * @details Contains inline implementations for Metropolis-Hastings MCMC integration,
 * including parallel volume estimation, burn-in, thinning, and mean computation.
 */
// MHintegrator.tpp
//
// Created by Giacomo Merlo on 15/01/26.
//

#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>

#include <omp.h>

#include "montecarlo/rng/rng_factory.hpp"
#include "montecarlo/integrators/MCintegrator.hpp"
#include "montecarlo/proposals/uniformProposal.hpp"

namespace mc::integrators {

/**
 * @brief Construct an MH-based integrator for a domain.
 * @tparam dim Dimensionality parameter.
 * @param d Reference to the integration domain.
 */
template <std::size_t dim>
MHMontecarloIntegrator<dim>::MHMontecarloIntegrator(const mc::domains::IntegrationDomain<dim>& d)
    : Integrator<dim>(d)
{}

/**
 * @brief Configure the MCMC sampler and volume estimation parameters.
 * @tparam dim Dimensionality parameter.
 * @param burn_in_ Number of initial MH samples to discard (warmup iterations).
 * @param thinning_ Thinning factor: keep every k-th sample (k = thinning_).
 * @param n_samples_volume_ Number of samples for Hit-or-Miss volume estimation.
 * @param deviation_ Proposal standard deviation (Gaussian random walk step size).
 * @param p_ Target density proportional to integrand: π(x) ∝ p(x).
 * @param x0_ Initial MH state (must be inside domain).
 * 
 * @throws std::invalid_argument If parameters are invalid:
 *         - thinning_ == 0
 *         - deviation_ ≤ 0
 *         - n_samples_volume_ == 0
 *         - x0_ ∉ domain
 * 
 * @details Must be called before integrate(). Validates all parameters.
 */
template <std::size_t dim>
void MHMontecarloIntegrator<dim>::setConfig(std::size_t burn_in_,
                                           std::size_t thinning_,
                                           std::size_t n_samples_volume_,
                                           double deviation_,
                                           Func p_,
                                           Point x0_)
{
    burn_in = burn_in_;
    thinning = thinning_;
    n_samples_volume = n_samples_volume_;
    deviation = deviation_;
    p = std::move(p_);
    x0 = x0_;
    configured = true;
    const double p0 = p(x0);

    if (thinning == 0) throw std::invalid_argument("thinning must be > 0");
    if (deviation <= 0.0) throw std::invalid_argument("deviation must be > 0");
    if (n_samples_volume == 0) throw std::invalid_argument("n_samples_volume must be > 0");
    if (!this->domain.isInside(x0)) throw std::invalid_argument("x0 must be inside the domain");
    if (!(p0 > 0.0) || !std::isfinite(p0)) throw std::invalid_argument("p(x0) must be finite and > 0 (choose x0 in support of p).");
}

/**
 * @brief Compute the integral using MH sampling combined with volume estimation.
 * @tparam dim Dimensionality parameter.
 * @param f Integrand function: ℝⁿ → ℝ.
 * @param n_samples Number of (post-burn-in, post-thinning) samples to keep.
 * @param proposal Unused (MH uses internal sampler with configured p_).
 * @param seed Random seed for reproducibility (used for both volume and MH sampling).
 * @return Estimated integral ∫_Ω f(x) dx.
 * 
 * @throws std::runtime_error If setConfig() was not called first.
 * @throws std::invalid_argument If n_samples ≤ 0.
 * @throws std::runtime_error If all sampled f(x) are non-finite.
 * 
 * @details Two-stage algorithm:
 * 1. **Volume Estimation**: Hit-or-Miss MC with n_samples_volume samples
 *    estimates V̂_Ω
 * 2. **Mean Estimation**: MH sampling with:
 *    - Burn-in: discard first burn_in iterations (warmup)
 *    - Thinning: keep every thinning-th sample (autocorrelation reduction)
 *    - Parallel: each thread runs independent MH chain
 * 3. **Integration**: ∫_Ω f dx ≈ V̂_Ω · (1/M) ∑ f(xᵢ)
 *    where M = number of kept samples
 */
template <std::size_t dim>
double MHMontecarloIntegrator<dim>::integrate(const Func& f,
                                             int n_samples,
                                             const mc::proposals::Proposal<dim>&,
                                             std::uint32_t seed)
{
    if (!configured)
        throw std::runtime_error("MHMontecarloIntegrator not configured. Call setConfig(...) first.");
    if (n_samples <= 0)
        throw std::invalid_argument("n_samples must be > 0");

    //Estimate Zp = ∫_Ω p(x) dx with a plain MC integrator (uniform samplI)
    mc::proposals::UniformProposal<dim> proposal(this->domain);
    mc::integrators::MontecarloIntegrator<dim> uni(this->domain);
    const double Zp_hat = uni.integrate(p, n_samples_volume, proposal, seed);

    long double sum = 0.0L;
    std::size_t kept = 0;

    const int T = omp_get_max_threads();
    const std::size_t base = static_cast<std::size_t>(n_samples) / static_cast<std::size_t>(T);
    const std::size_t rem  = static_cast<std::size_t>(n_samples) % static_cast<std::size_t>(T);

#pragma omp parallel reduction(+:sum, kept)
    {
        const int tid = omp_get_thread_num();
        const std::size_t n_local = base + (static_cast<std::size_t>(tid) < rem ? 1u : 0u);

        auto rng = mc::rng::make_engine_with_seed(std::optional<std::uint32_t>{seed},
                                                  static_cast<std::uint64_t>(tid));

        mc::mcmc::MetropolisHastingsSampler<dim> mh_local(this->domain, p, x0, deviation);
        Point x_local{};

        // Burn-in
        for (std::size_t i = 0; i < burn_in; ++i)
            (void)mh_local.next(rng);

        // Sampling + thinning
        for (std::size_t i = 0; i < n_local; ++i) {
            for (std::size_t t = 0; t < thinning; ++t)
                x_local = mh_local.next(rng);

            const double px = p(x_local);
            if (!(px > 0.0) || !std::isfinite(px))
                continue;

            const double fx = f(x_local);
            if (std::isfinite(fx)) {
                sum += static_cast<long double>(fx / px);
                ++kept;
            }
        }
    }

    if (kept == 0)
        throw std::runtime_error("All sampled f(x)/p(x) were invalid (p<=0 or non-finite).");

    const double mean_f_over_p = static_cast<double>(sum / static_cast<long double>(kept));

    //Integral estimate
    return Zp_hat * mean_f_over_p;
}

} //namespace mc::integrators
