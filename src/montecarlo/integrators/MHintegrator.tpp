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

namespace mc::integrators {

template <std::size_t dim>
MHMontecarloIntegrator<dim>::MHMontecarloIntegrator(const mc::domains::IntegrationDomain<dim>& d)
    : Integrator<dim>(d)
{}

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

    if (thinning == 0) throw std::invalid_argument("thinning must be > 0");
    if (deviation <= 0.0) throw std::invalid_argument("deviation must be > 0");
    if (n_samples_volume == 0) throw std::invalid_argument("n_samples_volume must be > 0");
    if (!this->domain.isInside(x0)) throw std::invalid_argument("x0 must be inside the domain");
}

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

    mc::estimators::VolumeEstimatorMC<dim> ve;
    const auto vol_hat = ve.estimate(this->domain, seed, n_samples_volume);

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

        for (std::size_t i = 0; i < burn_in; ++i)
            (void)mh_local.next(rng);

        for (std::size_t i = 0; i < n_local; ++i) {
            for (std::size_t t = 0; t < thinning; ++t)
                x_local = mh_local.next(rng);

            const double fx = f(x_local);
            if (std::isfinite(fx)) {
                sum += static_cast<long double>(fx);
                ++kept;
            }
        }
    }

    if (kept == 0)
        throw std::runtime_error("All sampled f(x) were non-finite.");

    const double mean_f = static_cast<double>(sum / static_cast<long double>(kept));
    return vol_hat.volume * mean_f;
}

} // namespace mc::integrators
