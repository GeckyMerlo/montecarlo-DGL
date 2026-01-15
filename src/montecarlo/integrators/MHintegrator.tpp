//
// Created by Giacomo Merlo on 15/01/26.
//

#include "../RngManager.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <omp.h>

template <std::size_t dim>
MHMontecarloIntegrator<dim>::MHMontecarloIntegrator(const IntegrationDomain<dim> &d)
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
                                             const Proposal<dim>&,
                                             std::uint32_t seed)
{
    if (!configured)
        throw std::runtime_error("MHMontecarloIntegrator not configured. Call setConfig(...) first.");
    if (n_samples <= 0)
        throw std::invalid_argument("n_samples must be > 0");

    VolumeEstimatorMC<dim> ve;
    const auto vol_hat = ve.estimate(this->domain, seed, n_samples_volume);
    std::cout << "Volume: " << vol_hat.volume << " +- " << 2 * vol_hat.stderr << "\n";

    long double sum = 0.0;
    std::size_t kept = 0;

    Point x{};

	const int T = omp_get_max_threads();

    const size_t base = n_samples / T;
    const size_t rem  = n_samples % T;

	RngManager rngs(seed);

	#pragma omp parallel reduction(+:sum, kept)
	{
		const int tid = omp_get_thread_num();
		const int n_local = base + (tid < rem ? 1 : 0);

		auto rng = rngs.make_rng(tid);
	    MetropolisHastingsSampler<dim> mh_local(this->domain, p, x0, deviation);

		for (std::size_t i = 0; i < burn_in; ++i)
        	(void)mh_local.next(rng);

    	for (int i = 0; i < n_local; ++i) {
        	for (std::size_t t = 0; t < thinning; ++t)
            	x = mh_local.next(rng);

        	const double fx = f(x);
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