/**
 * @file ISMeanEstimator.tpp
 * @brief ISMeanEstimator template implementation.
 * @details Implements importance sampling mean estimation with parallel sampling,
 * proposal weighting, and standard error computation.
 */

#ifndef MONTECARLO_DGL_ISMEANESTIMATOR_TPP
#define MONTECARLO_DGL_ISMEANESTIMATOR_TPP

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

#include <omp.h>

#include "montecarlo/rng/rng_factory.hpp"

namespace mc::estimators {

/**
 * @brief Estimates the mean of a function over a domain using importance sampling.
 * @tparam dim The dimensionality of the integration domain.
 * 
 * @param domain The integration domain constraint.
 * @param seed Random seed for reproducibility.
 * @param n_samples Number of samples to generate.
 * @param proposal The proposal distribution for importance sampling.
 * @param f The function to estimate the mean of.
 * 
 * @return ImportanceEstimate<dim> containing:
 *   - n_samples: Total samples generated
 *   - n_inside: Samples that fell within the domain
 *   - mean: Estimated mean value (sum of weighted samples / n_samples)
 *   - stderr: Standard error of the estimate
 * 
 * @throws std::invalid_argument If n_samples is zero.
 * 
 * @details
 * The estimator uses importance sampling with a proposal distribution q(x).
 * Each sample p from the proposal is weighted by f(p)/q(p). The algorithm:
 * 1. Samples from the proposal distribution
 * 2. Rejects samples outside the domain
 * 3. Computes weighted contribution f(p)/q(p)
 * 4. Accumulates first and second moments
 * 5. Computes sample variance for standard error
 * 
 * The computation is parallelized across available threads with each thread
 * maintaining independent accumulators to avoid synchronization overhead.
 * 
 * @note Uses OpenMP for parallel sampling.
 */
template <std::size_t dim>
ImportanceEstimate<dim> ISMeanEstimator<dim>::estimate(const mc::domains::IntegrationDomain<dim>& domain,
                                                      std::uint32_t seed,
                                                      std::size_t n_samples,
                                                      const mc::proposals::Proposal<dim>& proposal,
                                                      const std::function<double(const mc::geom::Point<dim>&)>& f) const
{
    if (n_samples == 0) throw std::invalid_argument("n_samples must be > 0");

    const int T = omp_get_max_threads();
    const std::size_t base = n_samples / static_cast<std::size_t>(T);
    const std::size_t rem  = n_samples % static_cast<std::size_t>(T);

    double sum = 0.0;
    double sum2 = 0.0;
    std::size_t inside_total = 0;

#pragma omp parallel for reduction(+:sum,sum2,inside_total)
    for (int tid = 0; tid < T; ++tid) {
        const std::size_t n_local = base + (static_cast<std::size_t>(tid) < rem ? 1u : 0u);
        auto rng = mc::rng::make_engine_with_seed(std::optional<std::uint32_t>{seed}, static_cast<std::uint64_t>(tid));

        for (std::size_t i = 0; i < n_local; ++i) {
            mc::geom::Point<dim> p = proposal.sample(rng);
            if (!domain.isInside(p)) continue;

            const double q = proposal.pdf(p);
            if (q <= 0.0) continue;

            const double term = f(p) / q;
            sum += term;
            sum2 += term * term;
            inside_total += 1;
        }
    }

    ImportanceEstimate<dim> out;
    out.n_samples = n_samples;
    out.n_inside  = inside_total;
    out.mean = sum / static_cast<double>(n_samples);

    const double e2  = sum2 / static_cast<double>(n_samples);
    const double var = std::max(0.0, e2 - out.mean * out.mean);
    out.stderr = std::sqrt(var / static_cast<double>(n_samples));

    return out;
}

} // namespace mc::estimators

#endif // MONTECARLO_DGL_ISMEANESTIMATOR_TPP