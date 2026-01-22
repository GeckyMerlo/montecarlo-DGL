// MCMeanEstimator.tpp
//
// Created by Giacomo Merlo on 14/01/26.
//

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

#include <omp.h>

#include "montecarlo/rng/rng_factory.hpp"

template <std::size_t dim>
MeanEstimate<dim> MCMeanEstimator<dim>::estimate(const IntegrationDomain<dim>& domain,
                                                std::uint32_t seed,
                                                std::size_t n_samples,
                                                const std::function<double(const Point<dim>&)>& f) const
{
    if (n_samples == 0) throw std::invalid_argument("n_samples must be > 0");

    auto bounds = domain.getBounds();
    for (std::size_t i = 0; i < dim; ++i) {
        dist[i] = std::uniform_real_distribution<double>(bounds[i].first, bounds[i].second);
    }

    const int T = omp_get_max_threads();
    const std::size_t base = n_samples / static_cast<std::size_t>(T);
    const std::size_t rem  = n_samples % static_cast<std::size_t>(T);

    double sum = 0.0;
    double sum2 = 0.0;
    std::size_t inside_total = 0;

#pragma omp parallel for reduction(+:sum,sum2,inside_total)
    for (int tid = 0; tid < T; ++tid) {
        const std::size_t n_local = base + (static_cast<std::size_t>(tid) < rem ? 1u : 0u);

        auto rng = mc::make_engine_with_seed(std::optional<std::uint32_t>{seed},
                                             static_cast<std::uint64_t>(tid));

        auto local_dist = dist;

        for (std::size_t k = 0; k < n_local; ++k) {
            geom::Point<dim> p;
            for (std::size_t j = 0; j < dim; ++j) {
                p[j] = local_dist[j](rng);
            }

            if (domain.isInside(p)) {
                const double term = f(p);
                sum += term;
                sum2 += term * term;
                inside_total += 1;
            }
        }
    }

    MeanEstimate<dim> out;
    out.n_samples = n_samples;
    out.n_inside  = inside_total;

    out.mean = sum / static_cast<double>(n_samples);
    const double e2  = sum2 / static_cast<double>(n_samples);
    const double var = std::max(0.0, e2 - out.mean * out.mean);
    out.stderr = std::sqrt(var / static_cast<double>(n_samples));

    return out;
}