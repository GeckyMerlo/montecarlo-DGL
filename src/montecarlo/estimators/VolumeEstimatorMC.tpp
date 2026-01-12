//
// Created by Giacomo Merlo on 12/01/26.
//
#include <omp.h>
#include "../rngManager.hpp"

template <std::size_t dim>
static geom::Point<dim>
sample_uniform_in_box(const geom::Bounds<dim>& b, std::mt19937& rng)
{
    geom::Point<dim> x;
    for (std::size_t k = 0; k < dim; ++k) {
        std::uniform_real_distribution<double> unif(b[k].first, b[k].second);
        x[k] = unif(rng);
    }
    return x;
}

template <std::size_t dim>
VolumeEstimate<dim>
VolumeEstimatorMC<dim>::estimate(const IntegrationDomain<dim>& domain,
                                 uint32_t seed,
                                 std::size_t n_samples) const
{
    if (n_samples == 0)
        throw std::invalid_argument("VolumeEstimatorMC: n_samples must be > 0.");

    const geom::Bounds<dim> bounds = domain.getBounds();
    const double boxV = domain.getBoxVolume();

    const int T = omp_get_max_threads();
    const size_t base = n_samples / T;
    const size_t rem  = n_samples % T;

    RngManager rngs(seed);
    std::size_t inside = 0;

    #pragma omp parallel for reduction(+:inside)
    for (int tid = 0; tid < T; ++tid) {
        const size_t n_local = base + (tid < rem ? 1 : 0);
        auto rng = rngs.make_rng(tid);
        
        for (size_t i = 0; i < n_local; ++i) {
            const auto x = sample_uniform_in_box<dim>(bounds, rng);
            if (domain.isInside(x)) ++inside;
        }
    }

    const double p = static_cast<double>(inside) / static_cast<double>(n_samples);

    // Bernoulli std error for p-hat
    const double var_p = (p * (1.0 - p)) / static_cast<double>(n_samples);
    const double se_p  = std::sqrt(var_p);

    VolumeEstimate<dim> out;
    out.n_samples = n_samples;
    out.inside_ratio = p;
    out.volume = boxV * p;
    out.stderr = boxV * se_p;
    return out;
}
