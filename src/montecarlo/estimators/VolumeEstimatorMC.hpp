/**
 * @file VolumeEstimatorMC.hpp
 * @brief Hit-or-Miss Monte Carlo volume estimation.
 * @author Giacomo Merlo
 * @date 12/01/26
 *
 * @details Estimates the volume of complex N-dimensional domains using
 * acceptance-rejection (hit-or-miss) sampling:
 * 
 * \f[
 * V_{\text{est}} = V_{\text{box}} \cdot \frac{N_{\text{hits}}}{N_{\text{total}}}
 * \f]
 * 
 * The standard error is computed from the Bernoulli distribution variance:
 * \f[
 * \sigma_V = V_{\text{box}} \cdot \sqrt{\frac{\hat{p}(1-\hat{p})}{N}}
 * \f]
 */

#ifndef MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP
#define MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <stdexcept>

#include "../domains/integration_domain.hpp"
#include "../geometry.hpp"

namespace mc::estimators {

/**
 * @struct VolumeEstimate
 * @brief Result of Monte Carlo volume estimation.
 * @tparam dim Dimensionality of the domain.
 */
template <std::size_t dim>
struct VolumeEstimate {
    /** Estimated domain volume: V̂ = V_box · p̂ */
    double volume = 0.0;
    /** Standard error in the volume estimate. */
    double stderr = 0.0;
    /** Fraction of samples inside domain: p̂ = N_hits / N_total */
    double inside_ratio = 0.0;
    /** Total number of samples evaluated. */
    std::size_t n_samples = 0;
};

/**
 * @class VolumeEstimatorMC
 * @brief Hit-or-Miss Monte Carlo volume estimator.
 * @tparam dim Dimensionality of the space.
 * 
 * @details Estimates the volume of an arbitrary domain Ω by sampling uniformly
 * in the bounding box and counting the fraction of samples inside Ω.
 * 
 * Algorithm:
 * 1. Sample N points uniformly from bounding box B
 * 2. Count hits: N_hits = |{x ∈ B : x ∈ Ω}|
 * 3. Estimate: V̂_Ω = V_B · (N_hits / N)
 * 4. Compute standard error from Bernoulli variance
 * 
 * @note This estimator is unbiased but can be slow for domains with small volume.
 * @note Used internally by integrators for two-stage integration (e.g., MH-based).
 */
template <std::size_t dim>
class VolumeEstimatorMC {
public:
    /**
     * @brief Estimate the volume of a domain using hit-or-miss sampling.
     * @param domain The integration domain whose volume to estimate.
     * @param seed Random seed for reproducibility.
     * @param n_samples Number of sample points to evaluate.
     * @return VolumeEstimate containing volume, stderr, hit ratio, and sample count.
     * @throws std::invalid_argument if n_samples == 0.
     * 
     * @details Uses OpenMP parallelization for efficient parallel sampling.
     * Each thread gets a deterministic RNG stream based on seed and thread ID.
     */
    VolumeEstimate<dim>
    estimate(const mc::domains::IntegrationDomain<dim>& domain,
             std::uint32_t seed,
             std::size_t n_samples) const;
};

} // namespace mc::estimators

#include "VolumeEstimatorMC.tpp"

#endif // MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP