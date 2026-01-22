/**
 * @file MCMeanEstimator.hpp
 * @brief Classic uniform Monte Carlo mean estimator.
 * @author Giacomo Merlo
 * @date 14/01/26
 *
 * @details Estimates the mean of a function over a domain using uniform sampling.
 * For integration: ∫_Ω f(x) dx ≈ V_Ω · (1/N) ∑ f(xᵢ) where xᵢ ~ Unif(box).
 */

#ifndef MONTECARLO_DGL_MCMEANESTIMATOR_HPP
#define MONTECARLO_DGL_MCMEANESTIMATOR_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>

namespace mc::estimators {

/**
 * @struct MeanEstimate
 * @brief Result of Monte Carlo mean estimation.
 * @tparam dim Dimensionality of the domain.
 */
template <std::size_t dim>
struct MeanEstimate {
    /** Estimated mean: μ̂ = (1/N) ∑ f(xᵢ) [only domain points] */
    double mean   = 0.0;
    /** Standard error (i.i.d. assumption): σ/√N */
    double stderr = 0.0;
    /** Total samples evaluated (both inside and outside domain). */
    std::size_t n_samples = 0;
    /** Samples that fell inside the domain (used for mean). */
    std::size_t n_inside  = 0;
};

/**
 * @class MCMeanEstimator
 * @brief Uniform Monte Carlo mean estimator.
 * @tparam dim Dimensionality of the space.
 * 
 * @details Estimates the mean of f over domain Ω by:
 * 1. Sampling N points uniformly in bounding box
 * 2. Filtering to points inside Ω: x₁, ..., x_M where M ≤ N
 * 3. Computing μ̂ = (1/N) ∑ᵢ f(xᵢ) · 𝟙[xᵢ ∈ Ω]
 * 4. Estimating variance: σ̂² = E[f²] - μ̂²
 * 5. Standard error: SE = σ̂/√N
 * 
 * The mean (before multiplying by volume) is the empirical average over
 * all N trials (treating out-of-domain points as contributing 0).
 * 
 * @note Uses parallel sampling with OpenMP for efficiency.
 */
template <std::size_t dim>
class MCMeanEstimator {
public:
    /**
     * @brief Estimate the mean of a function over a domain using uniform sampling.
     * @param domain The integration domain.
     * @param seed Random seed for reproducibility.
     * @param n_samples Total number of sample points to generate.
     * @param f Function to evaluate (only at points inside domain).
     * @return MeanEstimate with mean, stderr, and sample counts.
     * @throws std::invalid_argument if n_samples == 0.
     * 
     * @details Parallel sampling: each thread gets a deterministic RNG stream
     * and samples its portion of n_samples uniformly in the domain's bounding box.
     * Only points inside the domain contribute to the sum.
     */
    MeanEstimate<dim> estimate(const mc::domains::IntegrationDomain<dim>& domain,
                              std::uint32_t seed,
                              std::size_t n_samples,
                              const std::function<double(const mc::geom::Point<dim>&)>& f) const;

private:
    /** Per-dimension uniform distributions (initialized per call). */
    mutable std::array<std::uniform_real_distribution<double>, dim> dist{};
};

} // namespace mc::estimators

#include "MCMeanEstimator.tpp"

#endif // MONTECARLO_DGL_MCMEANESTIMATOR_HPP