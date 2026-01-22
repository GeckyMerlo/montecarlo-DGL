/**
 * @file ISMeanEstimator.hpp
 * @brief Importance sampling mean estimator.
 * @author Giacomo Merlo
 * @date 12/01/26
 *
 * @details Estimates the weighted mean of a function using importance sampling.
 * For integration: ∫_Ω f(x) dx = ∫ [f(x)/q(x)] · q(x) dx ≈ (1/N) ∑ [f(xᵢ)/q(xᵢ)]
 * where xᵢ ~ q (proposal distribution).
 */

#ifndef MONTECARLO_DGL_ISMEANESTIMATOR_HPP
#define MONTECARLO_DGL_ISMEANESTIMATOR_HPP

#include <cstddef>
#include <cstdint>
#include <functional>

namespace mc::estimators {

/**
 * @struct ImportanceEstimate
 * @brief Result of importance sampling mean estimation.
 * @tparam dim Dimensionality of the domain.
 */
template <std::size_t dim>
struct ImportanceEstimate {
    /** Estimated importance-weighted mean: μ̂ = (1/N) ∑ [f(xᵢ)/q(xᵢ)] */
    double mean = 0.0;
    /** Standard error of the weighted mean. */
    double stderr = 0.0;
    /** Total samples generated from proposal. */
    std::size_t n_samples = 0;
    /** Samples that fell inside domain (and had q > 0). */
    std::size_t n_inside = 0;
};

/**
 * @class ISMeanEstimator
 * @brief Importance sampling mean estimator with variance reduction.
 * @tparam dim Dimensionality of the space.
 * 
 * @details Estimates the mean using importance sampling:
 * 1. Sample N points from proposal q(x): x₁, ..., xₙ ~ q
 * 2. Filter to domain points with q > 0
 * 3. Compute weighted average: μ̂ = (1/N) ∑ᵢ [f(xᵢ)/q(xᵢ)]
 * 4. Estimate variance and standard error
 * 
 * **Variance Reduction**: If q(x) ≈ f(x), the weights [f(x)/q(x)] are roughly
 * constant, reducing variance compared to uniform sampling.
 * 
 * **Usage**: For integration over Ω:
 * ∫_Ω f(x) dx ≈ (mean result) [no volume factor—importance weight handles normalization]
 * 
 * @note Uses OpenMP for parallel sampling from proposal.
 * @note Proposal must support sample() and pdf() methods.
 */
template <std::size_t dim>
class ISMeanEstimator {
public:
    /**
     * @brief Estimate the weighted mean using importance sampling.
     * @param domain The integration domain.
     * @param seed Random seed for reproducibility.
     * @param n_samples Number of samples to draw from proposal.
     * @param proposal Sampling distribution q(x) (should resemble |f(x)|).
     * @param f Function to evaluate.
     * @return ImportanceEstimate with weighted mean, stderr, and sample counts.
     * @throws std::invalid_argument if n_samples == 0.
     * 
     * @details Algorithm:
     * - Each thread samples n_local points from proposal (parallel)
     * - Keeps only points inside domain with q(x) > 0
     * - Accumulates weighted sum ∑ [f(x)/q(x)] and sum of squares
     * - Returns (mean / N, stderr, n_inside)
     */
    ImportanceEstimate<dim> estimate(const mc::domains::IntegrationDomain<dim>& domain,
                                     std::uint32_t seed,
                                     std::size_t n_samples,
                                     const mc::proposals::Proposal<dim>& proposal,
                                     const std::function<double(const mc::geom::Point<dim>&)>& f) const;
};

} // namespace mc::estimators

#include "ISMeanEstimator.tpp"

#endif // MONTECARLO_DGL_ISMEANESTIMATOR_HPP