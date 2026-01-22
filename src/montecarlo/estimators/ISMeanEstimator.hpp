// ISMeanEstimator.hpp

/**
 * @file ISMeanEstimator.hpp
 * @brief Importance sampling mean estimator
 * @author Giacomo Merlo
 * @date 12/01/26
 *
 * Computes the mean of weighted samples for importance sampling integration.
 * Estimates ∫ f(x) dx = V * E_q[f(X)/q(X)] where q is the proposal and V is volume.
 */

#ifndef MONTECARLO_DGL_ISMEANESTIMATOR_HPP
#define MONTECARLO_DGL_ISMEANESTIMATOR_HPP

#include <cstddef>
#include <cstdint>
#include <functional>

template <std::size_t dim>
struct ImportanceEstimate {
    double mean = 0.0;
    double stderr = 0.0;
    std::size_t n_samples = 0;
    std::size_t n_inside = 0;
};

template <std::size_t dim>
class ISMeanEstimator {
public:
    ImportanceEstimate<dim> estimate(const IntegrationDomain<dim>& domain,
                                    std::uint32_t seed,
                                    std::size_t n_samples,
                                    const Proposal<dim>& proposal,
                                    const std::function<double(const Point<dim>&)>& f) const;
};

#include "ISMeanEstimator.tpp"

#endif // MONTECARLO_DGL_ISMEANESTIMATOR_HPP