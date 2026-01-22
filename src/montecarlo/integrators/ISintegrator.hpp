/**
 * @file ISintegrator.hpp
 * @brief Importance sampling Monte Carlo integration engine.
 * @details Implements variance-reduced integration by sampling from a custom
 * proposal distribution that approximates the integrand.
 */

#ifndef MONTECARLO_DGL_ISINTEGRATOR_HPP
#define MONTECARLO_DGL_ISINTEGRATOR_HPP

#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "../proposals/uniformProposal.hpp"
#include "../proposals/gaussianProposal.hpp"
#include "../proposals/mixtureProposal.hpp"
#include "../mcmc/metropolisHastingsSampler.hpp"
#include "../estimators/VolumeEstimatorMC.hpp"
#include "../estimators/ISMeanEstimator.hpp"
#include "../estimators/MCMeanEstimator.hpp"
#include "../geometry.hpp"
#include "integrator.hpp"
#include <functional>

namespace mc::integrators {

/**
 * @class ISMontecarloIntegrator
 * @brief Importance sampling Monte Carlo integrator.
 * @tparam dim Dimensionality of the integration domain.
 * 
 * @details Computes the integral using importance sampling:
 * ∫_Ω f(x) dx = ∫_Ω [f(x)/q(x)] · q(x) dx ≈ (1/N) ∑ᵢ [f(xᵢ)/q(xᵢ)]
 * 
 * where q(x) is a proposal distribution chosen to approximate f(x).
 * This reduces variance compared to uniform sampling when q resembles f.
 * 
 * @note The proposal should be chosen carefully to match the integrand's shape.
 */
template <std::size_t dim>
class ISMontecarloIntegrator : public Integrator<dim> {
public:

    /**
     * @brief Construct an importance sampling integrator.
     * @param d Reference to the integration domain.
     */
    explicit ISMontecarloIntegrator(const mc::domains::IntegrationDomain<dim> &d);

    /**
     * @brief Compute the integral using importance sampling.
     * @param f Integrand function: ℝⁿ → ℝ.
     * @param n_samples Number of samples drawn from the proposal.
     * @param proposal Custom sampling distribution q(x).
     *                 Should approximate f(x) to minimize variance.
     * @param seed Random seed for reproducibility.
     * @return Estimated integral ∫_Ω f(x) dx with reduced variance.
     * 
     * @details Uses ISMeanEstimator to compute weighted average:
     * (1/N) ∑ᵢ [f(xᵢ)/q(xᵢ)] where xᵢ ~ q.
     */
    double integrate(const std::function<double(const mc::geom::Point<dim>&)>& f,
                     int n_samples,
                     const mc::proposals::Proposal<dim>& proposal,
                     std::uint32_t seed) override;

};

} // namespace mc::integrators

#include "ISintegrator.tpp"

#endif //MONTECARLO_DGL_ISINTEGRATOR_HPP