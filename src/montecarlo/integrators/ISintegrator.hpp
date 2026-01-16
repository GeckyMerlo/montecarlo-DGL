//
// Created by Giacomo Merlo on 15/01/26.
//

#ifndef MONTECARLO_DGL_ISINTEGRATOR_HPP
#define MONTECARLO_DGL_ISINTEGRATOR_HPP

#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "../proposals/uniformProposal.hpp"
#include "../proposals/gaussianProposal.hpp"
#include "../mcmc/metropolisHastingsSampler.hpp"
#include "../estimators/VolumeEstimatorMC.hpp"
#include "../estimators/ISMeanEstimator.hpp"
#include "../estimators/MCMeanEstimator.hpp"
#include "../geometry.hpp"
#include "integrator.hpp"
#include <functional>

template <std::size_t dim>
class ISMontecarloIntegrator : public Integrator<dim> {
public:

    explicit ISMontecarloIntegrator(const IntegrationDomain<dim> &d);

    /**
     * @brief Importance sampling Monte Carlo integration
     * @param f Integrand function
     * @param n_samples Number of samples
     * @param proposal Custom sampling distribution (should approximate f)
     * @param seed Random seed for reproducibility
     * @return Estimated integral with reduced variance
     *
     * Samples from proposal distribution q(x) instead of uniform.
     * Computes ∫ f(x)/q(x) * q(x) dx using importance weights.
     */
    double integrate(const function<double(const Point<dim>&)>& f, int n_samples, const Proposal<dim>& proposal, uint32_t seed) override;

};

#include "ISintegrator.tpp"

#endif //MONTECARLO_DGL_ISINTEGRATOR_HPP