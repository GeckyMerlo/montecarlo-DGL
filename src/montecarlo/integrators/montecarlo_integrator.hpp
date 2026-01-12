#ifndef MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP
#define MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP

#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "../mcmc/metropolisHastingsSampler.hpp"
#include "../volumeEstimator/VolumeEstimatorMC.hpp"
#include "../geometry.hpp"
#include "integrator.hpp"
#include <functional>

template <std::size_t dim>
class MontecarloIntegrator : public Integrator<dim> {
public:
    explicit MontecarloIntegrator(const IntegrationDomain<dim> &d);

    // Calcola l'integrale di una funzione 'f' usando Monte Carlo
    double integrate(const std::function<double(const Point<dim>&)>& f, int n_samples);

    double integrate_importance(const std::function<double(const Point<dim>&)>& f, int n_samples, const Proposal<dim>& proposal, uint32_t seed);

    double integrate_with_mh(const std::function<double(const geom::Point<dim>&)>& f,
                                 const std::function<double(const geom::Point<dim>&)>& p,
                                 geom::Point<dim> x0,
                                 double deviation,
                                 std::mt19937& rng,
                                 std::size_t burn_in,
                                 std::size_t n_samples,
                                 std::size_t thinning,
                                 std::size_t n_samples_volume);
};

#include "montecarlo_integrator.tpp"

#endif // MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP
