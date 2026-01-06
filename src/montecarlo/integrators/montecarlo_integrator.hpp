#ifndef MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP
#define MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP

#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "integrator.hpp"
#include <functional>

template <std::size_t dim>
class MontecarloIntegrator : public Integrator<dim> {
public:
    explicit MontecarloIntegrator(const IntegrationDomain<dim> &d);

    // Calcola l'integrale di una funzione 'f' usando Monte Carlo
    double integrate(const std::function<double(const Point<dim>&)>& f, int n_samples);

    double integrate_importance(const std::function<double(const Point<dim>&)>& f, int n_samples, const Proposal<dim>& proposal, uint32_t seed);
};

#include "montecarlo_integrator.tpp"

#endif // MONTECARLO_1_MONTECARLO_INTEGRATOR_HPP
