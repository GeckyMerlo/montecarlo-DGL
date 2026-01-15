//
// Created by Giacomo Merlo on 15/01/26.
//
#ifndef MONTECARLO_DGL_MHINTEGRATOR_HPP
#define MONTECARLO_DGL_MHINTEGRATOR_HPP

#include "integrator.hpp"
#include "../domains/integration_domain.hpp"
#include "../proposals/proposal.hpp"
#include "../mcmc/metropolisHastingsSampler.hpp"
#include "../estimators/VolumeEstimatorMC.hpp"
#include "../geometry.hpp"

#include <functional>
#include <cstddef>
#include <cstdint>

template <std::size_t dim>
class MHMontecarloIntegrator : public Integrator<dim> {
public:
    using Point = geom::Point<dim>;
    using Func  = std::function<double(const Point&)>;

    explicit MHMontecarloIntegrator(const IntegrationDomain<dim> &d);

    double integrate(const Func& f,
                     int n_samples,
                     const Proposal<dim>& proposal,
                     std::uint32_t seed) override;

    void setConfig(std::size_t burn_in_,
                   std::size_t thinning_,
                   std::size_t n_samples_volume_,
                   double deviation_,
                   Func p_,
                   Point x0_);

private:
    bool configured = false;

    std::size_t burn_in = 0;
    std::size_t thinning = 1;
    std::size_t n_samples_volume = 0;
    double deviation = 1.0;

    std::function<double(const Point&)> p;
    geom::Point<dim> x0;
};

#include "MHintegrator.tpp"
#endif