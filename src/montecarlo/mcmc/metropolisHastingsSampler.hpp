//
// Created by Giacomo Merlo on 11/01/26.
//

#ifndef MONTECARLO_DGL_METROPOLISHASTINGSSAMPLER_HPP
#define MONTECARLO_DGL_METROPOLISHASTINGSSAMPLER_HPP

#include "../domains/integration_domain.hpp"
#include "../geometry.hpp"

#include <array>
#include <random>
#include <utility>
#include <functional>

template <size_t dim>
class MetropolisHastingsSampler
{
public:
    explicit MetropolisHastingsSampler(const IntegrationDomain<dim>& d,
                                            const std::function<double(const geom::Point<dim>&)>& p,
                                            geom::Point<dim> x0,
                                            double deviation);

    geom::Point<dim> next(std::mt19937& rng);

    double target_pdf(const geom::Point<dim>& x);

    double acceptance_rate() const;

private:
    const IntegrationDomain<dim>& domain;

    std::function<double(const geom::Point<dim>&)> target;

    geom::Point<dim> current;

    std::size_t n_steps = 0;
    std::size_t n_accept = 0;

    std::normal_distribution<double> rw_normal;
    std::uniform_real_distribution<double> uni;
};

#include "metropolisHastingsSampler.tpp"


#endif //MONTECARLO_DGL_METROPOLISHASTINGSSAMPLER_HPP