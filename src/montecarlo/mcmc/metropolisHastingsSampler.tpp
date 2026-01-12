//
// Created by Giacomo Merlo on 11/01/26.
//


#include <cmath>
#include <algorithm>
#include <limits>

template <size_t dim>
MetropolisHastingsSampler<dim>::MetropolisHastingsSampler(
    const IntegrationDomain<dim>& d,
    const std::function<double(const geom::Point<dim>&)>& p,
    geom::Point<dim> x0,
    double deviation)
  : domain(d)
  , target(p)
  , current(x0)
  , rw_normal(0.0, deviation)
  , uni(0.0, 1.0)
{
    if (!domain.isInside(current)) {
        throw std::invalid_argument("MetropolisHastingsSampler: x0 is outside the domain.");
    }
}

template <size_t dim>
geom::Point<dim>
MetropolisHastingsSampler<dim>::next(std::mt19937& rng)
{
    ++n_steps;

    //Genero y
    geom::Point<dim> y = current;
    for (std::size_t k = 0; k < dim; ++k)
        y[k] += rw_normal(rng);

    //Verifico che sia all'interno del dominio se no lo rifiuto subito
    /*
    if (!domain.isInside(y)) {
        return current;
    }

    Tolto perchè conviene dare in pasto al metodo una funzione che ritorni 0 fuori dal dominio
    target del tipo return domain.isInside(x) ? 1.0 : 0.0;
    */

    const double px = target(current);
    const double py = target(y);

    if (!std::isfinite(py) || py <= 0.0) {
        return current; // reject
    }

    double alpha = 0.0;


    alpha = std::min(1.0, py / px);


    if (uni(rng) < alpha) {
        current = y;
        ++n_accept;
    }
    return current;
}

template <size_t dim>
double
MetropolisHastingsSampler<dim>::target_pdf(const geom::Point<dim>& x)
{
    return target(x);
}

template <size_t dim>
double
MetropolisHastingsSampler<dim>::acceptance_rate() const
{
    return (n_steps == 0) ? 0.0 : static_cast<double>(n_accept) / static_cast<double>(n_steps);
}