#ifndef MONTECARLO_1_PROPOSAL_HPP
#define MONTECARLO_1_PROPOSAL_HPP

#include <random>
#include "../geometry.hpp"

template <size_t dim>
class Proposal
{
public:
    virtual ~Proposal() = default;

    // Sample a point according to p(x)
    virtual geom::Point<dim> sample(std::mt19937& rng) const = 0;

    // Evaluate p(x) (the PDF) at x
    virtual double pdf(const geom::Point<dim>& x) const = 0;
};

#endif