#ifndef MONTECARLO_1_UNIFORM_PROPOSAL_HPP
#define MONTECARLO_1_UNIFORM_PROPOSAL_HPP

#include "proposal.hpp"
#include "../domains/integration_domain.hpp"

#include <array>
#include <random>
#include <utility>

template <size_t dim>
class UniformProposal : public Proposal<dim>
{
public:
    explicit UniformProposal(const IntegrationDomain<dim>& d);

    geom::Point<dim> sample(std::mt19937& rng) const override;

    double pdf(const geom::Point<dim>&) const override;

private:
    const IntegrationDomain<dim>& domain;
    mutable std::array<std::uniform_real_distribution<double>, dim> dist{};
    double vol_box;
};

#include "uniformProposal.tpp"

#endif // MONTECARLO_1_UNIFORM_PROPOSAL_HPP