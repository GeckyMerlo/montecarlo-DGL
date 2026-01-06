#ifndef MONTECARLO_1_UNIFORM_PROPOSAL_TPP
#define MONTECARLO_1_UNIFORM_PROPOSAL_TPP

template <size_t dim>
UniformProposal<dim>::UniformProposal(const IntegrationDomain<dim>& d)
    : domain(d), vol_box(1.0)
{
    auto bounds = domain.getBounds();
    for (size_t i = 0; i < dim; ++i) {
        dist[i] = std::uniform_real_distribution<double>(
            bounds[i].first, bounds[i].second
        );
    }
    vol_box = domain.getBoxVolume();
}

template <size_t dim>
geom::Point<dim> UniformProposal<dim>::sample(std::mt19937& rng) const
{
    geom::Point<dim> x;
    for (size_t i = 0; i < dim; ++i) {
        x[i] = dist[i](rng);
    }
    return x;
}

template <size_t dim>
double UniformProposal<dim>::pdf(const geom::Point<dim>&) const
{
    return 1.0 / vol_box;
}

#endif // MONTECARLO_1_UNIFORM_PROPOSAL_TPP