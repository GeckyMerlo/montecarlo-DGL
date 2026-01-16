// mixtureProposal.tpp
#ifndef MONTECARLO_1_MIXTURE_PROPOSAL_TPP
#define MONTECARLO_1_MIXTURE_PROPOSAL_TPP

#include <numeric>   // std::accumulate
#include <cmath>     // std::isfinite

template <size_t dim>
void MixtureProposal<dim>::validateInputs(const std::vector<const Proposal<dim>*>& components,
                                          const std::vector<double>& weights)
{
    if (components.empty()) {
        throw std::invalid_argument("MixtureProposal: components must be non-empty.");
    }
    if (components.size() != weights.size()) {
        throw std::invalid_argument("MixtureProposal: components and weights must have the same size.");
    }

    for (std::size_t i = 0; i < components.size(); ++i) {
        if (components[i] == nullptr) {
            throw std::invalid_argument("MixtureProposal: component pointer cannot be null.");
        }
        if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
            throw std::invalid_argument("MixtureProposal: weights must be finite and >= 0.");
        }
    }

    const double sumw = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (!std::isfinite(sumw) || sumw <= 0.0) {
        throw std::invalid_argument("MixtureProposal: sum of weights must be > 0 and finite.");
    }
}

template <size_t dim>
std::vector<double> MixtureProposal<dim>::normalizeWeights(const std::vector<double>& weights)
{
    const double sumw = std::accumulate(weights.begin(), weights.end(), 0.0);

    std::vector<double> out;
    out.reserve(weights.size());

    for (double wi : weights) {
        out.push_back(wi / sumw);
    }
    return out;
}

template <size_t dim>
MixtureProposal<dim>::MixtureProposal(std::vector<const Proposal<dim>*> components,
                                      std::vector<double> weights)
    : comps(std::move(components))
{
    validateInputs(comps, weights);
    w = normalizeWeights(weights);

    // Build categorical distribution for sampling component indices.
    cat = std::discrete_distribution<std::size_t>(w.begin(), w.end());
}

template <size_t dim>
geom::Point<dim> MixtureProposal<dim>::sample(std::mt19937& rng) const
{
    const std::size_t k = cat(rng);
    return comps[k]->sample(rng);
}

template <size_t dim>
double MixtureProposal<dim>::pdf(const geom::Point<dim>& x) const
{
    double acc = 0.0;
    for (std::size_t k = 0; k < comps.size(); ++k) {
        acc += w[k] * comps[k]->pdf(x);
    }
    return acc;
}

#endif // MONTECARLO_1_MIXTURE_PROPOSAL_TPP