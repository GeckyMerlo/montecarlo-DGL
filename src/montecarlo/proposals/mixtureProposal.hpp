/**
 * @file mixtureProposal.hpp
 * @brief Mixture proposal distribution for importance sampling (non-owning raw pointers).
 *
 * A mixture proposal is defined as:
 *   q(x) = sum_{k=1..K} w_k * q_k(x)
 * where w_k >= 0 and sum_k w_k = 1, and each q_k is a Proposal<dim>.
 *
 * Sampling:
 *   1) draw k ~ Categorical(w)
 *   2) draw x ~ q_k
 *
 * PDF evaluation:
 *   q(x) = sum_k w_k * q_k.pdf(x)
 *
 * LIFETIME NOTE:
 * This class stores NON-OWNING pointers to proposal components.
 * The caller must guarantee that all component proposals outlive this MixtureProposal.
 */

#ifndef MONTECARLO_1_MIXTURE_PROPOSAL_HPP
#define MONTECARLO_1_MIXTURE_PROPOSAL_HPP

#include "proposal.hpp"

#include <vector>
#include <random>
#include <stdexcept>
#include <numeric>   // std::accumulate
#include <cmath>     // std::isfinite

template <size_t dim>
class MixtureProposal : public Proposal<dim>
{
public:
    /**
     * @brief Construct a mixture proposal from non-owning component pointers and weights.
     * @param components Vector of non-null pointers to Proposal<dim> components (q_k).
     * @param weights Vector of non-negative weights (w_k). Will be normalized to sum to 1.
     *
     * Requirements:
     * - components.size() == weights.size()
     * - components.size() > 0
     * - all weights >= 0 and sum(weights) > 0
     * - no null component pointers
     *
     * IMPORTANT: components are NOT owned by this class.
     */
    MixtureProposal(std::vector<const Proposal<dim>*> components,
                    std::vector<double> weights);

    /// @brief Sample from the mixture.
    geom::Point<dim> sample(std::mt19937& rng) const override;

    /// @brief Evaluate mixture PDF q(x) = sum_k w_k * q_k(x).
    double pdf(const geom::Point<dim>& x) const override;

    std::size_t numComponents() const noexcept { return comps.size(); }
    const std::vector<double>& getWeights() const noexcept { return w; }

private:
    std::vector<const Proposal<dim>*> comps; // non-owning pointers
    std::vector<double> w;                   // normalized weights

    mutable std::discrete_distribution<std::size_t> cat;

    static void validateInputs(const std::vector<const Proposal<dim>*>& components,
                               const std::vector<double>& weights);
    static std::vector<double> normalizeWeights(const std::vector<double>& weights);
};

#include "mixtureProposal.tpp"

#endif // MONTECARLO_1_MIXTURE_PROPOSAL_HPP