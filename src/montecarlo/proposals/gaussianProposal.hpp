/**
 * @file gaussianProposal.hpp
 * @brief Diagonal multivariate Gaussian proposal over a generic integration domain.
 *
 * This proposal draws samples from a multivariate Gaussian N(mu, diag(sigma^2))
 * and enforces the domain constraint via rejection sampling:
 *   - sample x ~ N(mu, diag(sigma^2))
 *   - if x ∉ D, resample
 *
 * IMPORTANT:
 * The true proposal induced by rejection is a TRUNCATED Gaussian:
 *   q_trunc(x) = phi(x) * 1_D(x) / Z,  where Z = P(X ∈ D).
 * Computing Z exactly is generally hard for generic domains.
 *
 * Therefore pdf() below returns the *unnormalized* density phi(x) * 1_D(x).
 * This is sufficient for self-normalized IS / MIS where the constant cancels out.
 */

#ifndef MONTECARLO_1_GAUSSIAN_PROPOSAL_HPP
#define MONTECARLO_1_GAUSSIAN_PROPOSAL_HPP

#include "proposal.hpp"
#include "../domains/integration_domain.hpp"

#include <array>
#include <cmath>
#include <vector>
#include <random>

/**
 * @brief Diagonal Gaussian proposal over a domain.
 * @tparam dim Dimensionality.
 *
 * Samples are Gaussian-distributed and then rejected until they lie in the domain.
 * pdf() returns the Gaussian density (diagonal covariance) times the domain indicator.
 */
template <size_t dim>
class GaussianProposal : public Proposal<dim>
{
public:
    /**
     * @brief Construct a diagonal Gaussian proposal over a domain.
     * @param d Integration domain
     * @param mean Mean vector (mu)
     * @param sigma Standard deviations per dimension (sigma_i > 0)
     *
     * The sampling uses rejection to ensure samples lie in the domain.
     */
    GaussianProposal(const IntegrationDomain<dim>& d,
                     const std::vector<double>& mean,
                     const std::vector<double>& sigma);

    /**
     * @brief Sample a point from the (truncated) Gaussian over the domain.
     * @param rng Random generator
     * @return A point in the domain
     *
     * Rejection sampling: keep drawing from N(mu, diag(sigma^2)) until inside domain.
     */
    geom::Point<dim> sample(std::mt19937& rng) const override;

    /**
     * @brief Evaluate the proposal "pdf" at x.
     * @param x Query point
     * @return phi(x; mu, diag(sigma^2)) if x in domain, otherwise 0.
     *
     * Note: this is UNNORMALIZED if you use rejection (missing 1/Z).
     * Use self-normalized IS/MIS or provide/estimate Z externally if needed.
     */
    double pdf(const geom::Point<dim>& x) const override;

private:
    const IntegrationDomain<dim>& domain;

    std::vector<double> mu{};
    std::vector<double> sig{};
    std::vector<double> inv_sig2{};   // 1/sigma^2 for each dimension
    double log_norm_const;                // log((2pi)^(-d/2) * prod(1/sigma_i))

    // One normal distribution per dimension (diagonal covariance).
    mutable std::array<std::normal_distribution<double>, dim> ndist{};
};

#include "gaussianProposal.tpp"

#endif // MONTECARLO_1_GAUSSIAN_PROPOSAL_HPP