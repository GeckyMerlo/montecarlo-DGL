#ifndef MONTECARLO_1_GAUSSIAN_PROPOSAL_TPP
#define MONTECARLO_1_GAUSSIAN_PROPOSAL_TPP

#include <stdexcept>

template <size_t dim>
GaussianProposal<dim>::GaussianProposal(const IntegrationDomain<dim>& d,
                                        const std::array<double, dim>& mean,
                                        const std::array<double, dim>& sigma)
    : domain(d), mu(mean), sig(sigma)
{
    // Validate sigmas and precompute constants.
    double sum_log_inv_sigma = 0.0;
    for (size_t i = 0; i < dim; ++i) {
        if (sig[i] <= 0.0) {
            throw std::invalid_argument("GaussianProposal: sigma must be > 0 for every dimension.");
        }
        inv_sig2[i] = 1.0 / (sig[i] * sig[i]);
        ndist[i] = std::normal_distribution<double>(mu[i], sig[i]);
        sum_log_inv_sigma += std::log(1.0 / sig[i]);
    }

    // log normalization constant for diagonal Gaussian:
    // phi(x) = (2pi)^(-d/2) * prod(1/sigma_i) * exp(-0.5 * sum((x_i-mu_i)^2 / sigma_i^2))
    const double log_2pi = std::log(2.0 * M_PI);
    log_norm_const = -0.5 * static_cast<double>(dim) * log_2pi + sum_log_inv_sigma;
}

template <size_t dim>
geom::Point<dim> GaussianProposal<dim>::sample(std::mt19937& rng) const
{
    geom::Point<dim> x;

    // Rejection sampling to enforce domain constraint.
    // If the domain is extremely small relative to the Gaussian mass,
    // this can be slow; in that case consider a better-adapted proposal.
    do {
        for (size_t i = 0; i < dim; ++i) {
            x[i] = ndist[i](rng);
        }
    } while (!domain.contains(x));

    return x;
}

template <size_t dim>
double GaussianProposal<dim>::pdf(const geom::Point<dim>& x) const
{
    if (!domain.contains(x)) {
        return 0.0;
    }

    // Compute log phi(x) for numerical stability.
    double quad = 0.0;
    for (size_t i = 0; i < dim; ++i) {
        const double diff = x[i] - mu[i];
        quad += diff * diff * inv_sig2[i];
    }

    const double log_phi = log_norm_const - 0.5 * quad;
    return std::exp(log_phi);
}

#endif // MONTECARLO_1_GAUSSIAN_PROPOSAL_TPP