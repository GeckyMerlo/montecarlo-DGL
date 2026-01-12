//
// Created by Giacomo Merlo on 12/01/26.
//

template <std::size_t dim>
static geom::Point<dim>

sample_uniform_in_box(const Bounds<dim>& b, std::mt19937& rng)
{
    geom::Point<dim> x;
    for (std::size_t k = 0; k < dim; ++k) {
        std::uniform_real_distribution<double> unif(b[k].first, b[k].second);
        x[k] = unif(rng);
    }
    return x;
}

template <std::size_t dim>
VolumeEstimate<dim>
VolumeEstimatorMC<dim>::estimate(const IntegrationDomain<dim>& domain,
                                 std::mt19937& rng,
                                 std::size_t n_samples) const
{
    if (n_samples == 0)
        throw std::invalid_argument("VolumeEstimatorMC: n_samples must be > 0.");

    const geom::Bounds<dim> bounds = domain.getBounds();
    const double boxV = domain.getBoxVolume();

    std::size_t inside = 0;
    for (std::size_t i = 0; i < n_samples; ++i) {
        const auto x = sample_uniform_in_box<dim>(bounds, rng);
        if (domain.isInside(x)) ++inside;
    }

    const double p = static_cast<double>(inside) / static_cast<double>(n_samples);

    // Bernoulli std error for p-hat
    const double var_p = (p * (1.0 - p)) / static_cast<double>(n_samples);
    const double se_p  = std::sqrt(var_p);

    VolumeEstimate<dim> out;
    out.n_samples = n_samples;
    out.inside_ratio = p;
    out.volume = boxV * p;
    out.stderr = boxV * se_p;
    return out;
}
