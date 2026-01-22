// MCintegrator.tpp
#include "integrator.hpp"
#include "../geometry.hpp"

#include <cstdint>
#include <functional>
#include <vector>

using namespace geom;
using namespace std;

template <size_t dim>
MontecarloIntegrator<dim>::MontecarloIntegrator(const IntegrationDomain<dim>& d)
    : Integrator<dim>(d) {}

template <size_t dim>
double MontecarloIntegrator<dim>::OLDintegrate(const function<double(const Point<dim>&)>& f,
                                               int n_samples)
{
    vector<Point<dim>> points = this->initializeRandomizer(n_samples);

    double sum = 0.0;
    for (const auto& p : points) {
        if (this->domain.isInside(p)) sum += f(p);
    }

    const double volume = this->domain.getBoxVolume();
    return (sum / static_cast<double>(n_samples)) * volume;
}

template <size_t dim>
double MontecarloIntegrator<dim>::integrate(const function<double(const Point<dim>&)>& f,
                                            int n_samples,
                                            const Proposal<dim>&,
                                            std::uint32_t seed)
{
    MCMeanEstimator<dim> mean_estimator;
    MeanEstimate<dim> mean_estimate =
        mean_estimator.estimate(this->domain, seed, static_cast<std::size_t>(n_samples), f);

    return mean_estimate.mean * this->domain.getBoxVolume();
}