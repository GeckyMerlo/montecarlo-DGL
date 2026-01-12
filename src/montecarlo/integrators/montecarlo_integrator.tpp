#include "integrator.hpp"
#include "../geometry.hpp"
#include <vector>
#include <fstream>
#include <iostream>

using namespace geom;
using namespace std;
// Costruttore
template <size_t dim>
MontecarloIntegrator<dim>::MontecarloIntegrator(const IntegrationDomain<dim> &d)
    : Integrator<dim>(d) {}

// Funzione di integrazione Monte Carlo
template <size_t dim>
double MontecarloIntegrator<dim>::integrate(
    const function<double(const Point<dim>&)> &f,
    int n_samples)
{
    // Genero n_samples punti casuali nel dominio
    vector<Point<dim>> points = this->initializeRandomizer(n_samples);

    // Somma dei valori della funzione nei punti generati
    double sum = 0.0;
    for (const auto& p : points) {
        if (this->domain.isInside(p)) {
                    sum += f(p);
        }
    }

    // Calcolo volume del dominio
    double volume = this->domain.getBoxVolume();

    // Restituisco l’integrale stimato
    return (sum / n_samples) * volume;
}

template <size_t dim>
double MontecarloIntegrator<dim>::integrate_importance(
    const function<double(const Point<dim>&)>& f,
    int n_samples,
    const Proposal<dim>& proposal,
    uint32_t seed)
{
    std::mt19937 rng(seed);
    double sum = 0.0;

    for (int i=0; i<n_samples; ++i) {
        Point<dim> p = proposal.sample(rng);

        if (this->domain.isInside(p)) {
            double q = proposal.pdf(p);
            if (q > 0.0) sum += f(p)/q;
        }
    }
    return sum / n_samples;
}

template <size_t dim>
double MontecarloIntegrator<dim>::integrate_with_mh(
                                 const std::function<double(const geom::Point<dim>&)>& f,
                                 const std::function<double(const geom::Point<dim>&)>& p,
                                 geom::Point<dim> x0,
                                 double deviation,
                                 std::mt19937& rng,
                                 std::size_t burn_in,
                                 std::size_t n_samples,
                                 std::size_t thinning,
                                 std::size_t n_samples_volume)
{
    if (n_samples == 0) throw std::invalid_argument("n_samples must be > 0");
    if (thinning == 0) throw std::invalid_argument("thinning must be > 0");
    if (!this->domain.isInside(x0)) throw std::invalid_argument("x0 must be inside the domain.");

    // 1) stima volume |D|
    VolumeEstimatorMC<dim> ve;
    const VolumeEstimate<dim> vol_hat = ve.estimate(this->domain, rng, n_samples_volume);
    std::cout << "Volume estimation results: " << vol_hat.volume << " +- " << 2 * vol_hat.stderr << std::endl;


    MetropolisHastingsSampler<dim> mh(this->domain, p, x0, deviation);

    // burn-in
    for (std::size_t i = 0; i < burn_in; ++i) {
        (void) mh.next(rng);
    }

    // campionamento con thinning
    long double sum = 0.0;
    std::size_t kept = 0;

    geom::Point<dim> x;
    for (std::size_t i = 0; i < n_samples; ++i) {
        for (std::size_t t = 0; t < thinning; ++t)
            x = mh.next(rng);
        const double fx = f(x);

        if (std::isfinite(fx)) {
            sum += static_cast<long double>(fx);
            ++kept;
        }
    }

    if (kept == 0)
        throw std::runtime_error("All sampled f(x) were non-finite. Cannot compute integral.");

    const double mean_f_hat = static_cast<double>(sum / static_cast<long double>(kept));

    // integrale ≈ |D| * media_uniforme(f)
    return vol_hat.volume * mean_f_hat;
}

