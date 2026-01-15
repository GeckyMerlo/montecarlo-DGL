#include "integrator.hpp"
#include "../geometry.hpp"
#include "../RngManager.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <omp.h>

using namespace geom;
using namespace std;
// Costruttore
template <size_t dim>
MontecarloIntegrator<dim>::MontecarloIntegrator(const IntegrationDomain<dim> &d)
    : Integrator<dim>(d) {}

// Funzione di integrazione Monte Carlo
template <size_t dim>
double MontecarloIntegrator<dim>::OLDintegrate(
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
double MontecarloIntegrator<dim>::integrate(
    const function<double(const Point<dim>&)>& f,
    int n_samples,
    const Proposal<dim>& proposal,
    uint32_t seed)
{
    MCMeanEstimator<dim> mean_estimator;
    MeanEstimate<dim> mean_estimate = mean_estimator.estimate(this->domain, seed, n_samples, proposal, f);
    return mean_estimate.mean * this->domain.getBoxVolume();
}



