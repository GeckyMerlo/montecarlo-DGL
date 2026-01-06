#include "integrator.hpp"
#include "../geometry.hpp"
#include <vector>
#include <fstream>

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

