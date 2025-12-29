#ifndef MONTECARLO_1_INTEGRATOR_HPP
#define MONTECARLO_1_INTEGRATOR_HPP

#include <random>
#include <functional>
#include <array>
#include <fstream>
#include <vector>

#include "../domains/integration_domain.hpp"
#include "../domains/hypersphere.hpp"
#include "../domains/hypercylinder.hpp"
#include "../domains/hyperrectangle.hpp"

using namespace std;
using namespace geom;

template <size_t dim>
class Integrator
{
protected:
    const IntegrationDomain<dim> &domain;

    vector<mt19937> randomizer;

    vector<Point<dim>> initializeRandomizer(int numbers)
    {
        // Inizializzo i seed per i generatori di numeri casuali
        seed_seq seq{1, 2, 3, 4, 5};
        vector<uint32_t> seeds(dim);
        seq.generate(seeds.begin(), seeds.end());

        // Creo 'dim' generatori indipendenti (uno per ogni dimensione)
        array<mt19937, dim> engines;
        for (size_t i = 0; i < dim; ++i)
            engines[i].seed(seeds[i]);

        // Creo un array di distribuzioni uniformi (una per ogni dimensione)
        array<uniform_real_distribution<double>, dim> distributions;
        for (size_t i = 0; i < dim; ++i)
        {
            auto bounds = this->domain.getBounds();
            distributions[i] = uniform_real_distribution<double>(bounds[i].first,
                                                                 bounds[i].second);
        }

        // Array di vettori: per ogni dimensione, un vettore di 'numbers' valori
        vector<Point<dim>> random_numbers;
        random_numbers.reserve(numbers);

        // Apro file di output
        std::ofstream outfile;
        if (typeid(domain) == typeid(Hypersphere<dim>))
        {
            outfile.open("hsphere_samples.dat");
        }
        else if (typeid(domain) == typeid(HyperCylinder<dim>))
        {
            outfile.open("cylinder_samples.dat");
        }
        else if (typeid(domain) == typeid(HyperRectangle<dim>))
        {
            outfile.open("hrectangle_samples.dat");
        }
        else
        {
            outfile.open("generic_samples.dat");
        }

        // Genero 'numbers' punti nel dominio, uno per riga
        for (int j = 0; j < numbers; ++j)
        {
            Point<dim> x;
            for (size_t i = 0; i < dim; ++i)
            {
                x[i] = distributions[i](engines[i]);
                outfile << x[i];
                if (i + 1 < dim)
                    outfile << ' ';
            }
            random_numbers.push_back(x);
            outfile << "\n";
        }

        return random_numbers;
    }

public:
    explicit Integrator(const IntegrationDomain<dim> &d) : domain(d) {}

    virtual ~Integrator() = default;
};

#endif // MONTECARLO_1_INTEGRATOR_HPP