// GA.cpp
#include "GA.hpp"
#include "../rng/rng_factory.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace mc{
namespace optim{

    GA::GA(const GAConfig& config)
        : m_config(config),
          m_global_best(Solution::make_worst(OptimizationMode::MINIMIZE)),
          m_rng(mc::rng::make_engine(100)) // stream_id=100 for GA
    {}

    void GA::setObjectiveFunction(ObjectiveFunction func) {
        m_func = std::move(func);
        m_initialized = false;
    }

    void GA::setCallback(StepCallback cb) {
        m_callback = std::move(cb);
    }

    void GA::setBounds(const Coordinates& lower, const Coordinates& upper) {
        if (lower.size() != upper.size()) {
            throw std::invalid_argument("Lower and Upper bounds must have the same dimension.");
        }
        m_lower_bounds = lower;
        m_upper_bounds = upper;
        m_initialized = false;
    }

    void GA::setMode(OptimizationMode mode) {
        m_mode = mode;
        m_global_best = Solution::make_worst(m_mode);
        m_initialized = false;
    }

    bool GA::isBetterFitness(Real a, Real b) const {
        return (m_mode == OptimizationMode::MINIMIZE) ? (a < b) : (a > b);
    }

    void GA::evaluate(Individual& ind) {
        ind.fitness = m_func(ind.genome);
        // Global best is updated serially outside of parallel regions
    }

    void GA::enforceBounds(Coordinates& x) {
        for (size_t i = 0; i < x.size(); ++i) {
            if (x[i] < m_lower_bounds[i]) x[i] = m_lower_bounds[i];
            else if (x[i] > m_upper_bounds[i]) x[i] = m_upper_bounds[i];
        }
    }

    void GA::initialize() {
        if (!m_func) throw std::runtime_error("Objective function not set.");
        if (m_lower_bounds.empty()) throw std::runtime_error("Bounds not set.");
        if (m_config.population_size == 0) throw std::runtime_error("Population size must be > 0.");
        if (m_config.elitism_count >= m_config.population_size)
            throw std::runtime_error("Elitism count must be < population size.");

        m_population.clear();
        m_population.resize(m_config.population_size);

        const size_t dim = m_lower_bounds.size();
        std::uniform_real_distribution<Real> u01(0.0, 1.0);

        m_generation = 0;
        m_global_best = Solution::make_worst(m_mode);

        // Genome init stays serial (uses shared RNG)
        for (auto& ind : m_population) {
            ind.genome.resize(dim);
            for (size_t i = 0; i < dim; ++i) {
                Real span = m_upper_bounds[i] - m_lower_bounds[i];
                ind.genome[i] = m_lower_bounds[i] + u01(m_rng) * span;
            }
            ind.fitness = 0; // will be evaluated below
        }

        // Fitness evaluation can be parallel
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(m_population.size()); ++i) {
            evaluate(m_population[static_cast<size_t>(i)]);
        }

        // Serial update of global best to ensure deterministic tie-breaking
        for (const auto& ind : m_population) {
            Solution sol{ind.genome, ind.fitness};
            if (sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = sol;
            }
        }

        m_initialized = true;
    }

    const GA::Individual& GA::tournamentSelect() {
        std::uniform_int_distribution<size_t> pick(0, m_population.size() - 1);

        const Individual* best = nullptr;
        for (size_t i = 0; i < m_config.tournament_k; ++i) {
            const Individual& cand = m_population[pick(m_rng)];
            if (!best || isBetterFitness(cand.fitness, best->fitness)) {
                best = &cand;
            }
        }
        return *best;
    }

    void GA::crossoverUniform(const Coordinates& p1, const Coordinates& p2,
                              Coordinates& c1, Coordinates& c2) {
        std::uniform_real_distribution<Real> u01(0.0, 1.0);
        const size_t dim = p1.size();

        c1 = p1;
        c2 = p2;

        for (size_t i = 0; i < dim; ++i) {
            if (u01(m_rng) < 0.5) {
                std::swap(c1[i], c2[i]);
            }
        }
    }

    void GA::mutateGaussian(Coordinates& x) {
        std::uniform_real_distribution<Real> u01(0.0, 1.0);
        std::normal_distribution<Real> n01(0.0, 1.0);

        for (size_t i = 0; i < x.size(); ++i) {
            if (u01(m_rng) < m_config.mutation_rate) {
                Real span = m_upper_bounds[i] - m_lower_bounds[i];
                Real sigma = m_config.mutation_sigma * span;
                x[i] += n01(m_rng) * sigma;
            }
        }
        enforceBounds(x);
    }

    void GA::step() {
        if (!m_initialized) initialize();

        std::sort(m_population.begin(), m_population.end(),
                  [&](const Individual& a, const Individual& b) {
                      return isBetterFitness(a.fitness, b.fitness);
                  });

        std::vector<Individual> next;
        next.reserve(m_population.size());

        for (size_t i = 0; i < m_config.elitism_count; ++i) {
            next.push_back(m_population[i]);
        }

        std::uniform_real_distribution<Real> u01(0.0, 1.0);

        // Selection + variation stays serial (uses shared RNG)
        while (next.size() < m_population.size()) {
            const Individual& p1 = tournamentSelect();
            const Individual& p2 = tournamentSelect();

            Individual c1, c2;
            c1.genome = p1.genome;
            c2.genome = p2.genome;

            if (u01(m_rng) < m_config.crossover_rate) {
                crossoverUniform(p1.genome, p2.genome, c1.genome, c2.genome);
            }

            mutateGaussian(c1.genome);
            mutateGaussian(c2.genome);

            c1.fitness = 0; // evaluated below
            next.push_back(std::move(c1));

            if (next.size() < m_population.size()) {
                c2.fitness = 0; // evaluated below
                next.push_back(std::move(c2));
            }
        }

        // Fitness evaluation can be parallel (skip elites)
        const size_t start = m_config.elitism_count;
        #pragma omp parallel for
        for (int i = static_cast<int>(start); i < static_cast<int>(next.size()); ++i) {
            evaluate(next[static_cast<size_t>(i)]);
        }

        m_population = std::move(next);
        ++m_generation;

        // Serial update of global best after population evolves
        for (const auto& ind : m_population) {
            Solution sol{ind.genome, ind.fitness};
            if (sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = sol;
            }
        }
    }

    Solution GA::optimize() {
        initialize();

        for (size_t gen = 0; gen < m_config.max_generations; ++gen) {
            step();

            if (m_callback) {
                m_callback(m_global_best, gen);
            }
        }
        return m_global_best;
    }

    Solution GA::getBestSolution() const {
        return m_global_best;
    }

} //namespace mc
} //namespace optim