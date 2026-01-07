// GA.hpp
#pragma once
#include "optimizer.hpp"
#include <vector>
#include <random>
#include <stdexcept>

namespace optimizers {

    /**
     * @brief Configuration parameters specific to the Genetic Algorithm.
     * Real-coded GA for continuous optimization.
     */
    struct GAConfig {
        size_t population_size = 80;
        size_t max_generations = 200;

        // Selection
        size_t tournament_k = 3;              // Tournament size

        // Variation operators
        Real crossover_rate = 0.9;            // Probability to crossover
        Real mutation_rate  = 0.1;            // Probability to mutate each gene
        Real mutation_sigma = 0.1;            // Stddev for Gaussian mutation (scaled by span)

        // Elitism
        size_t elitism_count = 1;             // Keep best N individuals each generation
    };

    /**
     * @brief Genetic Algorithm implementation.
     * Suitable for continuous, non-linear optimization.
     */
    class GA : public Optimizer {
    public:
        struct Individual {
            Coordinates genome;   // Candidate solution (same as params)
            Real fitness;         // Objective value
        };

        explicit GA(const GAConfig& config = GAConfig{});

        // Optimizer interface
        void setObjectiveFunction(ObjectiveFunction func) override;
        void setBounds(const Coordinates& lower, const Coordinates& upper) override;
        void setMode(OptimizationMode mode) override;
        void setCallback(StepCallback cb) override;

        Solution optimize() override;
        void step() override;
        [[nodiscard]] Solution getBestSolution() const override;

        [[nodiscard]] const std::vector<Individual>& getPopulation() const {
            return m_population;
        }

    private:
        void initialize();
        void evaluate(Individual& ind);
        void enforceBounds(Coordinates& x);

        // Selection + variation
        const Individual& tournamentSelect();
        void crossoverUniform(const Coordinates& p1, const Coordinates& p2,
                              Coordinates& c1, Coordinates& c2);
        void mutateGaussian(Coordinates& x);

        // Utils
        bool isBetterFitness(Real a, Real b) const;

        // Members
        GAConfig m_config;
        OptimizationMode m_mode = OptimizationMode::MINIMIZE;
        ObjectiveFunction m_func;

        Coordinates m_lower_bounds;
        Coordinates m_upper_bounds;

        std::vector<Individual> m_population;
        Solution m_global_best;

        std::mt19937 m_rng;
        bool m_initialized = false;

        StepCallback m_callback;
        size_t m_generation = 0;
    };

} // namespace optimizers