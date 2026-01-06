#pragma once
#include "optimizer.hpp"
#include <vector>
#include <random>

namespace optimizers {

    /**
     * @brief Configuration parameters specific to the PSO algorithm.
     */
    struct PSOConfig {
        size_t population_size = 50;  // Number of particles
        size_t max_iterations = 100;  // Stopping criterion

        // Coefficients for the velocity update equation:
        Real inertia_weight = 0.7;   // w: How much to keep previous velocity
        Real cognitive_coeff = 1.5;  // c1: Attraction to personal best
        Real social_coeff = 1.5;     // c2: Attraction to global best
    };

    /**
     * @brief Particle Swarm Optimization implementation.
     * Suitable for continuous, non-linear optimization problems.
     */
    class PSO : public Optimizer {
    public:
        // Moved Particle struct to public to allow external visualization
        struct Particle {
            Coordinates position;
            Coordinates velocity;

            Coordinates best_position; // Personal Best (pBest) location
            Real best_value;           // Personal Best (pBest) value
            Real current_value;        // Current value
        };

        // Constructor takes configuration, defaults provided in struct
        explicit PSO(const PSOConfig& config = PSOConfig{});

        // Override Interface methods
        void setObjectiveFunction(ObjectiveFunction func) override;
        void setBounds(const Coordinates& lower, const Coordinates& upper) override;
        void setMode(OptimizationMode mode) override;
        void setCallback(StepCallback cb) override;

        Solution optimize() override;
        void step() override;
        [[nodiscard]] Solution getBestSolution() const override;

        // Getter to access the swarm from outside (for plotting/debugging)
        [[nodiscard]] const std::vector<Particle>& getParticles() const {
            return m_swarm;
        }

    private:
        // Helper to initialize the swarm randomly within bounds
        void initialize();

        // Check if particles flew out of bounds and correct them
        void enforceBounds(Particle& p);

        // Member variables
        PSOConfig m_config;
        OptimizationMode m_mode = OptimizationMode::MINIMIZE;
        ObjectiveFunction m_func;

        Coordinates m_lower_bounds;
        Coordinates m_upper_bounds;

        std::vector<Particle> m_swarm;
        Solution m_global_best; // Global Best (gBest)

        // Random Number Generation (Mersenne Twister)
        std::mt19937 m_rng;
        bool m_initialized = false;

        StepCallback m_callback;
    };
}