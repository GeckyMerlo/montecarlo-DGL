#include "PSO.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <omp.h>
#include "../RngManager.hpp"

// Global seed defined in main.cpp
extern uint32_t GLOBAL_SEED;

namespace optimizers {

    PSO::PSO(const PSOConfig& config)
        : m_config(config),
          // Initialize global best to worst possible value based on default mode
          m_global_best(Solution::make_worst(OptimizationMode::MINIMIZE)),
          // Seed the random generator
          m_rng(GLOBAL_SEED)
    {}

    void PSO::setObjectiveFunction(ObjectiveFunction func) {
        m_func = std::move(func);
    }

    // --- Implementation of setCallback ---
    void PSO::setCallback(StepCallback cb) {
        m_callback = cb;
    }

    void PSO::setBounds(const Coordinates& lower, const Coordinates& upper) {
        if (lower.size() != upper.size()) {
            throw std::invalid_argument("Lower and Upper bounds must have the same dimension.");
        }
        m_lower_bounds = lower;
        m_upper_bounds = upper;
    }

    void PSO::setMode(OptimizationMode mode) {
        m_mode = mode;
        // Reset global best because "worst" definition changes
        m_global_best = Solution::make_worst(m_mode);
    }

    void PSO::initialize() {
        if (!m_func) throw std::runtime_error("Objective function not set.");
        if (m_lower_bounds.empty()) throw std::runtime_error("Bounds not set.");

        m_swarm.resize(m_config.population_size);
        size_t dim = m_lower_bounds.size();

        // Parallel initialization of particles
        #pragma omp parallel for schedule(dynamic)
        for (int p_idx = 0; p_idx < static_cast<int>(m_swarm.size()); ++p_idx) {
            auto& p = m_swarm[p_idx];
            
            // Local RNG for thread-safe random number generation
            RngManager local_rng(GLOBAL_SEED + p_idx);
            auto local_gen = local_rng.make_rng(0);
            std::uniform_real_distribution<Real> dist(0.0, 1.0);

            p.position.resize(dim);
            p.velocity.resize(dim);

            // Random initialization within bounds
            for(size_t i=0; i<dim; ++i) {
                Real span = m_upper_bounds[i] - m_lower_bounds[i];
                p.position[i] = m_lower_bounds[i] + dist(local_gen) * span;
                p.velocity[i] = (dist(local_gen) - 0.5) * span * 0.1;
            }

            // Evaluate initial fitness
            p.current_value = m_func(p.position);

            // Initial pBest is the current position
            p.best_position = p.position;
            p.best_value = p.current_value;
        }
        
        // Find global best (sequential, small overhead)
        for (const auto& p : m_swarm) {
            Solution current_sol = {p.position, p.current_value};
            if (current_sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = current_sol;
            }
        }
        
        m_initialized = true;
    }

    void PSO::step() {
        if (!m_initialized) initialize();

        size_t dim = m_lower_bounds.size();
        
        // Store old global best for update detection
        Solution old_global_best = m_global_best;

        // Parallel particle update
        #pragma omp parallel for schedule(dynamic)
        for (int p_idx = 0; p_idx < static_cast<int>(m_swarm.size()); ++p_idx) {
            auto& p = m_swarm[p_idx];
            
            // Local RNG for thread-safe random generation
            RngManager local_rng(GLOBAL_SEED + p_idx + 1000);
            auto local_gen = local_rng.make_rng(0);
            std::uniform_real_distribution<Real> r_dist(0.0, 1.0);

            for (size_t i = 0; i < dim; ++i) {
                Real r1 = r_dist(local_gen);
                Real r2 = r_dist(local_gen);

                // Velocity update
                Real cognitive_comp = m_config.cognitive_coeff * r1 * (p.best_position[i] - p.position[i]);
                Real social_comp    = m_config.social_coeff    * r2 * (old_global_best.params[i] - p.position[i]);

                p.velocity[i] = (m_config.inertia_weight * p.velocity[i]) + cognitive_comp + social_comp;

                // Position update
                p.position[i] += p.velocity[i];
            }

            // Keep particle inside bounds
            enforceBounds(p);

            // Evaluate new position
            Real new_val = m_func(p.position);
            p.current_value = new_val;

            // Update pBest
            Solution new_sol = {p.position, new_val};
            Solution pbest_sol = {p.best_position, p.best_value};

            if (new_sol.isBetterThan(pbest_sol, m_mode)) {
                p.best_value = new_val;
                p.best_position = p.position;
            }

            // Update gBest (with critical section for thread safety)
            if (new_sol.isBetterThan(m_global_best, m_mode)) {
                #pragma omp critical
                {
                    if (new_sol.isBetterThan(m_global_best, m_mode)) {
                        m_global_best = new_sol;
                    }
                }
            }
        }
    }

    void PSO::enforceBounds(Particle& p) {
        for (size_t i = 0; i < p.position.size(); ++i) {
            if (p.position[i] < m_lower_bounds[i]) {
                p.position[i] = m_lower_bounds[i];
                p.velocity[i] *= -0.5; // Bounce back with damping
            } else if (p.position[i] > m_upper_bounds[i]) {
                p.position[i] = m_upper_bounds[i];
                p.velocity[i] *= -0.5; // Bounce back
            }
        }
    }

    // --- Optimize loop with Callback ---
    Solution PSO::optimize() {
        initialize();
        for(size_t i=0; i<m_config.max_iterations; ++i) {
            step();

            // If the user provided a callback, execute it passing current best and iteration
            if (m_callback) {
                m_callback(m_global_best, i);
            }
        }
        return m_global_best;
    }

    Solution PSO::getBestSolution() const {
        return m_global_best;
    }
}