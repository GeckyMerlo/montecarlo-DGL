#include "PSO.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdint>

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

        std::uniform_real_distribution<Real> dist(0.0, 1.0);

        for (auto& p : m_swarm) {
            p.position.resize(dim);
            p.velocity.resize(dim);

            // Random initialization within bounds
            for(size_t i=0; i<dim; ++i) {
                Real span = m_upper_bounds[i] - m_lower_bounds[i];
                p.position[i] = m_lower_bounds[i] + dist(m_rng) * span;
                p.velocity[i] = (dist(m_rng) - 0.5) * span * 0.1; // Small initial random velocity
            }

            // Evaluate initial fitness
            p.current_value = m_func(p.position);

            // Initial pBest is the current position
            p.best_position = p.position;
            p.best_value = p.current_value;

            // Update Global Best if this particle is the new champion
            Solution current_sol = {p.position, p.current_value};
            if (current_sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = current_sol;
            }
        }
        m_initialized = true;
    }

    void PSO::step() {
        if (!m_initialized) initialize();

        std::uniform_real_distribution<Real> r_dist(0.0, 1.0);
        size_t dim = m_lower_bounds.size();

        for (auto& p : m_swarm) {
            for (size_t i = 0; i < dim; ++i) {
                Real r1 = r_dist(m_rng); // Random cognitive factor
                Real r2 = r_dist(m_rng); // Random social factor

                // --- Velocity Update Equation ---
                // v[t+1] = w*v[t] + c1*r1*(pBest - x) + c2*r2*(gBest - x)
                Real cognitive_comp = m_config.cognitive_coeff * r1 * (p.best_position[i] - p.position[i]);
                Real social_comp    = m_config.social_coeff    * r2 * (m_global_best.params[i] - p.position[i]);

                p.velocity[i] = (m_config.inertia_weight * p.velocity[i]) + cognitive_comp + social_comp;

                // --- Position Update ---
                p.position[i] += p.velocity[i];
            }

            // Keep particle inside the box
            enforceBounds(p);

            // --- Evaluation ---
            Real new_val = m_func(p.position);
            p.current_value = new_val;

            // --- pBest Update ---
            Solution new_sol = {p.position, new_val};
            Solution pbest_sol = {p.best_position, p.best_value};

            if (new_sol.isBetterThan(pbest_sol, m_mode)) {
                p.best_value = new_val;
                p.best_position = p.position;
            }

            // --- gBest Update ---
            if (new_sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = new_sol;
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