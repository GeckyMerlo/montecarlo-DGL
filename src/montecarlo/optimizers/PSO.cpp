/**
 * @file PSO.cpp
 * @brief Particle Swarm Optimization implementation
 */

#include "PSO.hpp"
#include "../rng/rng_factory.hpp"
#include <stdexcept>
#include <cstdint>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace mc{
namespace optim{

    PSO::PSO(const PSOConfig& config)
        : m_config(config),
          m_global_best(Solution::make_worst(OptimizationMode::MINIMIZE))
    {}

    void PSO::setObjectiveFunction(ObjectiveFunction func) {
        m_func = std::move(func);
        m_initialized = false;
    }

    void PSO::setCallback(StepCallback cb) {
        m_callback = std::move(cb);
    }

    void PSO::setBounds(const Coordinates& lower, const Coordinates& upper) {
        if (lower.size() != upper.size()) {
            throw std::invalid_argument("Lower and Upper bounds must have the same dimension.");
        }
        m_lower_bounds = lower;
        m_upper_bounds = upper;
        m_initialized = false;
    }

    void PSO::setMode(OptimizationMode mode) {
        m_mode = mode;
        m_global_best = Solution::make_worst(m_mode);
        m_initialized = false;
    }

    void PSO::initialize() {
        if (!m_func) throw std::runtime_error("Objective function not set.");
        if (m_lower_bounds.empty()) throw std::runtime_error("Bounds not set.");

        m_swarm.resize(m_config.population_size);
        const size_t dim = m_lower_bounds.size();
        m_current_iter = 0; // Reset iteration counter for deterministic seeding

        // Init positions/velocities and evaluate
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int p_idx = 0; p_idx < static_cast<int>(m_swarm.size()); ++p_idx) {
            auto& p = m_swarm[static_cast<size_t>(p_idx)];

            // Deterministic RNG decoupled from thread scheduling.
            // Unique stream per particle.
            auto local_gen = mc::rng::make_engine(1000ULL + static_cast<std::uint64_t>(p_idx));
            std::uniform_real_distribution<Real> dist(0.0, 1.0);

            p.position.resize(dim);
            p.velocity.resize(dim);

            for (size_t i = 0; i < dim; ++i) {
                const Real span = m_upper_bounds[i] - m_lower_bounds[i];
                p.position[i] = m_lower_bounds[i] + dist(local_gen) * span;
                p.velocity[i] = (dist(local_gen) - 0.5) * span * 0.1;
            }

            p.current_value = m_func(p.position);
            p.best_position = p.position;
            p.best_value = p.current_value;
        }

        // Find global best (small, keep it simple)
        m_global_best = Solution::make_worst(m_mode);
        for (const auto& p : m_swarm) {
            Solution current_sol{p.position, p.current_value};
            if (current_sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = current_sol;
            }
        }

        m_initialized = true;
    }

    void PSO::step() {
        if (!m_initialized) initialize();

        const size_t dim = m_lower_bounds.size();
        const Solution old_global_best = m_global_best;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int p_idx = 0; p_idx < static_cast<int>(m_swarm.size()); ++p_idx) {
            auto& p = m_swarm[static_cast<size_t>(p_idx)];

            // Deterministic RNG: mix iteration count into stream id so values change each step.
            std::uint64_t stream_id = 2000ULL +
                                      (static_cast<std::uint64_t>(m_current_iter) * m_config.population_size) +
                                      static_cast<std::uint64_t>(p_idx);
            auto local_gen = mc::rng::make_engine(stream_id);
            std::uniform_real_distribution<Real> r_dist(0.0, 1.0);

            for (size_t i = 0; i < dim; ++i) {
                const Real r1 = r_dist(local_gen);
                const Real r2 = r_dist(local_gen);

                const Real cognitive_comp =
                    m_config.cognitive_coeff * r1 * (p.best_position[i] - p.position[i]);
                const Real social_comp =
                    m_config.social_coeff * r2 * (old_global_best.params[i] - p.position[i]);

                p.velocity[i] =
                    (m_config.inertia_weight * p.velocity[i]) + cognitive_comp + social_comp;

                p.position[i] += p.velocity[i];
            }

            enforceBounds(p);

            const Real new_val = m_func(p.position);
            p.current_value = new_val;

            Solution new_sol{p.position, new_val};
            Solution pbest_sol{p.best_position, p.best_value};

            if (new_sol.isBetterThan(pbest_sol, m_mode)) {
                p.best_value = new_val;
                p.best_position = p.position;
            }

            // Global best update moved to a serial post-processing section for determinism
        }

        // Serial update of global best to ensure deterministic tie-breaking
        for (const auto& p : m_swarm) {
            Solution sol{p.best_position, p.best_value};
            if (sol.isBetterThan(m_global_best, m_mode)) {
                m_global_best = sol;
            }
        }

        // Advance iteration counter for next step's RNG seeding
        ++m_current_iter;
    }

    void PSO::enforceBounds(Particle& p) {
        for (size_t i = 0; i < p.position.size(); ++i) {
            if (p.position[i] < m_lower_bounds[i]) {
                p.position[i] = m_lower_bounds[i];
                p.velocity[i] *= -0.5;
            } else if (p.position[i] > m_upper_bounds[i]) {
                p.position[i] = m_upper_bounds[i];
                p.velocity[i] *= -0.5;
            }
        }
    }

    Solution PSO::optimize() {
        initialize();
        for (size_t i = 0; i < m_config.max_iterations; ++i) {
            step();
            if (m_callback) {
                m_callback(m_global_best, i);
            }
        }
        return m_global_best;
    }

    Solution PSO::getBestSolution() const {
        return m_global_best;
    }

} //namespace mc
} //namespace optim