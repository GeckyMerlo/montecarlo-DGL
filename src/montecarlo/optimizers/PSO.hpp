#pragma once
#include "optimizer.hpp"
#include <vector>
#include <random>

namespace mc{
namespace optim{

    struct PSOConfig {
        size_t population_size = 50;
        size_t max_iterations = 100;

        Real inertia_weight = 0.7;
        Real cognitive_coeff = 1.5;
        Real social_coeff = 1.5;
    };

    class PSO : public Optimizer {
    public:
        struct Particle {
            Coordinates position;
            Coordinates velocity;

            Coordinates best_position;
            Real best_value;
            Real current_value;
        };

        explicit PSO(const PSOConfig& config = PSOConfig{});

        void setObjectiveFunction(ObjectiveFunction func) override;
        void setBounds(const Coordinates& lower, const Coordinates& upper) override;
        void setMode(OptimizationMode mode) override;
        void setCallback(StepCallback cb) override;

        Solution optimize() override;
        void step() override;
        [[nodiscard]] Solution getBestSolution() const override;

        [[nodiscard]] const std::vector<Particle>& getParticles() const {
            return m_swarm;
        }

    private:
        void initialize();
        void enforceBounds(Particle& p);

        PSOConfig m_config;
        OptimizationMode m_mode = OptimizationMode::MINIMIZE;
        ObjectiveFunction m_func;

        Coordinates m_lower_bounds;
        Coordinates m_upper_bounds;

        std::vector<Particle> m_swarm;
        Solution m_global_best;

        bool m_initialized = false;

        // Iteration counter used to derive unique RNG stream IDs per step
        size_t m_current_iter = 0;

        StepCallback m_callback;
    };
} //namespace mc
} //namespace optim