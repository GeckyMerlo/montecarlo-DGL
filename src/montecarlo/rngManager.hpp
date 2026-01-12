//
// Created by Giacomo Merlo on 12/01/26.
//

#ifndef MONTECARLO_DGL_RNGMANAGER_HPP
#define MONTECARLO_DGL_RNGMANAGER_HPP

#include <random>
#include <cstdint>

class RngManager {
public:
    explicit RngManager(uint64_t seed) : master_seed(seed) {}

    std::mt19937 make_rng(int thread_id, int run_id = 0) const {
        std::seed_seq seq{
            static_cast<uint32_t>(master_seed),
            static_cast<uint32_t>(thread_id),
            static_cast<uint32_t>(run_id)
        };
        return std::mt19937(seq);
    }

private:
    uint64_t master_seed;
};

#endif //MONTECARLO_DGL_RNGMANAGER_HPP