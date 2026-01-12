//
// Created by Giacomo Merlo on 12/01/26.
//

#ifndef MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP
#define MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP

#include <random>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include "../domains/integration_domain.hpp"
#include "../geometry.hpp"

template <std::size_t dim>
struct VolumeEstimate {
    double volume = 0.0;     // estimated |D|
    double stderr = 0.0;     // standard error estimate
    double inside_ratio = 0.0; // p-hat
    std::size_t n_samples = 0;
};

template <std::size_t dim>
class VolumeEstimatorMC {
public:
    VolumeEstimate<dim>
    estimate(const IntegrationDomain<dim>& domain,
             std::mt19937& rng,
             std::size_t n_samples) const;
};

#include "VolumeEstimatorMC.tpp"

#endif //MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP