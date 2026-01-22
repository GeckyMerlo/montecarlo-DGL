// VolumeEstimatorMC.hpp
/**
 * @file VolumeEstimatorMC.hpp
 * @brief Hit-or-Miss Monte Carlo volume estimation
 * @author Giacomo Merlo
 * @date 12/01/26
 *
 * Estimates the volume of complex domains using acceptance-rejection sampling:
 * \f[
 * V_{\text{est}} = V_{\text{box}} \cdot \frac{N_{\text{hits}}}{N_{\text{total}}}
 * \f]
 */

#ifndef MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP
#define MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <stdexcept>

#include "../domains/integration_domain.hpp"
#include "../geometry.hpp"

template <std::size_t dim>
struct VolumeEstimate {
    double volume = 0.0;        ///< Estimated |D| (domain volume)
    double stderr = 0.0;        ///< Standard error estimate
    double inside_ratio = 0.0;  ///< Fraction of samples inside (p-hat)
    std::size_t n_samples = 0;  ///< Total samples used
};

template <std::size_t dim>
class VolumeEstimatorMC {
public:
    VolumeEstimate<dim>
    estimate(const IntegrationDomain<dim>& domain,
             std::uint32_t seed,
             std::size_t n_samples) const;
};

#include "VolumeEstimatorMC.tpp"

#endif // MONTECARLO_DGL_VOLUMEESTIMATORMC_HPP