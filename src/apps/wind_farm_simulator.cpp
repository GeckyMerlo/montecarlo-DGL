/**
 * @file wind_farm_simulator.cpp
 * @brief Wind Farm Layout Optimization using Hybrid Monte Carlo + PSO
 * 
 * @details This application optimizes the placement of N wind turbines in a 2D domain
 * using Particle Swarm Optimization (PSO). The objective function uses Monte Carlo
 * simulation to estimate average energy production under varying wind conditions.
 * 
 * **Problem Setup:**
 * - 2D domain: rectangular wind farm area
 * - Optimization variables: (x, y) positions for each turbine
 * - Objective: maximize average energy production considering wake effects
 * 
 * **Physical Model:**
 * - Wind speed sampled from Weibull distribution (realistic wind modeling)
 * - Energy production: Power ∝ v³ (cubic law for wind turbines)
 * - Wake interference: massive penalty if turbines are too close (< 50m)
 * 
 * **Outputs:**
 * - Console: optimization progress and final layout
 * - results.dat: optimal turbine coordinates
 * - plot_script.gp: Gnuplot visualization script
 * 
 * @see PSO, RngManager
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <limits>
#include <functional>
#include <chrono>
#include <omp.h>

// Include montecarlo library headers
#include "../montecarlo/rng/rng_global.hpp"
#include "../montecarlo/rng/rng_factory.hpp"
#include "../montecarlo/optimizers/PSO.hpp"
#include "../montecarlo/optimizers/types.hpp"

using namespace optimizers;

// =============================================================================
// CONFIGURATION CONSTANTS
// =============================================================================

/// Number of wind turbines to optimize
constexpr size_t NUM_TURBINES = 15;

/// Wind farm domain size (meters)
constexpr double FARM_WIDTH = 1000.0;   // x range: [0, FARM_WIDTH]
constexpr double FARM_HEIGHT = 1000.0;  // y range: [0, FARM_HEIGHT]

/// Minimum distance between turbines (meters) - wake/collision threshold
constexpr double MIN_TURBINE_DISTANCE = 50.0;

/// Penalty coefficient for turbines too close together
constexpr double PROXIMITY_PENALTY = 1e8;

/// Number of Monte Carlo samples for wind speed simulation
constexpr size_t MC_WIND_SAMPLES = 500;

/// Weibull distribution parameters (typical for wind)
/// Shape parameter k (around 2 for Rayleigh-like distribution)
constexpr double WEIBULL_SHAPE = 2.0;
/// Scale parameter λ (mean wind speed ~7-10 m/s typical)
constexpr double WEIBULL_SCALE = 8.0;

/// Turbine rated power coefficient (simplified)
constexpr double POWER_COEFFICIENT = 0.4;  // Betz limit ~0.593
/// Air density (kg/m³)
constexpr double AIR_DENSITY = 1.225;
/// Rotor swept area (m²) - assuming 50m diameter rotor
constexpr double ROTOR_AREA = M_PI * 25.0 * 25.0;

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

/**
 * @brief Calculate Euclidean distance between two turbines in 2D
 */
inline double turbineDistance(double x1, double y1, double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief Parse turbine positions from PSO coordinates vector
 * @param coords Flat vector of [x0, y0, x1, y1, ..., xN-1, yN-1]
 * @param x Output vector of x coordinates
 * @param y Output vector of y coordinates
 */
inline void extractTurbinePositions(const Coordinates& coords, 
                                     std::vector<double>& x, 
                                     std::vector<double>& y) {
    size_t n = coords.size() / 2;
    x.resize(n);
    y.resize(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = coords[2 * i];
        y[i] = coords[2 * i + 1];
    }
}

/**
 * @brief Check all pairwise distances and compute proximity penalty
 * @param x X-coordinates of turbines
 * @param y Y-coordinates of turbines
 * @return Total penalty (0 if all distances >= MIN_TURBINE_DISTANCE)
 */
double computeProximityPenalty(const std::vector<double>& x, 
                                const std::vector<double>& y) {
    double penalty = 0.0;
    size_t n = x.size();
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dist = turbineDistance(x[i], y[i], x[j], y[j]);
            if (dist < MIN_TURBINE_DISTANCE) {
                // Quadratic penalty: closer turbines = higher penalty
                double violation = MIN_TURBINE_DISTANCE - dist;
                penalty += PROXIMITY_PENALTY * violation * violation;
            }
        }
    }
    return penalty;
}

/**
 * @brief Calculate power output for a given wind speed
 * Power = 0.5 * ρ * A * Cp * v³
 */
inline double windPower(double wind_speed) {
    return 0.5 * AIR_DENSITY * ROTOR_AREA * POWER_COEFFICIENT * 
           std::pow(wind_speed, 3.0);
}

/**
 * @brief Sample wind speed from Weibull distribution
 * @param rng Random number generator
 * @return Wind speed in m/s
 */
inline double sampleWeibullWindSpeed(std::mt19937& rng) {
    std::weibull_distribution<double> weibull(WEIBULL_SHAPE, WEIBULL_SCALE);
    return weibull(rng);
}

/**
 * @brief Apply simplified wake effect based on distance
 * Turbines downwind experience reduced wind speed
 * @param base_speed Original wind speed
 * @param distance Distance from upwind turbine
 * @return Reduced wind speed
 */
inline double applyWakeEffect(double base_speed, double distance) {
    // Jensen wake model (simplified)
    // Wake recovery increases with distance
    if (distance < MIN_TURBINE_DISTANCE) {
        return base_speed * 0.3; // Severe wake loss
    }
    double wake_decay = 0.04; // Wake expansion coefficient
    double recovery = 1.0 - 0.3 / (1.0 + wake_decay * distance / 50.0);
    return base_speed * std::min(1.0, recovery);
}

// =============================================================================
// OBJECTIVE FUNCTION
// =============================================================================

/**
 * @brief Objective function for wind farm optimization
 * 
 * Uses Monte Carlo simulation to estimate average energy production:
 * 1. Sample multiple wind speeds from Weibull distribution
 * 2. For each sample, compute total farm power output
 * 3. Average over all samples
 * 4. Apply proximity penalty for constraint violations
 * 
 * @param coords Flat vector of turbine positions [x0, y0, x1, y1, ...]
 * @param thread_seed Seed offset for thread-safe RNG
 * @return Negative average power (negative because PSO minimizes)
 */
double windFarmObjective(const Coordinates& coords, std::uint64_t thread_seed) {
    std::vector<double> x, y;
    extractTurbinePositions(coords, x, y);
    size_t n = x.size();
    
    // Check for proximity constraint violations
    double penalty = computeProximityPenalty(x, y);
    if (penalty > 0) {
        // Return large negative value (since we're maximizing, PSO minimizes -energy)
        return penalty;
    }
    
    // Monte Carlo simulation for wind speed
    auto rng = mc::make_engine_with_seed(static_cast<std::uint32_t>(thread_seed), 0);
    
    double total_energy = 0.0;
    
    for (size_t sample = 0; sample < MC_WIND_SAMPLES; ++sample) {
        // Sample a wind speed for this Monte Carlo iteration
        double wind_speed = sampleWeibullWindSpeed(rng);
        
        // Random wind direction (0 to 2π)
        std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
        double wind_direction = angle_dist(rng);
        
        // Calculate total power for all turbines considering wake effects
        double sample_power = 0.0;
        
        for (size_t i = 0; i < n; ++i) {
            double turbine_wind = wind_speed;
            
            // Check wake effects from all other turbines
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                
                // Check if turbine j is upwind of turbine i
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dist = std::sqrt(dx * dx + dy * dy);
                
                // Project onto wind direction
                double wind_proj = dx * std::cos(wind_direction) + 
                                   dy * std::sin(wind_direction);
                
                // If turbine j is upwind (positive projection), apply wake
                if (wind_proj > 0 && dist < 500.0) {
                    // Wake effect decreases with distance
                    turbine_wind = applyWakeEffect(turbine_wind, dist);
                }
            }
            
            // Add power from this turbine
            sample_power += windPower(turbine_wind);
        }
        
        total_energy += sample_power;
    }
    
    // Average energy over all Monte Carlo samples
    double avg_energy = total_energy / static_cast<double>(MC_WIND_SAMPLES);
    
    // Return negative because PSO minimizes (we want to maximize energy)
    return -avg_energy;
}

// =============================================================================
// OUTPUT AND VISUALIZATION
// =============================================================================

/**
 * @brief Write optimized turbine positions to data file
 */
void writeResultsFile(const std::string& filename, 
                       const std::vector<double>& x,
                       const std::vector<double>& y) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create " << filename << std::endl;
        return;
    }
    
    out << "# Wind Farm Turbine Positions\n";
    out << "# x_coord  y_coord  turbine_id\n";
    
    for (size_t i = 0; i < x.size(); ++i) {
        out << std::fixed << std::setprecision(2) 
            << x[i] << "  " << y[i] << "  " << (i + 1) << "\n";
    }
    
    out.close();
    std::cout << "[INFO] Results written to " << filename << std::endl;
}

/**
 * @brief Generate Gnuplot script for wind farm visualization
 */
void writePlotScript(const std::string& filename,
                      const std::string& data_file) {
    std::ofstream gp(filename);
    if (!gp.is_open()) {
        std::cerr << "Error: Cannot create " << filename << std::endl;
        return;
    }
    
    gp << "# Gnuplot script for Wind Farm Visualization\n";
    gp << "# Generated by wind_farm_simulator.cpp\n\n";
    
    // Output settings
    gp << "set terminal pngcairo size 1200,1000 enhanced font 'Arial,12'\n";
    gp << "set output 'wind_farm_layout.png'\n\n";
    
    // Title and labels
    gp << "set title 'Wind Farm Optimization Layout' font 'Arial,16'\n";
    gp << "set xlabel 'X Position (m)' font 'Arial,12'\n";
    gp << "set ylabel 'Y Position (m)' font 'Arial,12'\n\n";
    
    // Grid and styling
    gp << "set grid\n";
    gp << "set key outside right top\n\n";
    
    // Set axis ranges
    gp << "set xrange [-50:" << (FARM_WIDTH + 50) << "]\n";
    gp << "set yrange [-50:" << (FARM_HEIGHT + 50) << "]\n";
    gp << "set size ratio 1\n\n";
    
    // Draw farm boundary
    gp << "# Draw farm boundary\n";
    gp << "set object 1 rect from 0,0 to " << FARM_WIDTH << "," << FARM_HEIGHT;
    gp << " fs empty border lc rgb 'black' lw 2\n\n";
    
    // Define triangle marker for turbines (pointing up - north wind direction)
    gp << "# Use triangles to represent turbines (wind direction indicator)\n";
    gp << "# Point type 9 = filled triangle up in many Gnuplot versions\n\n";
    
    // Plot command
    gp << "plot '" << data_file << "' using 1:2 with points \\\n";
    gp << "     pt 9 ps 3 lc rgb '#2E86AB' title 'Wind Turbines', \\\n";
    gp << "     '" << data_file << "' using 1:2:3 with labels \\\n";
    gp << "     offset 0,1.5 font 'Arial,10' notitle\n\n";
    
    // Interactive version (for terminal display)
    gp << "# For interactive viewing, uncomment below:\n";
    gp << "# set terminal qt\n";
    gp << "# replot\n";
    gp << "# pause -1 'Press Enter to close'\n";
    
    gp.close();
    std::cout << "[INFO] Gnuplot script written to " << filename << std::endl;
}

/**
 * @brief Print optimization summary to console
 */
void printSummary(const Solution& best, 
                   const std::vector<double>& x,
                   const std::vector<double>& y) {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║          WIND FARM OPTIMIZATION - FINAL RESULTS              ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════╣\n";
    
    std::cout << "║ Optimized Turbine Positions:                                 ║\n";
    std::cout << "║                                                              ║\n";
    
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << "║   Turbine " << std::setw(2) << (i + 1) 
                  << ":  x = " << std::fixed << std::setprecision(2) << std::setw(8) << x[i]
                  << " m,  y = " << std::setw(8) << y[i] << " m               ║\n";
    }
    
    std::cout << "║                                                              ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════╣\n";
    
    // Compute average power (negate because we stored -energy)
    double avg_power = -best.value;
    double avg_power_kw = avg_power / 1000.0;
    double avg_power_mw = avg_power / 1e6;
    
    std::cout << "║ Total Average Power Output: " << std::setw(12) << std::setprecision(2) 
              << avg_power_kw << " kW                  ║\n";
    std::cout << "║                             " << std::setw(12) << std::setprecision(4)
              << avg_power_mw << " MW                  ║\n";
    
    // Compute minimum pairwise distance
    double min_dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = i + 1; j < x.size(); ++j) {
            double d = turbineDistance(x[i], y[i], x[j], y[j]);
            min_dist = std::min(min_dist, d);
        }
    }
    
    std::cout << "║ Minimum Turbine Separation: " << std::setw(12) << std::setprecision(2)
              << min_dist << " m                   ║\n";
    std::cout << "║ Required Minimum Distance:  " << std::setw(12) << MIN_TURBINE_DISTANCE 
              << " m                   ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║       WIND FARM LAYOUT OPTIMIZATION                          ║\n";
    std::cout << "║       Hybrid Monte Carlo + Particle Swarm Optimization       ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    // Set global seed for reproducibility
    mc::set_global_seed(42u);
    std::cout << "[INFO] Global RNG seed set to 42\n";
    
    // Configuration info
    std::cout << "[INFO] Number of turbines: " << NUM_TURBINES << "\n";
    std::cout << "[INFO] Farm size: " << FARM_WIDTH << "m x " << FARM_HEIGHT << "m\n";
    std::cout << "[INFO] Min turbine distance: " << MIN_TURBINE_DISTANCE << "m\n";
    std::cout << "[INFO] Monte Carlo samples per evaluation: " << MC_WIND_SAMPLES << "\n";
    std::cout << "[INFO] Weibull wind parameters: k=" << WEIBULL_SHAPE 
              << ", λ=" << WEIBULL_SCALE << " m/s\n\n";
    
    // Define search space bounds
    // Each turbine has (x, y), so total dimensions = 2 * NUM_TURBINES
    size_t total_dims = 2 * NUM_TURBINES;
    
    Coordinates lower_bounds(total_dims);
    Coordinates upper_bounds(total_dims);
    
    for (size_t i = 0; i < NUM_TURBINES; ++i) {
        // x coordinate bounds
        lower_bounds[2 * i] = 0.0;
        upper_bounds[2 * i] = FARM_WIDTH;
        // y coordinate bounds
        lower_bounds[2 * i + 1] = 0.0;
        upper_bounds[2 * i + 1] = FARM_HEIGHT;
    }
    
    // Configure PSO
    PSOConfig config;
    config.population_size = 60;    // Number of particles
    config.max_iterations = 150;    // Maximum iterations
    config.inertia_weight = 0.6;    // Velocity inertia
    config.cognitive_coeff = 1.8;   // Personal best attraction
    config.social_coeff = 2.0;      // Global best attraction
    
    std::cout << "[INFO] PSO Configuration:\n";
    std::cout << "       - Population size: " << config.population_size << "\n";
    std::cout << "       - Max iterations: " << config.max_iterations << "\n";
    std::cout << "       - Inertia weight: " << config.inertia_weight << "\n";
    std::cout << "       - Cognitive coefficient: " << config.cognitive_coeff << "\n";
    std::cout << "       - Social coefficient: " << config.social_coeff << "\n\n";
    
    // Create PSO optimizer
    PSO optimizer(config);
    
    // Set bounds
    optimizer.setBounds(lower_bounds, upper_bounds);
    
    // Set optimization mode (minimize because objective returns -energy)
    optimizer.setMode(OptimizationMode::MINIMIZE);
    
    // Create thread-safe objective function wrapper
    // Each evaluation gets a unique seed based on parameters hash
    auto objectiveWrapper = [](const Coordinates& coords) -> Real {
        // Create a simple hash from coordinates for reproducible seeding
        std::uint64_t hash = 0;
        for (size_t i = 0; i < coords.size(); ++i) {
            hash ^= std::hash<double>{}(coords[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        // Add thread ID for thread safety
        #ifdef _OPENMP
        hash += static_cast<std::uint64_t>(omp_get_thread_num()) * 1000000;
        #endif
        
        return windFarmObjective(coords, hash);
    };
    
    optimizer.setObjectiveFunction(objectiveWrapper);
    
    // Set callback for progress reporting
    size_t callback_counter = 0;
    optimizer.setCallback([&callback_counter](const Solution& best, size_t iteration) {
        if (iteration % 10 == 0 || iteration == 1) {
            double power_mw = -best.value / 1e6;
            std::cout << "[ITER " << std::setw(4) << iteration << "] "
                      << "Best power: " << std::fixed << std::setprecision(4) 
                      << power_mw << " MW\n";
        }
        callback_counter++;
    });
    
    // Run optimization
    std::cout << "[INFO] Starting optimization...\n\n";
    
    auto start_time = std::chrono::high_resolution_clock::now();
    Solution best = optimizer.optimize();
    auto end_time = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\n[INFO] Optimization completed in " 
              << duration.count() / 1000.0 << " seconds\n";
    
    // Extract final positions
    std::vector<double> opt_x, opt_y;
    extractTurbinePositions(best.params, opt_x, opt_y);
    
    // Print summary
    printSummary(best, opt_x, opt_y);
    
    // Write output files
    std::string results_file = "results.dat";
    std::string plot_script = "plot_script.gp";
    
    writeResultsFile(results_file, opt_x, opt_y);
    writePlotScript(plot_script, results_file);
    
    std::cout << "\n[INFO] To visualize results, run:\n";
    std::cout << "       gnuplot " << plot_script << "\n";
    std::cout << "       Then open 'wind_farm_layout.png'\n\n";
    
    return 0;
}
