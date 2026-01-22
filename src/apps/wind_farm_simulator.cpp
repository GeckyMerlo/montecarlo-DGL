/**
 * @file wind_farm_simulator.cpp
 * @brief Wind Farm Layout Optimization using Hybrid MH-Monte Carlo + PSO/GA (comparison)
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
#include <cstdint>
#include <array>
#include <omp.h>

// --- Montecarlo library headers ---
#include "../montecarlo/rng/rng_global.hpp"
#include "../montecarlo/rng/rng_factory.hpp"

#include "../montecarlo/optimizers/PSO.hpp"
#include "../montecarlo/optimizers/GA.hpp"
#include "../montecarlo/optimizers/types.hpp"

// IMPORTANT: include proposals BEFORE MHintegrator.hpp
#include "../montecarlo/proposals/proposal.hpp"
#include "../montecarlo/proposals/gaussianProposal.hpp"

#include "../montecarlo/integrators/MHintegrator.hpp"

#include "../montecarlo/domains/integration_domain.hpp"
#include "../montecarlo/domains/hyperrectangle.hpp"
#include "../montecarlo/geometry.hpp"

// =============================================================================
// CONFIGURATION CONSTANTS
// =============================================================================

constexpr size_t NUM_TURBINES = 15;

constexpr double FARM_WIDTH  = 1000.0;
constexpr double FARM_HEIGHT = 1000.0;

constexpr double MIN_TURBINE_DISTANCE = 50.0;
constexpr double PROXIMITY_PENALTY    = 1e8;

// Weibull wind parameters
constexpr double WEIBULL_SHAPE = 2.0;
constexpr double WEIBULL_SCALE = 8.0;

// Turbine model constants
constexpr double POWER_COEFFICIENT = 0.4;
constexpr double AIR_DENSITY       = 1.225;
constexpr double ROTOR_AREA        = M_PI * 25.0 * 25.0;

// --- MH integration configuration ---
constexpr int    MH_SAMPLES         = 1500;
constexpr size_t MH_BURN_IN         = 400;
constexpr size_t MH_THINNING        = 2;
constexpr size_t MH_SAMPLES_FOR_VOL = 2000;

// Proposal sigmas (in physical coordinates v, theta)
constexpr double PROPOSAL_SIGMA_V     = 2.5;   // m/s
constexpr double PROPOSAL_SIGMA_THETA = 0.6;   // rad

// Wind integration bounds (physical)
constexpr double WIND_SPEED_MIN  = 0.0;
constexpr double WIND_SPEED_MAX  = 40.0;
constexpr double WIND_THETA_MIN  = 0.0;
constexpr double WIND_THETA_MAX  = 2.0 * M_PI;

// Centers for mapping from centered domain coords to physical coords
constexpr double WIND_SPEED_CENTER = 0.5 * (WIND_SPEED_MIN + WIND_SPEED_MAX); // 20
constexpr double WIND_THETA_CENTER = 0.5 * (WIND_THETA_MIN + WIND_THETA_MAX); // pi

// =============================================================================
// UTILITY: WIND + GEOMETRY
// =============================================================================

inline double turbineDistance(double x1, double y1, double x2, double y2) {
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

inline void extractTurbinePositions(const mc::optim::Coordinates& coords,
                                    std::vector<double>& x,
                                    std::vector<double>& y) {
    const size_t n = coords.size() / 2;
    x.resize(n);
    y.resize(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = coords[2 * i];
        y[i] = coords[2 * i + 1];
    }
}

double computeProximityPenalty(const std::vector<double>& x,
                               const std::vector<double>& y) {
    double penalty = 0.0;
    const size_t n = x.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double dist = turbineDistance(x[i], y[i], x[j], y[j]);
            if (dist < MIN_TURBINE_DISTANCE) {
                const double violation = MIN_TURBINE_DISTANCE - dist;
                penalty += PROXIMITY_PENALTY * violation * violation;
            }
        }
    }
    return penalty;
}

inline double windPower(double wind_speed) {
    return 0.5 * AIR_DENSITY * ROTOR_AREA * POWER_COEFFICIENT *
           std::pow(wind_speed, 3.0);
}

inline double applyWakeEffect(double base_speed, double distance) {
    if (distance < MIN_TURBINE_DISTANCE) {
        return base_speed * 0.3;
    }
    const double wake_decay = 0.04;
    const double recovery = 1.0 - 0.3 / (1.0 + wake_decay * distance / 50.0);
    return base_speed * std::min(1.0, recovery);
}

// Weibull pdf (unnormalized OK for MH)
inline double weibull_pdf(double v) {
    if (v < 0.0) return 0.0;
    const double k = WEIBULL_SHAPE;
    const double l = WEIBULL_SCALE;
    const double a = (k / l);
    const double x = v / l;
    return a * std::pow(x, k - 1.0) * std::exp(-std::pow(x, k));
}

// Map centered-domain point w -> physical (v,theta)
inline void mapToPhysicalWind(const mc::geom::Point<2>& w, double& v, double& theta) {
    v     = w[0] + WIND_SPEED_CENTER;
    theta = w[1] + WIND_THETA_CENTER;
}

// =============================================================================
// DOMAIN BUILDER (HyperRectangle takes dims, centered at origin)
// =============================================================================

static mc::domains::IntegrationDomain<2>* buildWindDomainOwned() {
    // Domain extents (dims) for centered box:
    // w0 in [-dv/2, +dv/2], w1 in [-dtheta/2, +dtheta/2]
    std::array<double, 2> dims{
        WIND_SPEED_MAX - WIND_SPEED_MIN,
        WIND_THETA_MAX - WIND_THETA_MIN
    };
    return new mc::domains::HyperRectangle<2>(dims);
}

// =============================================================================
// POWER FOR GIVEN (v,theta) AND LAYOUT
// =============================================================================

static inline double farmPowerGivenWind(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       double wind_speed,
                                       double wind_direction) {
    const size_t n = x.size();
    double sample_power = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double turbine_wind = wind_speed;

        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;

            const double dx   = x[i] - x[j];
            const double dy   = y[i] - y[j];
            const double dist = std::sqrt(dx * dx + dy * dy);

            const double wind_proj = dx * std::cos(wind_direction) +
                                     dy * std::sin(wind_direction);

            if (wind_proj > 0.0 && dist < 500.0) {
                turbine_wind = applyWakeEffect(turbine_wind, dist);
            }
        }

        sample_power += windPower(turbine_wind);
    }

    return sample_power;
}

// =============================================================================
// MH-BASED EXPECTATION ESTIMATION
// =============================================================================

static double estimateAveragePowerMH(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    std::uint32_t seed32) {
    using Point2 = mc::geom::Point<2>;
    using Integrator2 = mc::integrators::MHMontecarloIntegrator<2>;

    mc::domains::IntegrationDomain<2>* domain = buildWindDomainOwned();
    Integrator2 mh(*domain);

    // Target density p(w) = p(v(w),theta(w)) with theta uniform -> proportional to weibull(v)
    auto p = [](const Point2& w) -> double {
        double v, theta;
        mapToPhysicalWind(w, v, theta);

        (void)theta; // uniform -> constant
        if (v < WIND_SPEED_MIN || v > WIND_SPEED_MAX) return 0.0;
        // theta should be in [0,2pi] by construction (domain mapping)
        return weibull_pdf(v);
    };

    // Initial point in centered coordinates: choose physical (v=WEIBULL_SCALE, theta=pi)
    Point2 x0;
    x0[0] = WEIBULL_SCALE - WIND_SPEED_CENTER; // shift to centered coords
    x0[1] = M_PI - WIND_THETA_CENTER;          // = 0

    // deviation_ parameter: keep 1.0 unless your MH uses it internally
    const double mh_deviation_param = 1.0;

    mh.setConfig(MH_BURN_IN,
                 MH_THINNING,
                 MH_SAMPLES_FOR_VOL,
                 mh_deviation_param,
                 p,
                 x0);

    // GaussianProposal samples in centered coordinates w
    // mean in centered coords for (v=WEIBULL_SCALE, theta=pi)
    std::vector<double> mu{
        WEIBULL_SCALE - WIND_SPEED_CENTER,
        M_PI - WIND_THETA_CENTER
    };
    std::vector<double> sigma{
        PROPOSAL_SIGMA_V,
        PROPOSAL_SIGMA_THETA
    };

    mc::proposals::GaussianProposal<2> proposal(*domain, mu, sigma);

    auto f = [&](const Point2& w) -> double {
        double v, theta;
        mapToPhysicalWind(w, v, theta);
        return farmPowerGivenWind(x, y, v, theta);
    };

    const double avg_power = mh.integrate(f, MH_SAMPLES, proposal, seed32);

    delete domain;
    return avg_power;
}

// =============================================================================
// OBJECTIVE FUNCTION (PSO/GA unchanged)
// =============================================================================

inline std::uint64_t hashCoordsForSeed(const mc::optim::Coordinates& coords) {
    std::uint64_t hash = 0;
    for (size_t i = 0; i < coords.size(); ++i) {
        hash ^= std::hash<double>{}(coords[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
#ifdef _OPENMP
    hash += static_cast<std::uint64_t>(omp_get_thread_num()) * 1000000ULL;
#endif
    return hash;
}

double windFarmObjective(const mc::optim::Coordinates& coords, std::uint64_t thread_seed) {
    std::vector<double> x, y;
    extractTurbinePositions(coords, x, y);

    const double penalty = computeProximityPenalty(x, y);
    if (penalty > 0.0) return penalty;

    const std::uint32_t seed32 = static_cast<std::uint32_t>(thread_seed);
    const double avg_power = estimateAveragePowerMH(x, y, seed32);

    return -avg_power;
}

// =============================================================================
// OUTPUT AND VISUALIZATION
// =============================================================================

void writeResultsFile(const std::string& filename,
                      const std::vector<double>& x,
                      const std::vector<double>& y) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create " << filename << "\n";
        return;
    }

    out << "# Wind Farm Turbine Positions\n";
    out << "# x_coord  y_coord  turbine_id\n";
    for (size_t i = 0; i < x.size(); ++i) {
        out << std::fixed << std::setprecision(2)
            << x[i] << "  " << y[i] << "  " << (i + 1) << "\n";
    }
    std::cout << "[INFO] Results written to " << filename << "\n";
}

void writePlotScript(const std::string& filename,
                     const std::string& data_file,
                     const std::string& output_png) {
    std::ofstream gp(filename);
    if (!gp.is_open()) {
        std::cerr << "Error: Cannot create " << filename << "\n";
        return;
    }

    gp << "set terminal pngcairo size 1200,1000 enhanced font 'Arial,12'\n";
    gp << "set output '" << output_png << "'\n\n";
    gp << "set title 'Wind Farm Optimization Layout'\n";
    gp << "set xlabel 'X Position (m)'\n";
    gp << "set ylabel 'Y Position (m)'\n";
    gp << "set grid\n";
    gp << "set key outside right top\n";
    gp << "set xrange [-50:" << (FARM_WIDTH + 50) << "]\n";
    gp << "set yrange [-50:" << (FARM_HEIGHT + 50) << "]\n";
    gp << "set size ratio 1\n\n";
    gp << "set object 1 rect from 0,0 to " << FARM_WIDTH << "," << FARM_HEIGHT
       << " fs empty border lc rgb 'black' lw 2\n\n";

    gp << "plot '" << data_file << "' using 1:2 with points pt 9 ps 3 title 'Wind Turbines',\\\n"
          "     '" << data_file << "' using 1:2:3 with labels offset 0,1.5 notitle\n";

    std::cout << "[INFO] Gnuplot script written to " << filename << "\n";
}

void printSummary(const std::string& name,
                  const mc::optim::Solution& best,
                  const std::vector<double>& x,
                  const std::vector<double>& y,
                  double seconds) {
    const double avg_power_mw = (-best.value) / 1e6;

    double min_dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = i + 1; j < x.size(); ++j) {
            min_dist = std::min(min_dist, turbineDistance(x[i], y[i], x[j], y[j]));
        }
    }

    std::cout << "\n=== " << name << " SUMMARY ===\n";
    std::cout << "Time: " << std::fixed << std::setprecision(3) << seconds << " s\n";
    std::cout << "Best avg power: " << std::fixed << std::setprecision(6) << avg_power_mw << " MW\n";
    std::cout << "Min distance: " << std::fixed << std::setprecision(2) << min_dist << " m\n";
}

// =============================================================================
// BENCH HARNESS (PSO/GA API unchanged)
// =============================================================================

struct RunResult {
    std::string name;
    mc::optim::Solution best;
    std::vector<double> x, y;
    double seconds = 0.0;
};

template <typename OptimizerT>
RunResult runOptimizer(const std::string& name,
                       OptimizerT& optimizer,
                       const mc::optim::Coordinates& lower_bounds,
                       const mc::optim::Coordinates& upper_bounds,
                       const std::string& results_file,
                       const std::string& plot_script,
                       const std::string& output_png)
{
    optimizer.setBounds(lower_bounds, upper_bounds);
    optimizer.setMode(mc::optim::OptimizationMode::MINIMIZE);

    optimizer.setObjectiveFunction([](const mc::optim::Coordinates& coords) -> mc::optim::Real {
        return windFarmObjective(coords, hashCoordsForSeed(coords));
    });

    optimizer.setCallback([name](const mc::optim::Solution& best, size_t it) {
        if (it % 10 == 0 || it == 1) {
            const double power_mw = -best.value / 1e6;
            std::cout << "[" << name << " | ITER " << std::setw(4) << it << "] "
                      << "Best power: " << std::fixed << std::setprecision(6)
                      << power_mw << " MW\n";
        }
    });

    std::cout << "\n[INFO] Starting " << name << " optimization...\n\n";

    const auto t0 = std::chrono::high_resolution_clock::now();
    mc::optim::Solution best = optimizer.optimize();
    const auto t1 = std::chrono::high_resolution_clock::now();

    const double seconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;

    std::vector<double> x, y;
    extractTurbinePositions(best.params, x, y);

    writeResultsFile(results_file, x, y);
    writePlotScript(plot_script, results_file, output_png);

    return RunResult{ name, best, x, y, seconds };
}

// =============================================================================
// MAIN
// =============================================================================

int main() {
    std::cout << "\n=== WIND FARM: MH + PSO vs GA ===\n";

    mc::rng::set_global_seed(42u);
    std::cout << "[INFO] Global RNG seed set to 42\n";

    // Layout bounds
    const size_t total_dims = 2 * NUM_TURBINES;
    mc::optim::Coordinates lower_bounds(total_dims);
    mc::optim::Coordinates upper_bounds(total_dims);

    for (size_t i = 0; i < NUM_TURBINES; ++i) {
        lower_bounds[2 * i]     = 0.0;
        upper_bounds[2 * i]     = FARM_WIDTH;
        lower_bounds[2 * i + 1] = 0.0;
        upper_bounds[2 * i + 1] = FARM_HEIGHT;
    }

    // --- PSO (UNCHANGED API) ---
    mc::optim::PSOConfig pso_cfg;
    pso_cfg.population_size = 60;
    pso_cfg.max_iterations  = 150;
    pso_cfg.inertia_weight  = 0.6;
    pso_cfg.cognitive_coeff = 1.8;
    pso_cfg.social_coeff    = 2.0;
    mc::optim::PSO pso(pso_cfg);

    // --- GA (UNCHANGED API) ---
    mc::optim::GAConfig ga_cfg;
    ga_cfg.population_size = 80;
    ga_cfg.max_generations = 200;
    ga_cfg.tournament_k    = 3;
    ga_cfg.crossover_rate  = 0.9;
    ga_cfg.mutation_rate   = 0.10;
    ga_cfg.mutation_sigma  = 25.0;  // meters
    ga_cfg.elitism_count   = 2;
    mc::optim::GA ga(ga_cfg);

    // Run both
    RunResult pso_res = runOptimizer("PSO", pso, lower_bounds, upper_bounds,
                                     "results_pso.dat", "plot_pso.gp", "wind_farm_layout_pso.png");

    RunResult ga_res  = runOptimizer("GA", ga, lower_bounds, upper_bounds,
                                     "results_ga.dat", "plot_ga.gp", "wind_farm_layout_ga.png");

    // Summaries
    printSummary("PSO", pso_res.best, pso_res.x, pso_res.y, pso_res.seconds);
    printSummary("GA",  ga_res.best,  ga_res.x,  ga_res.y,  ga_res.seconds);

    const double pso_mw = -pso_res.best.value / 1e6;
    const double ga_mw  = -ga_res.best.value  / 1e6;

    std::cout << "\n=== COMPARISON ===\n";
    std::cout << "PSO: " << std::fixed << std::setprecision(6) << pso_mw
              << " MW  | time " << std::setprecision(2) << pso_res.seconds << " s\n";
    std::cout << "GA : " << std::fixed << std::setprecision(6) << ga_mw
              << " MW  | time " << std::setprecision(2) << ga_res.seconds << " s\n";
    std::cout << "Winner (power): " << ((ga_mw > pso_mw) ? "GA" : "PSO") << "\n";

    std::cout << "\n[INFO] Plots:\n";
    std::cout << "  gnuplot plot_pso.gp  -> wind_farm_layout_pso.png\n";
    std::cout << "  gnuplot plot_ga.gp   -> wind_farm_layout_ga.png\n\n";

    return 0;
}