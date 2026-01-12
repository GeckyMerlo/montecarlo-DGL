//
// Created by domenico on 12/3/25.
//

#include "apps/benchmarks.hpp"
#include <montecarlo/utils/muParserXInterface.hpp>
#include <cmath>
#include <cstdint>

using namespace MuParserInterface;

namespace opt = optimizers;

// --- Global Configuration ---
const std::vector<size_t> n_samples_vector = {10'000, 50'000, 100'000, 500'000, 1'000'000};
unsigned int n_threads;

// --- Utility Functions ---

void saveResults(const std::string &filename, const std::vector<results> &results, const std::string &function_expr) {
    std::ofstream outfile;
    outfile.open(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to create results file " << filename << std::endl;
        return;
    }

    // Header: Function description and columns
    outfile << "Function: " << function_expr << "\n";
    outfile << "Number of points\tIntegration Result\tDuration (ms)\n";

    for (const auto &result : results) {
        outfile << result.n_samples << "\t"
                << result.integration_result << "\t"
                << result.duration << "\n";
    }
    outfile.close();
}

/**
 * GENERIC execution loop.
 * Runs the integration for increasing sample sizes, measures time, plots results, and saves to file.
 */
template <size_t dim, typename Func>
void executeBenchmark(const std::string& title,
                      const std::string& filename,
                      MontecarloIntegrator<dim>& integrator,
                      const IntegrationDomain<dim>& domain, // Needed for plotting
                      Func&& f,
                      bool useGnuplot,
                      const std::string& rawDataFile,
                      const std::string& functionExpr)      // Function string for display/save
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << title << ":" << std::endl;
    std::cout << "Integrating Function: " << functionExpr << std::endl << std::endl;

    std::vector<results> testResults;

    for (size_t n_i : n_samples_vector) {
        auto startTimer = std::chrono::high_resolution_clock::now();

        // 1. Integration (Data points are written to rawDataFile by the integrator)
        double result = integrator.integrate(f, n_i);

        auto endTimer = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer - startTimer);

        auto startTimer2 = std::chrono::high_resolution_clock::now();

        // 1. Integration (Data points are written to rawDataFile by the integrator)
        double result2 = integrator.integrate_importance(f, n_i, UniformProposal<dim>(domain), 12345);

        auto endTimer2 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer2 - startTimer2);

        // Store results
        results newLine;
        newLine.n_samples = n_i;
        newLine.integration_result = std::to_string(result);
        newLine.duration = std::to_string(duration.count());
        testResults.push_back(newLine);

        // Print to console
        std::cout << std::setw(12) << n_i << " Normal Samples"
                  << std::setw(18) << std::fixed << std::setprecision(6) << result
                  << std::setw(18) << duration.count() << "ms" << std::endl;

        std::cout << std::setw(12) << n_i << " Important Samples"
                  << std::setw(18) << std::fixed << std::setprecision(6) << result2
                  << std::setw(18) << duration2.count() << "ms" << std::endl;

        // 2. Plotting (Inside the loop to show progress for each sample size)
        if (useGnuplot) {
            // Plot Domain Geometry (Inside/Outside points)
            createGnuplotScript(rawDataFile, domain, n_i);

            // Plot Function Value (f(x))
            createFunctionGnuplotScript(rawDataFile, domain, f, n_i);
        }
    }

    std::cout << "\nSaved txt file: " << filename << "\n";
    saveResults(filename, testResults, functionExpr);
}

// --- Domain-Specific Benchmark Wrappers ---

template <typename Func>
void runCircleBenchmark(Func f, const std::string& modeLabel, bool useGnuplot, const std::string& funcStr) {
    Hypersphere<2> circle(5.0);
    MontecarloIntegrator<2> integrator(circle);
    std::string title = "2D Circle Integration (Radius 5) [" + modeLabel + "]";
    std::string filename = "resultsCircle_" + modeLabel + ".txt";
    std::string dataFile = "hsphere_samples.dat"; // File name must match what Integrator writes

    executeBenchmark(title, filename, integrator, circle, f, useGnuplot, dataFile, funcStr);
}

template <typename Func>
void runSphereBenchmark(Func f, const std::string& modeLabel, bool useGnuplot, const std::string& funcStr) {
    double radius = 10.0;
    Hypersphere<4> sphere(radius);
    MontecarloIntegrator<4> integrator(sphere);
    std::string title = "4D Hypersphere Integration [" + modeLabel + "]";
    std::string filename = "resultsSphere4D_" + modeLabel + ".txt";
    std::string dataFile = "hsphere_samples.dat";

    executeBenchmark(title, filename, integrator, sphere, f, useGnuplot, dataFile, funcStr);
}

template <typename Func>
void runRectBenchmark(Func f, const std::string& modeLabel, bool useGnuplot, const std::string& funcStr) {
    std::array<double, 4> sides = {10.0, 5.0, 10.0, 5.0};
    HyperRectangle<4> rectangle(sides);
    MontecarloIntegrator<4> integrator(rectangle);
    std::string title = "4D HyperRectangle Integration [" + modeLabel + "]";
    std::string filename = "resultsRectangle4D_" + modeLabel + ".txt";
    std::string dataFile = "hrectangle_samples.dat";

    executeBenchmark(title, filename, integrator, rectangle, f, useGnuplot, dataFile, funcStr);
}

template <typename Func>
void runCylinderBenchmark(Func f, const std::string& modeLabel, bool useGnuplot, const std::string& funcStr) {
    double radius = 5.0;
    double height = 10.0;
    HyperCylinder<4> cylinder(radius, height);
    MontecarloIntegrator<4> integrator(cylinder);
    std::string title = "4D HyperCylinder Integration [" + modeLabel + "]";
    std::string filename = "resultsCylinder4D_" + modeLabel + ".txt";
    std::string dataFile = "cylinder_samples.dat";

    executeBenchmark(title, filename, integrator, cylinder, f, useGnuplot, dataFile, funcStr);
}


// --- Specific Implementations (Hardcoded vs Parser) ---

// 1. HARDCODED (C++ Lambda)
// Manually defining the function string description for the output file/console.

void circleIntegration(bool useGnuplot) {
    auto f = [](const Point<2> &x) { return x[0] * x[0] - x[1] * x[1]; };
    std::string funcStr = "x[0]^2 - x[1]^2";
    runCircleBenchmark(f, "Hardcoded", useGnuplot, funcStr);
}

void sphereIntegration(bool useGnuplot) {
    auto f = [](const Point<4> &x) { return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]; };
    std::string funcStr = "x[0]^2 + x[1]^2 + x[2]^2 + x[3]^2";
    runSphereBenchmark(f, "Hardcoded", useGnuplot, funcStr);
}

void rectangularIntegration(bool useGnuplot) {
    auto f = [](const Point<4> &x) { return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]; };
    std::string funcStr = "x[0]^2 + x[1]^2 + x[2]^2 + x[3]^2";
    runRectBenchmark(f, "Hardcoded", useGnuplot, funcStr);
}

void cylinderIntegration(bool useGnuplot) {
    auto f = [](const Point<4> &x) { return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]; };
    std::string funcStr = "x[0]^2 + x[1]^2 + x[2]^2 + x[3]^2";
    runCylinderBenchmark(f, "Hardcoded", useGnuplot, funcStr);
}

// 2. PARSER (muParserX)
// Using the 'expr' string directly read from file.

void circleIntegrationParser(const std::string& expr, bool useGnuplot) {
    muParserXInterface<2, Point<2>> parser(expr);
    runCircleBenchmark(parser, "Parser", useGnuplot, expr);
}

void sphereIntegrationParser(const std::string& expr, bool useGnuplot) {
    muParserXInterface<4, Point<4>> parser(expr);
    runSphereBenchmark(parser, "Parser", useGnuplot, expr);
}

void rectangularIntegrationParser(const std::string& expr, bool useGnuplot) {
    muParserXInterface<4, Point<4>> parser(expr);
    runRectBenchmark(parser, "Parser", useGnuplot, expr);
}

void cylinderIntegrationParser(const std::string& expr, bool useGnuplot) {
    muParserXInterface<4, Point<4>> parser(expr);
    runCylinderBenchmark(parser, "Parser", useGnuplot, expr);
}


// --- Main Entry Points ---

void runBenchmarks(bool useGnuplot) {
    n_threads = std::thread::hardware_concurrency();
    if (n_threads == 0) n_threads = 16;

    circleIntegration(useGnuplot);
    sphereIntegration(useGnuplot);
    rectangularIntegration(useGnuplot);
    cylinderIntegration(useGnuplot);

}

void runBenchmarksMH(bool useGnuplot) {
    n_threads = std::thread::hardware_concurrency();
    if (n_threads == 0) n_threads = 16;

    Hypersphere<2> domain(10.0);

    std::function<double(const geom::Point<2>&)> indicator =
    [&domain](const geom::Point<2>& x) -> double {
        return domain.isInside(x) ? 1.0 : 0.0;
    };

    const double deviation = 0.15;
    const std::size_t burn_in = 20'000;
    const std::size_t thinning = 10;
    const std::size_t n_samples = 1'000'000;
    const std::size_t n_samples_volume = 200'000;

    geom::Point<2> x0;

    std::function<double(const geom::Point<2>&)> f = [](const geom::Point<2>& x) -> double {
        return x[0]*x[0] + x[1]*x[1];
    };

    MontecarloIntegrator<2> integrator(domain);

    std::cout << "Running Benchmarks" << std::endl;
    std::cout << "Metropolis Hastings result: " << integrator.integrate_with_mh(f, indicator, x0, deviation, GLOBAL_SEED, burn_in, n_samples, thinning, n_samples_volume)  << std::endl;
    std::cout << "Montecarlo result: " << integrator.integrate(f, n_samples) << std::endl;
    std::cout << "Exact result: " << 5000 * M_PI << std::endl;
}

void runBenchmarks(const std::string& expression, bool useGnuplot) {
    n_threads = std::thread::hardware_concurrency();
    if (n_threads == 0) n_threads = 16;

    // The expression is already printed in executeBenchmark, but printing it here once is fine too
    // std::cout << "Running benchmarks with expression: " << expression << std::endl;

    circleIntegrationParser(expression, useGnuplot);
    sphereIntegrationParser(expression, useGnuplot);
    rectangularIntegrationParser(expression, useGnuplot);
    cylinderIntegrationParser(expression, useGnuplot);
}

// --- Optimization Helper Functions ---

/**
 * @brief Test 1: Sphere Function
 * Objective: Minimize f(x, y) = x^2 + y^2
 * Global Minimum: 0 at [0, 0]
 */
void runSphereTest(opt::PSO& pso, const opt::Coordinates& lower, const opt::Coordinates& upper) {
    std::cout << "Optimization Problem: Minimize Sphere Function in 2D" << std::endl;
    std::cout << "Search Space: [-10, 10] per dimension" << std::endl;
    std::cout << "Running optimizer..." << std::endl;

    // Define the objective function: Sphere Function f(x) = sum(x_i^2)
    opt::ObjectiveFunction sphere_function = [](const opt::Coordinates& coords) {
        opt::Real sum = 0.0;
        for (auto val : coords) {
            sum += val * val;
        }
        return sum;
    };

    // Set up the optimizer
    pso.setBounds(lower, upper);
    pso.setObjectiveFunction(sphere_function);
    pso.setMode(opt::OptimizationMode::MINIMIZE);

    // Set a callback to print progress every 10 iterations
    pso.setCallback([](const opt::Solution& current_best, size_t iteration) {
        if (iteration == 0 || iteration % 10 == 0) {
            std::cout << "[Sphere Test | Step " << std::setw(3) << iteration << "] "
                      << "Best Value: " << std::scientific << std::setprecision(5)
                      << current_best.value << std::defaultfloat << std::endl;
        }
    });

    // Execute
    auto start = std::chrono::high_resolution_clock::now();
    opt::Solution best_sol = pso.optimize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    // Report
    std::cout << "\nOptimization Completed in " << duration.count() << " ms." << std::endl;
    std::cout << "Best Value Found: " << std::fixed << std::setprecision(10) << best_sol.value << std::endl;
    std::cout << "Best Position: [ ";
    for (auto val : best_sol.params) {
        std::cout << std::fixed << std::setprecision(5) << val << " ";
    }
    std::cout << "]" << std::endl;
}

/**
 * @brief Test 2: Boundary Constraint Test
 * Objective: Minimize f(x, y) = x + y
 * This function is a constant inclined plane. There is no local minimum inside
 * the domain (constant gradient). Particles must push to the extreme lower limit.
 * Expected Minimum: -20.0 at [-10, -10]
 */
void runBoundaryTest(opt::PSO& pso, const opt::Coordinates& lower, const opt::Coordinates& upper) {
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "TEST 2: Boundary Constraint Test (Linear Plane)" << std::endl;
    std::cout << "Objective: f(x,y) = x + y (Minimization)" << std::endl;
    std::cout << "Search Space: [-10, 10] per dimension" << std::endl;
    std::cout << "Expected Result: -20.0 at [-10.0, -10.0]" << std::endl;
    std::cout << "Running optimizer..." << std::endl;

    // Define the linear function
    opt::ObjectiveFunction plane_function = [](const opt::Coordinates& coords) {
        opt::Real sum = 0.0;
        for (auto val : coords) {
            sum += val;
        }
        return sum;
    };

    // Set up the optimizer (reusing the instance)
    pso.setBounds(lower, upper);
    pso.setObjectiveFunction(plane_function);
    // Mode is already MINIMIZE, but good practice to ensure
    pso.setMode(opt::OptimizationMode::MINIMIZE);

    // Callback to print progress
    pso.setCallback([](const opt::Solution& current_best, size_t iteration) {
        if (iteration % 20 == 0) {
            std::cout << "[Boundary Test | Step " << std::setw(3) << iteration << "] "
                      << "Val: " << std::fixed << std::setprecision(4) << current_best.value << std::endl;
        }
    });

    // Execute
    auto start = std::chrono::high_resolution_clock::now();
    opt::Solution best_sol = pso.optimize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    // Report
    std::cout << "Optimization Completed in " << duration.count() << " ms." << std::endl;
    std::cout << "Best Value Found: " << std::fixed << std::setprecision(5) << best_sol.value << std::endl;
    std::cout << "Best Position: [ ";
    for (auto val : best_sol.params) {
        std::cout << std::fixed << std::setprecision(5) << val << " ";
    }
    std::cout << "]" << std::endl;

    // Verification
    if (std::abs(best_sol.value - (-20.0)) < 1e-3) {
        std::cout << ">> SUCCESS: Boundary minimum found correctly!" << std::endl;
    } else {
        std::cout << ">> WARNING: Did not reach the exact boundary." << std::endl;
    }
}

/**
 * @brief Test 3: High-Dimensional Rastrigin Function
 * Objective: Minimize f(x) = 10n + sum(x_i^2 - 10cos(2*pi*x_i))
 * Domain: [-5.12, 5.12]
 * Global Minimum: 0.0 at x = [0, 0, ..., 0]
 *
 * Why it is hard:
 * This function creates a grid of local minima. In high dimensions (e.g., 10D),
 * simplistic optimizers get stuck in a local valley instead of finding the global 0.
 * A successful run requires a good balance of exploration and exploitation.
 */
void runRastriginTest(opt::PSO& pso, int dim) {
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "TEST 3: High-Dimensional Stress Test (Rastrigin Function)" << std::endl;
    std::cout << "Dimension: " << dim << "D" << std::endl;
    std::cout << "Search Space: [-5.12, 5.12] per dimension" << std::endl;
    std::cout << "Goal: Find global minimum 0.0 (avoiding local traps)" << std::endl;

    // 1. Define the Rastrigin function
    opt::ObjectiveFunction rastrigin_func = [dim](const opt::Coordinates& coords) {
        double sum = 0.0;
        double A = 10.0;
        // M_PI is standard in <cmath>, if missing use 3.14159265358979323846
        double pi = 3.14159265358979323846;

        for (auto x : coords) {
            sum += (x * x) - (A * std::cos(2.0 * pi * x));
        }
        return A * dim + sum;
    };

    // 2. Define Bounds for N dimensions
    opt::Coordinates lower(dim, -5.12);
    opt::Coordinates upper(dim, 5.12);

    // 3. Configure PSO specifically for a harder problem
    // We need more particles and more time to explore 10 dimensions
    opt::PSOConfig hard_config;
    hard_config.population_size = 500; // Increased from 50 - with population 100 only finds local minima
    hard_config.max_iterations = 1000; // Increased from 100
    hard_config.inertia_weight = 0.729; // Classic "constriction factor" value
    hard_config.cognitive_coeff = 1.49;
    hard_config.social_coeff = 1.49;

    // We create a new local instance to not mess up the previous config
    opt::PSO local_pso(hard_config);

    local_pso.setBounds(lower, upper);
    local_pso.setObjectiveFunction(rastrigin_func);
    local_pso.setMode(opt::OptimizationMode::MINIMIZE);

    // Minimal callback to show we are alive
    local_pso.setCallback([](const opt::Solution& sol, size_t i) {
        if (i % 100 == 0) {
            std::cout << "[Rastrigin " << i << "] Best: "
                      << std::scientific << std::setprecision(4) << sol.value
                      << std::defaultfloat << std::endl;
        }
    });

    std::cout << "Running optimizer (this might take longer)..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    opt::Solution best_sol = local_pso.optimize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Optimization Completed in " << duration.count() << " ms." << std::endl;
    std::cout << "Best Value Found: " << std::fixed << std::setprecision(5) << best_sol.value << std::endl;

    // Validation: It's very hard to get exactly 0.0 in 10D without huge resources.
    // Anything below 1.0 is considered a "good" result for a basic PSO.
    // Anything close to 0.0 is excellent.
    if (best_sol.value < 1e-2) {
        std::cout << ">> SUCCESS: Global minimum found!" << std::endl;
    } else if (best_sol.value < 5.0) {
        std::cout << ">> ACCEPTABLE: Found a good local minimum, but not global." << std::endl;
    } else {
        std::cout << ">> FAIL: Stuck in a high local minimum." << std::endl;
    }
}

void runVisualPSOBenchmark() {
    std::cout << "===========================================" << std::endl;
    std::cout << "   Visual PSO Benchmark (2D Animation)" << std::endl;
    std::cout << "===========================================" << std::endl;

    // 1. Configuration
    opt::PSOConfig config;
    config.population_size = 60; // 30 particles
    config.max_iterations = 100;  // 50 frames for animation
    opt::PSO pso(config);

    // 2. Objective Function (Rastrigin 2D)
    // Global minimum is at (0,0)
    auto rastrigin = [](const opt::Coordinates& x) {
        double A = 10.0;
        double sum = 0.0;
        // M_PI is typically defined in <cmath>, ensure it is available or use 3.14...
        double pi = 3.14159265358979323846;
        for (double val : x) sum += val*val - A*std::cos(2*pi*val);
        return 2*A + sum;
    };

    pso.setObjectiveFunction(rastrigin);
    pso.setBounds({-5.12, -5.12}, {5.12, 5.12});
    pso.setMode(opt::OptimizationMode::MINIMIZE);

    // 3. Prepare Plotting
    std::string baseName = "pso_vis";
    std::string gridFile = "pso_grid.dat";

    std::cout << "Generating background grid (heatmap)..." << std::endl;
    // Save the static background (function landscape)
    saveFunctionGrid(gridFile, rastrigin, -5.12, 5.12, -5.12, 5.12, 100);

    // 4. Set Callback to save each frame
    pso.setCallback([&](const opt::Solution&, size_t iter) {
        // Use the public getter to access particle positions
        saveSwarmFrame(baseName, iter, pso.getParticles());
        std::cout << "Saved frame " << iter << "/" << config.max_iterations << "\r" << std::flush;
    });

    // 5. Run Optimization
    std::cout << "Running optimization..." << std::endl;
    pso.optimize();
    std::cout << "\nOptimization finished." << std::endl;

    // 6. Launch Animation
    std::cout << "Launching Gnuplot animation..." << std::endl;
    createPSOAnimationScript("run_pso.gp", gridFile, baseName, config.max_iterations, "PSO Rastrigin 2D");
}

void runVisualPSO3DBenchmark() {
    std::cout << "===========================================" << std::endl;
    std::cout << "   Visual PSO Benchmark (3D Animation)" << std::endl;
    std::cout << "===========================================" << std::endl;

    // 1. Configuration
    opt::PSOConfig config;
    config.population_size = 100;
    config.max_iterations = 150;
    opt::PSO pso(config);

    // 2. Objective: 3D Rastrigin Function
    auto rastrigin3D = [](const opt::Coordinates& x) {
        double sum = 0.0;
        double A = 10.0;
        double pi = 3.14159265358979323846;
        for (double val : x) {
            sum += val * val - A * std::cos(2 * pi * val);
        }
        return 3.0 * A + sum;
    };

    pso.setObjectiveFunction(rastrigin3D);

    // Bounds [-5.12, 5.12]
    double min_b = -5.12;
    double max_b = 5.12;
    pso.setBounds({min_b, min_b, min_b}, {max_b, max_b, max_b});
    pso.setMode(opt::OptimizationMode::MINIMIZE);

    // 3. Setup Visualization
    std::string baseName = "pso_vis_3d";
    std::string slicesFile = "pso_slices_3d.dat"; // [NEW] File for wall heatmaps

    // [NEW] Generate the 3D Slices (Walls)
    std::cout << "Generating 3D function slices (walls)..." << std::endl;
    saveFunctionSlices3D(slicesFile, rastrigin3D, min_b, max_b, 50); // 50 = resolution

    // 4. Callback
    pso.setCallback([&](const opt::Solution&, size_t iter) {
        saveSwarmFrame(baseName, iter, pso.getParticles());
        if (iter % 10 == 0) std::cout << "Generating Frame " << iter << "/" << config.max_iterations << "\r" << std::flush;
    });

    // 5. Run
    std::cout << "Running 3D optimization..." << std::endl;
    pso.optimize();
    std::cout << "\nOptimization finished." << std::endl;

    // 6. Launch 3D Animation (Pass slicesFile)
    std::cout << "Launching Gnuplot 3D animation..." << std::endl;
    createPSOAnimationScript3D("run_pso_3d.gp", slicesFile, baseName, config.max_iterations, "PSO 3D Rastrigin", min_b, max_b);
}

// --- Main Optimization Benchmark Entry Point ---

void runOptimizationBenchmarksPSO() {
    std::cout << "=========================================" << std::endl;
    std::cout << "   Particle Swarm Optimization (PSO) Benchmark" << std::endl;
    std::cout << "=========================================" << std::endl;

    // 1. Configuration for PSO
    opt::PSOConfig config;
    config.population_size = 50;   // Number of particles
    config.max_iterations = 100;   // Iterations
    config.inertia_weight = 0.7;
    config.cognitive_coeff = 1.5;
    config.social_coeff = 1.5;

    // 2. Instantiate PSO
    opt::PSO pso(config);

    // 3. Define Shared Bounds [-10, 10]
    opt::Coordinates lower_bounds = {-10.0, -10.0};
    opt::Coordinates upper_bounds = {10.0, 10.0};

    try {
        // Run Test 1
        runSphereTest(pso, lower_bounds, upper_bounds);

        // Run Test 2
        runBoundaryTest(pso, lower_bounds, upper_bounds);

        // Run Test 3: 10-Dimensional Rastrigin
        runRastriginTest(pso, 10);

        // Run visual test:
        runVisualPSOBenchmark();

        // Run 3D visual test
        runVisualPSO3DBenchmark();

    } catch (const std::exception& e) {
        std::cerr << "Optimization failed: " << e.what() << std::endl;
    }
}

// --- Genetic Algorithm (GA) Benchmark ---

// Names file frames: baseName_000.dat, baseName_001.dat, ...
static std::string makeFrameName(const std::string& baseName, size_t iter) {
    std::ostringstream oss;
    oss << baseName << "_iter_" << iter << ".dat";
    return oss.str();
}

static void savePopulationFrame2D(const std::string& baseName, size_t iter,
                                  const std::vector<opt::GA::Individual>& pop) {
    std::string filename = makeFrameName(baseName, iter);
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Warning: Could not create file " << filename << std::endl;
        return;
    }
    for (const auto& ind : pop) {
        if (ind.genome.size() >= 2) {
            out << ind.genome[0] << " " << ind.genome[1] << "\n";
        }
    }
    out.close();
}

static void savePopulationFrame3D(const std::string& baseName, size_t iter,
                                  const std::vector<opt::GA::Individual>& pop) {
    std::string filename = makeFrameName(baseName, iter);
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Warning: Could not create file " << filename << std::endl;
        return;
    }
    for (const auto& ind : pop) {
        if (ind.genome.size() >= 3) {
            out << ind.genome[0] << " " << ind.genome[1] << " " << ind.genome[2] << "\n";
        }
    }
    out.close();
}

void runSphereTest(opt::GA& ga, const opt::Coordinates& lower, const opt::Coordinates& upper) {
    std::cout << "Optimization Problem: Minimize Sphere Function in 2D" << std::endl;
    std::cout << "Search Space: [-10, 10] per dimension" << std::endl;
    std::cout << "Running optimizer..." << std::endl;

    opt::ObjectiveFunction sphere_function = [](const opt::Coordinates& coords) {
        opt::Real sum = 0.0;
        for (auto val : coords) sum += val * val;
        return sum;
    };

    ga.setBounds(lower, upper);
    ga.setObjectiveFunction(sphere_function);
    ga.setMode(opt::OptimizationMode::MINIMIZE);

    ga.setCallback([](const opt::Solution& current_best, size_t iteration) {
        if (iteration == 0 || iteration % 10 == 0) {
            std::cout << "[Sphere Test | Step " << std::setw(3) << iteration << "] "
                      << "Best Value: " << std::scientific << std::setprecision(5)
                      << current_best.value << std::defaultfloat << std::endl;
        }
    });

    auto start = std::chrono::high_resolution_clock::now();
    opt::Solution best_sol = ga.optimize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "\nOptimization Completed in " << duration.count() << " ms." << std::endl;
    std::cout << "Best Value Found: " << std::fixed << std::setprecision(10) << best_sol.value << std::endl;
    std::cout << "Best Position: [ ";
    for (auto val : best_sol.params) {
        std::cout << std::fixed << std::setprecision(5) << val << " ";
    }
    std::cout << "]" << std::endl;
}

void runBoundaryTest(opt::GA& ga, const opt::Coordinates& lower, const opt::Coordinates& upper) {
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "TEST 2: Boundary Constraint Test (Linear Plane)" << std::endl;
    std::cout << "Objective: f(x,y) = x + y (Minimization)" << std::endl;
    std::cout << "Search Space: [-10, 10] per dimension" << std::endl;
    std::cout << "Expected Result: -20.0 at [-10.0, -10.0]" << std::endl;
    std::cout << "Running optimizer..." << std::endl;

    opt::ObjectiveFunction plane_function = [](const opt::Coordinates& coords) {
        opt::Real sum = 0.0;
        for (auto val : coords) sum += val;
        return sum;
    };

    ga.setBounds(lower, upper);
    ga.setObjectiveFunction(plane_function);
    ga.setMode(opt::OptimizationMode::MINIMIZE);

    ga.setCallback([](const opt::Solution& current_best, size_t iteration) {
        if (iteration % 20 == 0) {
            std::cout << "[Boundary Test | Step " << std::setw(3) << iteration << "] "
                      << "Val: " << std::fixed << std::setprecision(4) << current_best.value << std::endl;
        }
    });

    auto start = std::chrono::high_resolution_clock::now();
    opt::Solution best_sol = ga.optimize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Optimization Completed in " << duration.count() << " ms." << std::endl;
    std::cout << "Best Value Found: " << std::fixed << std::setprecision(5) << best_sol.value << std::endl;
    std::cout << "Best Position: [ ";
    for (auto val : best_sol.params) {
        std::cout << std::fixed << std::setprecision(5) << val << " ";
    }
    std::cout << "]" << std::endl;

    if (std::abs(best_sol.value - (-20.0)) < 1e-3) {
        std::cout << ">> SUCCESS: Boundary minimum found correctly!" << std::endl;
    } else {
        std::cout << ">> WARNING: Did not reach the exact boundary." << std::endl;
    }
}

void runRastriginTest(opt::GA& /*ga*/, int dim) {
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "TEST 3: High-Dimensional Stress Test (Rastrigin Function)" << std::endl;
    std::cout << "Dimension: " << dim << "D" << std::endl;
    std::cout << "Search Space: [-5.12, 5.12] per dimension" << std::endl;
    std::cout << "Goal: Find global minimum 0.0 (avoiding local traps)" << std::endl;

    opt::ObjectiveFunction rastrigin_func = [dim](const opt::Coordinates& coords) {
        double sum = 0.0;
        double A = 10.0;
        double pi = 3.14159265358979323846;
        for (auto x : coords) sum += (x * x) - (A * std::cos(2.0 * pi * x));
        return A * dim + sum;
    };

    opt::Coordinates lower(dim, -5.12);
    opt::Coordinates upper(dim,  5.12);

    // Hard configuration for GA (similar to PSO's hard config)
    opt::GAConfig hard_config;
    hard_config.population_size = 800;    // More individuals = more exploration
    hard_config.max_generations = 1200;   // More time
    hard_config.tournament_k    = 4;
    hard_config.crossover_rate  = 0.90;
    hard_config.mutation_rate   = 0.15;
    hard_config.mutation_sigma  = 0.08;
    hard_config.elitism_count   = 3;

    opt::GA local_ga(hard_config);
    local_ga.setBounds(lower, upper);
    local_ga.setObjectiveFunction(rastrigin_func);
    local_ga.setMode(opt::OptimizationMode::MINIMIZE);

    local_ga.setCallback([](const opt::Solution& sol, size_t i) {
        if (i % 100 == 0) {
            std::cout << "[Rastrigin " << i << "] Best: "
                      << std::scientific << std::setprecision(4) << sol.value
                      << std::defaultfloat << std::endl;
        }
    });

    std::cout << "Running optimizer (this might take longer)..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    opt::Solution best_sol = local_ga.optimize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Optimization Completed in " << duration.count() << " ms." << std::endl;
    std::cout << "Best Value Found: " << std::fixed << std::setprecision(5) << best_sol.value << std::endl;

    if (best_sol.value < 1e-2) {
        std::cout << ">> SUCCESS: Global minimum found!" << std::endl;
    } else if (best_sol.value < 5.0) {
        std::cout << ">> ACCEPTABLE: Found a good local minimum, but not global." << std::endl;
    } else {
        std::cout << ">> FAIL: Stuck in a high local minimum." << std::endl;
    }
}

void runVisualGABenchmark() {
    std::cout << "===========================================" << std::endl;
    std::cout << "   Visual GA Benchmark (2D Animation)" << std::endl;
    std::cout << "===========================================" << std::endl;

    opt::GAConfig config;
    config.population_size = 80;
    config.max_generations = 120; // frames
    config.tournament_k    = 3;
    config.crossover_rate  = 0.9;
    config.mutation_rate   = 0.10;
    config.mutation_sigma  = 0.05;
    config.elitism_count   = 2;

    opt::GA ga(config);

    auto rastrigin = [](const opt::Coordinates& x) {
        double A = 10.0;
        double sum = 0.0;
        double pi = 3.14159265358979323846;
        for (double val : x) sum += val*val - A*std::cos(2*pi*val);
        return 2*A + sum;
    };

    ga.setObjectiveFunction(rastrigin);
    ga.setBounds({-5.12, -5.12}, {5.12, 5.12});
    ga.setMode(opt::OptimizationMode::MINIMIZE);

    std::string baseName = "ga_vis";
    std::string gridFile = "ga_grid.dat";

    std::cout << "Generating background grid (heatmap)..." << std::endl;
    saveFunctionGrid(gridFile, rastrigin, -5.12, 5.12, -5.12, 5.12, 100);

    ga.setCallback([&](const opt::Solution&, size_t gen) {
        savePopulationFrame2D(baseName, gen, ga.getPopulation());
        std::cout << "Saved frame " << gen << "/" << config.max_generations << "\r" << std::flush;
    });

    std::cout << "Running optimization..." << std::endl;
    ga.optimize();
    std::cout << "\nOptimization finished." << std::endl;

    std::cout << "Launching Gnuplot animation..." << std::endl;
    // riuso lo script PSO: basta cambiare nome file gp e titolo
    createPSOAnimationScript("run_ga.gp", gridFile, baseName, config.max_generations, "GA Rastrigin 2D");
}

void runVisualGA3DBenchmark() {
    std::cout << "===========================================" << std::endl;
    std::cout << "   Visual GA Benchmark (3D Animation)" << std::endl;
    std::cout << "===========================================" << std::endl;

    opt::GAConfig config;
    config.population_size = 120;
    config.max_generations = 160;
    config.tournament_k    = 3;
    config.crossover_rate  = 0.9;
    config.mutation_rate   = 0.10;
    config.mutation_sigma  = 0.05;
    config.elitism_count   = 2;

    opt::GA ga(config);

    auto rastrigin3D = [](const opt::Coordinates& x) {
        double sum = 0.0;
        double A = 10.0;
        double pi = 3.14159265358979323846;
        for (double val : x) sum += val * val - A * std::cos(2 * pi * val);
        return 3.0 * A + sum;
    };

    double min_b = -5.12;
    double max_b =  5.12;

    ga.setObjectiveFunction(rastrigin3D);
    ga.setBounds({min_b, min_b, min_b}, {max_b, max_b, max_b});
    ga.setMode(opt::OptimizationMode::MINIMIZE);

    std::string baseName   = "ga_vis_3d";
    std::string slicesFile = "ga_slices_3d.dat";

    std::cout << "Generating 3D function slices (walls)..." << std::endl;
    saveFunctionSlices3D(slicesFile, rastrigin3D, min_b, max_b, 50);

    ga.setCallback([&](const opt::Solution&, size_t gen) {
        savePopulationFrame3D(baseName, gen, ga.getPopulation());
        if (gen % 10 == 0) {
            std::cout << "Generating Frame " << gen << "/" << config.max_generations << "\r" << std::flush;
        }
    });

    std::cout << "Running 3D optimization..." << std::endl;
    ga.optimize();
    std::cout << "\nOptimization finished." << std::endl;

    std::cout << "Launching Gnuplot 3D animation..." << std::endl;
    createPSOAnimationScript3D("run_ga_3d.gp", slicesFile, baseName, config.max_generations,
                              "GA 3D Rastrigin", min_b, max_b);
}

void runOptimizationBenchmarksGA() {
    std::cout << "=========================================" << std::endl;
    std::cout << "   Genetic Algorithm (GA) Benchmark" << std::endl;
    std::cout << "=========================================" << std::endl;

    opt::GAConfig config;
    config.population_size = 80;
    config.max_generations = 200;
    config.tournament_k    = 3;
    config.crossover_rate  = 0.9;
    config.mutation_rate   = 0.10;
    config.mutation_sigma  = 0.05;
    config.elitism_count   = 2;

    opt::GA ga(config);

    opt::Coordinates lower_bounds = {-10.0, -10.0};
    opt::Coordinates upper_bounds = { 10.0,  10.0};

    try {
        runSphereTest(ga, lower_bounds, upper_bounds);
        runBoundaryTest(ga, lower_bounds, upper_bounds);
        runRastriginTest(ga, 10);

        // Visual specchio
        runVisualGABenchmark();
        runVisualGA3DBenchmark();

    } catch (const std::exception& e) {
        std::cerr << "Optimization failed: " << e.what() << std::endl;
    }
}

