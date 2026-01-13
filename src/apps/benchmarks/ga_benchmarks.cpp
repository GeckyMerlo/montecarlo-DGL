//
// GA (Genetic Algorithm) benchmarks
//

#include "apps/benchmarks.hpp"
#include <cmath>
#include <cstdint>
#include <sstream>
#include <filesystem>

namespace opt = optimizers;

// --- GA Helper Functions ---

// Names file frames: baseName_iter_0.dat, baseName_iter_1.dat, ...
static std::string makeFrameName(const std::string& baseName, size_t iter) {
    std::ostringstream oss;
    oss << baseName << "_iter_" << iter << ".dat";
    return oss.str();
}

static void savePopulationFrame2D(const std::string& baseName, size_t iter,
                                  const std::vector<opt::GA::Individual>& pop) {
    // Create ga_frames subdirectory
    std::string dir = "./ga_frames";
    try {
        std::filesystem::create_directories(dir);
    } catch (const std::exception&) {}
    
    std::string filename = dir + "/" + makeFrameName(baseName, iter);
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
    // Create ga_frames subdirectory
    std::string dir = "./ga_frames";
    try {
        std::filesystem::create_directories(dir);
    } catch (const std::exception&) {}
    
    std::string filename = dir + "/" + makeFrameName(baseName, iter);
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

// --- GA Test Functions ---

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
