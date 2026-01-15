//
// Created by domenico on 12/3/25.
//

#ifndef MONTECARLO_1_BENCHMARKS_HPP
#define MONTECARLO_1_BENCHMARKS_HPP

#include <montecarlo/domains/hypercylinder.hpp>
#include <montecarlo/domains/hyperrectangle.hpp>
#include <montecarlo/domains/hypersphere.hpp>
#include <montecarlo/geometry.hpp>
#include <montecarlo/integrators/MCintegrator.hpp>
#include <montecarlo/integrators/MHintegrator.hpp>
#include <montecarlo/integrators/ISintegrator.hpp>
#include <montecarlo/proposals/uniformProposal.hpp>
#include <montecarlo/utils/plotter.hpp>
#include <montecarlo/optimizers/PSO.hpp>
#include <montecarlo/optimizers/GA.hpp>

#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <functional>
#include <thread>
#include <cmath>
#include <cstdint>

// Global seed for all random number generation
extern uint32_t GLOBAL_SEED;

// Global benchmark configuration
extern const std::vector<size_t> n_samples_vector;
extern unsigned int n_threads;

// Struct used to write through the save results function the results in a file of name file_name
struct results {
    size_t n_samples;
    std::string integration_result;
    std::string duration;
};

void saveResults(const std::string &filename, const std::vector<results> &results, const std::string &function_expr);

// Functions that contain both the domains and the functions to integrate over them
void uniDimIntegration();
void circleIntegration();
void sphereIntegration();
void rectangularIntegration();
void cylinderIntegration();
void parallelepipedIntegration();
void fiveDimIntegration();
void fourDimIntegration();
void eightDimIntegration();
void twelveDimIntegration();

// --- Monte Carlo Integration Benchmarks ---
// (Implemented in benchmarks/integration_benchmarks.cpp)
void runBenchmarks(bool useGnuplot);
void runBenchmarks(const std::string& expression, bool useGnuplot);
void runBenchmarksMH(bool useGnuplot);

// --- PSO Optimization Benchmarks ---
// (Implemented in benchmarks/pso_benchmarks.cpp)
void runOptimizationBenchmarksPSO();

// --- GA Optimization Benchmarks ---
// (Implemented in benchmarks/ga_benchmarks.cpp)
void runOptimizationBenchmarksGA();


#endif //MONTECARLO_1_BENCHMARKS_HPP