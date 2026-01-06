//
// Created by domenico on 12/3/25.
//

#ifndef MONTECARLO_1_BENCHMARKS_HPP
#define MONTECARLO_1_BENCHMARKS_HPP

#include <montecarlo/domains/hypercylinder.hpp>
#include <montecarlo/domains/hyperrectangle.hpp>
#include <montecarlo/domains/hypersphere.hpp>
#include <montecarlo/geometry.hpp>
#include <montecarlo/integrators/montecarlo_integrator.hpp>
#include <montecarlo/proposals/uniformProposal.hpp>
#include <montecarlo/utils/plotter.hpp>

#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <functional>
#include <thread>
#include <cmath>

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

// Function that call all the others
void runBenchmarks(bool useGnuplot);

// executes all benchmarks on parsed function.txt's function
void runBenchmarks(const std::string& expression,bool useGnuplot);


#endif //MONTECARLO_1_BENCHMARKS_HPP