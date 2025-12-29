//
// Created by domenico on 12/3/25.
//
#include "benchmarks.hpp"
#include "utils/muParserXInterface.hpp"

using namespace MuParserInterface;

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

        // Store results
        results newLine;
        newLine.n_samples = n_i;
        newLine.integration_result = std::to_string(result);
        newLine.duration = std::to_string(duration.count());
        testResults.push_back(newLine);

        // Print to console
        std::cout << std::setw(12) << n_i << " Samples"
                  << std::setw(18) << std::fixed << std::setprecision(6) << result
                  << std::setw(18) << duration.count() << "ms" << std::endl;

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