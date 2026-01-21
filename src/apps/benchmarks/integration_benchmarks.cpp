//
// Integration benchmarks for Monte Carlo methods
//

#include "apps/benchmarks.hpp"
#include <montecarlo/integrators/ISintegrator.hpp>
#include <montecarlo/utils/muParserXInterface.hpp>
#include <cmath>
#include <cstdint>

using namespace MuParserInterface;

// --- Helper for formatted console output ---

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
    // --- Helper for console output ---
    auto printLine = [&](std::size_t n,
                     const std::string& label,
                     double result,
                     long time_ms)
    {
        std::cout << std::setw(8)  << n << " | "
                << std::setw(30) << label << " | "
                << std::setw(20) << std::fixed << std::setprecision(6) << result << " | "
                << std::setw(5)  << time_ms << " ms\n";
    };

    std::cout << "------------------------------------------------" << std::endl;
    std::cout << title << ":" << std::endl;
    std::cout << "Integrating Function: " << functionExpr << std::endl << std::endl;

    std::vector<results> testResults;

	ISMontecarloIntegrator<dim> isIntegrator(domain);
    UniformProposal<dim> uprop(domain);
    std::vector<double> init_mean(dim, 0.0);
    std::vector<double> init_sigma(dim, 2.5);
    auto bounds = domain.getBounds();
    for (size_t i = 0; i < dim; ++i) {
        init_mean[i]  = 0.5 * (bounds[i].first + bounds[i].second);
        init_sigma[i] = (bounds[i].second - bounds[i].first) / 3.0; // oppure /2.0
    }
    GaussianProposal<dim> gprop(domain, init_mean, init_sigma);
    MixtureProposal<dim> mix({&uprop, &gprop}, {0.5, 0.5});

    // Header table for console output
    std::cout << std::string(107, '-') << '\n';
    std::cout << std::setw(8)  << "Samples" << " | "
            << std::setw(30) << "Method"  << " | "
            << std::setw(20) << "Result"  << " | "
            << std::setw(6)  << "Time"    << '\n';
    std::cout << std::string(107, '-') << '\n';

    for (size_t n_i : n_samples_vector) {
        // 1. Normal MC (Data points are written to rawDataFile by the integrator)
        auto startTimer1 = std::chrono::high_resolution_clock::now();
        double result1 = integrator.OLDintegrate(f, n_i);
        auto endTimer1 = std::chrono::high_resolution_clock::now();
        auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer1 - startTimer1);

        // 2. Uniform IS (Data points are written to rawDataFile by the integrator)
        auto startTimer2 = std::chrono::high_resolution_clock::now();
        double result2 = isIntegrator.integrate(f, n_i, uprop, 12345);
        auto endTimer2 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer2 - startTimer2);

        // 3. Gaussian IS (Data points are written to rawDataFile by the integrator)
        auto startTimer3 = std::chrono::high_resolution_clock::now();
        double result3 = isIntegrator.integrate(f, n_i, gprop, 12345);
        auto endTimer3 = std::chrono::high_resolution_clock::now();
        auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer3 - startTimer3);

        // 4. Mixture IS (Data points are written to rawDataFile by the integrator)
        auto startTimer4 = std::chrono::high_resolution_clock::now();
        double result4 = isIntegrator.integrate(f, n_i, mix, 12345);
        auto endTimer4 = std::chrono::high_resolution_clock::now();
        auto duration4 = std::chrono::duration_cast<std::chrono::milliseconds>(endTimer4 - startTimer4);

        // Store results
        results newLine;
        newLine.n_samples = n_i;
        newLine.integration_result = std::to_string(result1);
        newLine.duration = std::to_string(duration1.count());
        testResults.push_back(newLine); 

        // Console Output
        printLine(n_i, "Normal Sampling",              result1, duration1.count());
        printLine(n_i, "Uniform Importance Sampling",  result2, duration2.count());
        printLine(n_i, "Gaussian Importance Sampling", result3, duration3.count());
        printLine(n_i, "Mixture Importance Sampling",  result4, duration4.count());
        std::cout << std::string(107, '-') << '\n';

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
	MHMontecarloIntegrator<2> mhintegrator(domain);
 	UniformProposal<2> dummy_proposal(domain);

	mhintegrator.setConfig(
        burn_in,
        thinning,
        n_samples_volume,
        deviation,
        indicator,
        x0
    );

    std::cout << "Running Benchmarks" << std::endl;
    std::cout << "Metropolis Hastings result: " << mhintegrator.integrate(f, static_cast<int>(n_samples), dummy_proposal, std::random_device{}()) << std::endl;
    std::cout << "Montecarlo result: " << integrator.OLDintegrate(f, n_samples) << std::endl;
    std::cout << "Exact result: " << 5000 * M_PI << std::endl;
}

void runBenchmarks(const std::string& expression, bool useGnuplot) {
    n_threads = std::thread::hardware_concurrency();
    if (n_threads == 0) n_threads = 16;

    circleIntegrationParser(expression, useGnuplot);
    sphereIntegrationParser(expression, useGnuplot);
    rectangularIntegrationParser(expression, useGnuplot);
    cylinderIntegrationParser(expression, useGnuplot);
}
