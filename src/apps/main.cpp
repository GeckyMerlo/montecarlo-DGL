#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <array>
#include <vector>
#include <stdexcept>

#include "apps/benchmarks.hpp"
#include <montecarlo/utils/plotter.hpp>
#include <montecarlo/integrators/montecarlo_integrator.hpp>
#include <montecarlo/domains/polytope.hpp>
#include <montecarlo/geometry.hpp>

// Path to the file containing the function.
#define FUNCTION_FILE "../function.txt"

constexpr int dim = 2;


// --- FUNCTION PROTOTYPES ---
std::string readFunctionFromFile(const std::string& filename);

template <int dim>
std::vector<Point<dim>> read_points_from_file(const std::string& filename);

template <int dim>
void read_normals_and_offsets_from_qhull_n(
    const std::string& filename,
    std::vector<std::array<double, dim>>& normals,
    std::vector<double>& offsets);


// --- MAIN ---

int main() {
    // Closes already open gnuplot windows
    closeGnuplotWindows();

    std::cout << "===========================================" << std::endl;
    std::cout << "   Monte Carlo Integration Benchmarks" << std::endl;
    std::cout << "===========================================" << std::endl;

    // Mode choice
    std::cout << "Select mode:" << std::endl;
    std::cout << "1. Use function from file (function.txt) - Uses Parser (Slower)" << std::endl;
    std::cout << "2. Use hardcoded function - Uses C++ Lambda (Faster)" << std::endl;
    std::cout << "3. Do you want to use a polytope con covex hull (IT: Inviluppo Convesso)" << std::endl;
    std::cout << "Choice: ";

    int choice;
    if (!(std::cin >> choice)) {
        std::cerr << "Invalid input." << std::endl;
        return 1;
    }
    // Clears buffer until next line
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Choice of visualizing with gnuPlot
    std::cout << "Enable Gnuplot visualization for results? (y/n): ";
    char gpChoice;
    std::cin >> gpChoice;
    // Clears buffer
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    bool useGnuplot = (gpChoice == 'y' || gpChoice == 'Y');

    if (choice == 1) {
        try {
            std::string expression = readFunctionFromFile(FUNCTION_FILE);
            std::cout << "\nLoaded expression: " << expression << std::endl;
            std::cout << "Starting PARSER benchmarks..." << std::endl;
            // Passa il flag gnuplot
            runBenchmarks(expression, useGnuplot);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    } else if (choice == 2) {
        std::cout << "\nStarting HARDCODED benchmarks..." << std::endl;
        // Sends file to gnuPlot
        runBenchmarks(useGnuplot);
    }else if (choice == 3) {
        std::cout << "\nReading Points, Normals and Offsets..." << std::endl;

        std::vector<geom::Point<dim>> points = read_points_from_file<dim>("../points.txt");


        std::vector<std::array<double, dim>> normals;
        std::vector<double> offsets;

        read_normals_and_offsets_from_qhull_n<dim>("../hull.txt", normals, offsets);

        PolyTope<dim> polytope(points, normals, offsets);
        MontecarloIntegrator<dim> integrator(polytope);
        /*
        auto f_const = [](const Point<3>& p) {
            return 1.0;
        };

        auto f_linear = [](const Point<3>& p) {
            return p[0] + p[1] + p[2];   // x + y + z
        };

        auto f_quad = [](const Point<3>& p) {
            return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
        };
        std::cout << integrator.integrate(f_const, 1000000) << std::endl;
        std::cout << "Expected: 1.0" << std::endl;
        std::cout << integrator.integrate(f_quad, 1000000) << std::endl;
        std::cout << "Expected: 1.0" << std::endl;
        std::cout << integrator.integrate(f_linear, 1000000) << std::endl;
        std::cout << "Expected: 1.5" << std::endl;
        */


        // f(x,y) = 1
        auto f_const = [](const Point<2>& p) {
            return 1.0;
        };

        // f(x,y) = x
        auto f_x = [](const Point<2>& p) {
            return p[0];
        };

        // f(x,y) = y
        auto f_y = [](const Point<2>& p) {
            return p[1];
        };

        double I_const = integrator.integrate(f_const, 1000000);
        double I_x     = integrator.integrate(f_x,     1000000);
        double I_y     = integrator.integrate(f_y,     1000000);

        std::cout << "Integral f=1   ≈ " << I_const << "  (exact: " << 3*std::sqrt(3)/2 << ")\n";
        std::cout << "Integral f=x   ≈ " << I_x     << "  (exact: 0)\n";
        std::cout << "Integral f=y   ≈ " << I_y     << "  (exact: 0)\n";

    }else {
        std::cerr << "Invalid choice selected." << std::endl;
        return 1;
    }

    std::cout << "\nAll benchmarks completed." << std::endl;
    if (useGnuplot) {
        std::cout << "Check opened Gnuplot windows. Press Enter in console or close windows to finish fully if needed." << std::endl;
    }

    return 0;
}

// --- FUNCTION IMPLEMENTATIONS ---

std::string readFunctionFromFile(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open function file at: " + filename +
                                 "\nMake sure the file exists in the repository root.");
    }

    std::string expression;
    if (!std::getline(file, expression)) {
        throw std::runtime_error("File is empty: " + filename);
    }

    if (!expression.empty() && expression.back() == '\r') {
        expression.pop_back();
    }

    file.close();

    if (expression.empty()) {
        throw std::runtime_error("Expression in file is empty: " + filename);
    }

    return expression;
}

template <int dim>
std::vector<Point<dim>> read_points_from_file(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::size_t num_points = 0;
    std::size_t file_dim   = 0;

    // Prima riga: <num_points> <dim>
    in >> num_points >> file_dim;

    if (!in.good()) {
        throw std::runtime_error("Error reading header from file: " + filename);
    }

    if (file_dim != static_cast<std::size_t>(dim)) {
        throw std::runtime_error(
            "Dimension mismatch: file has dim = " + std::to_string(file_dim) +
            " but template expects dim = " + std::to_string(dim));
    }

    std::vector<Point<dim>> points;
    points.reserve(num_points);

    for (std::size_t i = 0; i < num_points; ++i) {
        Point<dim> p;
        for (int k = 0; k < dim; ++k) {
            if (!(in >> p[k])) {
                throw std::runtime_error(
                    "Error reading coordinate " + std::to_string(k) +
                    " of point " + std::to_string(i) +
                    " from file: " + filename);
            }
        }
        points.push_back(p);
    }

    return points;
}

template <int dim>
void read_normals_and_offsets_from_qhull_n(
    const std::string& filename,
    std::vector<std::array<double, dim>>& normals,
    std::vector<double>& offsets)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open normals file: " + filename);
    }

    std::size_t file_dim = 0;
    std::size_t num_facets = 0;

    // Legge dimensione e numero di facce (possono essere sulla stessa riga o su due righe)
    in >> file_dim >> num_facets;
    if (!in.good()) {
        throw std::runtime_error("Error reading header (dim, num_facets) from: " + filename);
    }

    if (file_dim != dim + 1) {
        throw std::runtime_error(
            "Dimension mismatch in normals file: file has dim = " +
            std::to_string(file_dim) + " but template expects dim = " +
            std::to_string(dim));
    }

    normals.clear();
    offsets.clear();
    normals.reserve(num_facets);
    offsets.reserve(num_facets);

    for (std::size_t f = 0; f < num_facets; ++f) {
        std::array<double, dim> n{};
        double d = 0.0;

        // Legge i componenti della normale
        for (std::size_t k = 0; k < dim; ++k) {
            if (!(in >> n[k])) {
                throw std::runtime_error(
                    "Error reading normal component " + std::to_string(k) +
                    " for facet " + std::to_string(f) +
                    " from: " + filename);
            }
        }

        // Legge l'offset d (nel piano n·x + d = 0)
        if (!(in >> d)) {
            throw std::runtime_error(
                "Error reading offset d for facet " + std::to_string(f) +
                " from: " + filename);
        }

        double b = -d;

        normals.push_back(n);
        offsets.push_back(b);
    }
}