//
// Created by domenico on 11/27/25.
//

#ifndef MONTECARLO_1_PLOTTER_HPP
#define MONTECARLO_1_PLOTTER_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib> // std::system
#include "../domains/integration_domain.hpp"

/**
 * Utility to close all currently open Gnuplot windows.
 */
inline void closeGnuplotWindows() {
    std::system("pkill -f gnuplot > /dev/null 2>&1");
}

inline std::string formatTitle(std::string name) {
    std::string target = "_samples";
    size_t pos = name.find(target);
    if (pos != std::string::npos) name.erase(pos, target.length());
    for (auto &c : name) if (c == '_') c = ' ';
    return name;
}

/**
 * 1. DOMAIN GEOMETRY PLOT (Visualizes Inside vs Outside points)
 * Creates a unique .dat file and a gnuplot script for the current sample size.
 */
template <size_t dim>
inline void createGnuplotScript(const std::string& tempRawDataFile,
                                const IntegrationDomain<dim>& domain,
                                size_t currentSamples) {

    if (dim > 3) return; // Cannot plot 4D+ geometry easily

    // Prepare unique filenames
    std::string baseName = tempRawDataFile;
    size_t lastDot = baseName.find_last_of(".");
    if (lastDot != std::string::npos) baseName = baseName.substr(0, lastDot);
    size_t lastSlash = baseName.find_last_of("/\\");
    if (lastSlash != std::string::npos) baseName = baseName.substr(lastSlash + 1);

    std::string uniqueID = baseName + "_geom_" + std::to_string(currentSamples);
    std::string uniqueDataFile = uniqueID + ".dat";
    std::string scriptName = "plot_" + uniqueID + ".gp";

    // Read temp raw data and write unique file with Inside/Outside status
    std::ifstream inFile(tempRawDataFile);
    if (!inFile.is_open()) return;

    std::vector<geom::Point<dim>> points;
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        geom::Point<dim> p;
        for (size_t i = 0; i < dim; ++i) ss >> p[i];
        points.push_back(p);
    }
    inFile.close();

    std::ofstream outFile(uniqueDataFile);
    for (const auto& p : points) {
        for (size_t i = 0; i < dim; ++i) outFile << p[i] << " ";
        bool isInside = domain.isInside(p);
        outFile << (isInside ? 1 : 0) << "\n";
    }
    outFile.close();

    // Create Gnuplot script
    std::ofstream gp(scriptName);
    if (!gp.is_open()) return;

    gp << "set grid\n";
    gp << "set title 'Domain Geometry (" << dim << "D) - N=" << currentSamples << "'\n";
    gp << "set key outside\n";

    // Set Axis Ranges
    auto bounds = domain.getBounds();
    double margin = 0.1;
    auto setRange = [&](int axisIdx, std::string axisName) {
        if (axisIdx < dim) {
            double min = bounds[axisIdx].first;
            double max = bounds[axisIdx].second;
            double span = max - min;
            gp << "set " << axisName << "range [" << (min - span * margin) << ":" << (max + span * margin) << "]\n";
        }
    };
    setRange(0, "x");
    if (dim >= 2) setRange(1, "y");
    if (dim >= 3) setRange(2, "z");

    int statusCol = dim + 1; // The added status column index (1-based)

    if (dim == 1) {
        gp << "set xlabel 'X'\n unset ytics\n set yrange [-1:1]\n";
        gp << "plot '" << uniqueDataFile << "' u 1:($" << statusCol << "==0?0:1/0) w p pt 7 ps 0.4 lc rgb 'blue' t 'Out', \\\n";
        gp << "     '" << uniqueDataFile << "' u 1:($" << statusCol << "==1?0:1/0) w p pt 7 ps 0.4 lc rgb 'red' t 'In'\n";
    }
    else if (dim == 2) {
        gp << "set size square\n set xlabel 'X'\n set ylabel 'Y'\n";
        gp << "plot '" << uniqueDataFile << "' u 1:($" << statusCol << "==0?$2:1/0) w p pt 7 ps 0.4 lc rgb 'blue' t 'Out', \\\n";
        gp << "     '" << uniqueDataFile << "' u 1:($" << statusCol << "==1?$2:1/0) w p pt 7 ps 0.4 lc rgb 'red' t 'In'\n";
    }
    else if (dim == 3) {
        gp << "set view equal xyz\n set xlabel 'X'\n set ylabel 'Y'\n set zlabel 'Z'\n set view 60, 30\n";
        gp << "splot '" << uniqueDataFile << "' u 1:2:($" << statusCol << "==0?$3:1/0) w p pt 7 ps 0.4 lc rgb 'blue' t 'Out', \\\n";
        gp << "      '" << uniqueDataFile << "' u 1:2:($" << statusCol << "==1?$3:1/0) w p pt 7 ps 0.4 lc rgb 'red' t 'In'\n";
    }

    gp << "pause mouse close\n"; // Keep window open
    gp.close();

    // Execute Gnuplot in background
    std::string command = "gnuplot " + scriptName + " > /dev/null 2>&1 &";
    std::system(command.c_str());
}

/**
 * 2. PLOT FUNCTION VALUE (Visualizes f(x))
 * - Uses SOLID GREEN color for 1D and 2D domains.
 * - Uses PALETTE colors (heatmap) only for 3D domains (where color represents the 4th dimension).
 */
template <size_t dim, typename Func>
inline void createFunctionGnuplotScript(const std::string& tempRawDataFile,
                                        const IntegrationDomain<dim>& domain,
                                        const Func& func,
                                        size_t currentSamples) {
    if (dim > 3) return;

    std::string baseName = tempRawDataFile;
    size_t lastDot = baseName.find_last_of(".");
    if (lastDot != std::string::npos) baseName = baseName.substr(0, lastDot);
    size_t lastSlash = baseName.find_last_of("/\\");
    if (lastSlash != std::string::npos) baseName = baseName.substr(lastSlash + 1);

    // Unique filename for function values
    std::string uniqueID = baseName + "_func_" + std::to_string(currentSamples);
    std::string uniqueDataFile = uniqueID + ".dat";
    std::string scriptName = "plot_" + uniqueID + ".gp";

    // Read raw data, evaluate f(x) ONLY for points inside the domain, and write to file
    std::ifstream inFile(tempRawDataFile);
    if (!inFile.is_open()) return;

    std::ofstream outFile(uniqueDataFile);

    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        geom::Point<dim> p;
        for (size_t i = 0; i < dim; ++i) ss >> p[i];

        if (domain.isInside(p)) {
            double val = func(p);
            for(size_t i=0; i<dim; ++i) outFile << p[i] << " ";
            outFile << val << "\n";
        }
    }
    inFile.close();
    outFile.close();

    // Generate Gnuplot script
    std::ofstream gp(scriptName);
    if (!gp.is_open()) return;

    gp << "set grid\n";
    gp << "set title 'Function Value (" << dim << "D) - N=" << currentSamples << "'\n";

    // Configure Color Palette logic
    if (dim == 3) {
        // Only 3D Domain (4D total visualization) uses heatmap
        gp << "set palette rgbformulae 33,13,10\n";
        gp << "set colorbox\n";
    } else {
        // 1D and 2D Domains use solid green
        gp << "unset colorbox\n";
    }

    int valCol = dim + 1; // Column containing function value

    if (dim == 1) {
        gp << "set xlabel 'X'\n set ylabel 'f(X)'\n";
        // 1D Domain -> 2D Plot (x, f(x)) -> Solid Green
        gp << "plot '" << uniqueDataFile << "' u 1:" << valCol << " w p pt 7 ps 0.4 lc rgb 'green' t 'f(x)'\n";
    }
    else if (dim == 2) {
        gp << "set size square\n set view 60, 30\n";
        gp << "set xlabel 'X'\n set ylabel 'Y'\n set zlabel 'f(X,Y)'\n";
        // 2D Domain -> 3D Plot (x, y, z=f(x,y)) -> Solid Green
        gp << "splot '" << uniqueDataFile << "' u 1:2:" << valCol << " w p pt 7 ps 0.4 lc rgb 'green' t 'f(x,y)'\n";
    }
    else if (dim == 3) {
        gp << "set view equal xyz\n set view 60, 30\n";
        gp << "set xlabel 'X'\n set ylabel 'Y'\n set zlabel 'Z'\n";
        // 3D Domain -> 4D Visualization (x, y, z, color=f(x,y,z)) -> Palette
        gp << "splot '" << uniqueDataFile << "' u 1:2:3:" << valCol << " w p pt 7 ps 0.4 lc palette t 'f(x,y,z)'\n";
    }

    gp << "pause mouse close\n";
    gp.close();

    std::string command = "gnuplot " + scriptName + " > /dev/null 2>&1 &";
    std::system(command.c_str());
}

#endif //MONTECARLO_1_PLOTTER_HPP