/**
 * @file drone_optimization_ellipsoid_strength.cpp
 * @brief Drone arm CM targeting + vertical bending-stress proxy with an ellipsoidal hole.
 *
 * @details
 * This example demonstrates the full power of the Monte Carlo library stack:
 * - Complex 3D domain with CSG union (arm + motor + optional cabin polytope)
 * - Stochastic objective via Monte Carlo integration (mass, CM, second-moment proxies)
 * - Deterministic parameter-based seeding via RngManager for reproducible MC noise
 * - Non-smooth constraints handled by penalties (hole must remain inside the body)
 * - Comparison of PSO vs GA optimization strategies
 *
 * Optimization variables (ellipsoidal hole):
 *   params = [hx, hy, hz, a, b, c, yaw, pitch, roll]
 * where (hx,hy,hz) is the hole center, (a,b,c) are semi-axes, and yaw/pitch/roll define rotation.
 *
 * Objectives:
 *   1) Minimize CM error from target (1,0,0)
 *   2) Minimize a vertical bending-stress proxy for a cantilever-like arm
 *
 * Scalar objective (single-objective optimization):
 *   f = w_cm * cm_error + w_sigma * sigma_proxy + penalties
 *
 * @see PSO, GA, MontecarloIntegrator, MCMeanEstimator, VolumeEstimatorMC, RngManager
 */

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <iomanip>
#include <limits>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <fstream>
#include <filesystem>

#include <omp.h>

// Library headers
#include "../montecarlo/geometry.hpp"
#include "../montecarlo/RngManager.hpp"
#include "../montecarlo/domains/integration_domain.hpp"
#include "../montecarlo/domains/hyperrectangle.hpp"
#include "../montecarlo/domains/hypercylinder.hpp"
#include "../montecarlo/domains/polytope.hpp"
#include "../montecarlo/integrators/MCintegrator.hpp"
#include "../montecarlo/estimators/VolumeEstimatorMC.hpp"
#include "../montecarlo/estimators/MCMeanEstimator.hpp"
#include "../montecarlo/proposals/uniformProposal.hpp"
#include "../montecarlo/optimizers/PSO.hpp"
#include "../montecarlo/optimizers/GA.hpp"
#include "../montecarlo/utils/plotter.hpp"

using namespace geom;
using namespace optimizers;

/// Global seed for deterministic behavior (also used by optimizers)
uint32_t GLOBAL_SEED = 12345;

/// Problem geometry dimension
constexpr size_t DIM = 3;

// =====================================================================================
// 0) Helpers: deterministic hashing, 3x3 matrix algebra for ellipsoid rotation
// =====================================================================================

/**
 * @brief Deterministic hash of parameter vector for reproducible MC noise.
 * @details Ensures that identical parameters produce identical random sequences,
 *          making the stochastic objective function stationary for PSO/GA.
 */
static uint32_t hashParams(const Coordinates& params) {
    uint32_t seed = 0;
    std::hash<double> hasher;
    for (const auto& p : params) {
        seed ^= static_cast<uint32_t>(hasher(p)) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
    }
    return seed;
}

/// 3x3 matrix for rotation operations
struct Mat3 {
    double m[3][3]{};
};

/// Matrix multiplication: C = A * B
static Mat3 matMul(const Mat3& A, const Mat3& B) {
    Mat3 C;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            C.m[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) 
                C.m[i][j] += A.m[i][k] * B.m[k][j];
        }
    return C;
}

/// Matrix transpose
static Mat3 matTranspose(const Mat3& A) {
    Mat3 T;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            T.m[i][j] = A.m[j][i];
    return T;
}

/// Matrix-vector multiplication: result = A * v
static std::array<double, 3> matVec(const Mat3& A, const std::array<double, 3>& v) {
    return {
        A.m[0][0]*v[0] + A.m[0][1]*v[1] + A.m[0][2]*v[2],
        A.m[1][0]*v[0] + A.m[1][1]*v[1] + A.m[1][2]*v[2],
        A.m[2][0]*v[0] + A.m[2][1]*v[1] + A.m[2][2]*v[2]
    };
}

/**
 * @brief Build rotation matrix from Euler angles (Z-Y-X convention).
 * @param yaw   Rotation around Z axis
 * @param pitch Rotation around Y axis
 * @param roll  Rotation around X axis
 * @return Combined rotation matrix R = Rz * Ry * Rx
 */
static Mat3 rotationZYX(double yaw, double pitch, double roll) {
    const double cy = std::cos(yaw),   sy = std::sin(yaw);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double cr = std::cos(roll),  sr = std::sin(roll);

    Mat3 Rz{{ {cy, -sy, 0.0}, {sy, cy, 0.0}, {0.0, 0.0, 1.0} }};
    Mat3 Ry{{ {cp, 0.0, sp},  {0.0, 1.0, 0.0}, {-sp, 0.0, cp} }};
    Mat3 Rx{{ {1.0, 0.0, 0.0}, {0.0, cr, -sr}, {0.0, sr, cr} }};

    return matMul(matMul(Rz, Ry), Rx);
}

/**
 * @brief Precomputed ellipsoid data for fast inside-tests.
 * @details Stores inverse squared semi-axes and transposed rotation matrix
 *          to avoid recomputation during millions of point tests.
 */
struct EllipsoidCache {
    double hx, hy, hz;               // Ellipsoid center
    double inv_a2, inv_b2, inv_c2;   // 1/a², 1/b², 1/c²
    Mat3 Rt;                         // R^T (transpose of rotation)
};

/// Build ellipsoid cache from parameter vector
static EllipsoidCache buildEllipsoidCache(const Coordinates& params) {
    EllipsoidCache E;
    E.hx = params[0]; 
    E.hy = params[1]; 
    E.hz = params[2];

    const double a = params[3], b = params[4], c = params[5];
    E.inv_a2 = 1.0 / (a * a);
    E.inv_b2 = 1.0 / (b * b);
    E.inv_c2 = 1.0 / (c * c);

    E.Rt = matTranspose(rotationZYX(params[6], params[7], params[8]));
    return E;
}

/// Fast ellipsoid inside-test using precomputed cache
static inline bool isInsideEllipsoid(const Point<DIM>& p, const EllipsoidCache& E) {
    const double dx = p[0] - E.hx;
    const double dy = p[1] - E.hy;
    const double dz = p[2] - E.hz;

    // Transform to ellipsoid local coords: q = R^T * (p - center)
    const double qx = E.Rt.m[0][0]*dx + E.Rt.m[0][1]*dy + E.Rt.m[0][2]*dz;
    const double qy = E.Rt.m[1][0]*dx + E.Rt.m[1][1]*dy + E.Rt.m[1][2]*dz;
    const double qz = E.Rt.m[2][0]*dx + E.Rt.m[2][1]*dy + E.Rt.m[2][2]*dz;

    // Ellipsoid equation: (qx/a)² + (qy/b)² + (qz/c)² ≤ 1
    const double u = qx*qx*E.inv_a2 + qy*qy*E.inv_b2 + qz*qz*E.inv_c2;
    return u <= 1.0;
}

// =====================================================================================
// 1) PolyTope file readers (cabin geometry)
// =====================================================================================

template <int dim>
static std::vector<Point<dim>> readPointsFromFile(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) 
        throw std::runtime_error("Cannot open file: " + filename);

    std::size_t num_points = 0, file_dim = 0;
    in >> num_points >> file_dim;
    if (!in.good()) 
        throw std::runtime_error("Error reading header from: " + filename);

    if (file_dim != static_cast<std::size_t>(dim))
        throw std::runtime_error("Dimension mismatch in " + filename);

    std::vector<Point<dim>> points;
    points.reserve(num_points);

    for (std::size_t i = 0; i < num_points; ++i) {
        Point<dim> p;
        for (int k = 0; k < dim; ++k) {
            if (!(in >> p[k]))
                throw std::runtime_error("Error reading point " + std::to_string(i) + " from " + filename);
        }
        points.push_back(p);
    }
    return points;
}

template <int dim>
static void readNormalsAndOffsets(const std::string& filename,
                                   std::vector<std::array<double, dim>>& normals,
                                   std::vector<double>& offsets) {
    std::ifstream in(filename);
    if (!in.is_open()) 
        throw std::runtime_error("Cannot open normals file: " + filename);

    std::size_t file_dim = 0, num_facets = 0;
    in >> file_dim >> num_facets;
    if (!in.good()) 
        throw std::runtime_error("Error reading header from: " + filename);

    if (file_dim != static_cast<std::size_t>(dim + 1))
        throw std::runtime_error("Dimension mismatch in normals file: " + filename);

    normals.clear();
    offsets.clear();
    normals.reserve(num_facets);
    offsets.reserve(num_facets);

    for (std::size_t f = 0; f < num_facets; ++f) {
        std::array<double, dim> n{};
        double d = 0.0;

        for (std::size_t k = 0; k < dim; ++k) {
            if (!(in >> n[k]))
                throw std::runtime_error("Error reading normal from: " + filename);
        }
        if (!(in >> d)) 
            throw std::runtime_error("Error reading offset from: " + filename);

        normals.push_back(n);
        offsets.push_back(-d);  // Convert from n·x + d ≤ 0 to n·x ≤ -d
    }
}

// =====================================================================================
// 2) DroneArmDomain: CSG union of Arm + Motor + optional Cabin (PolyTope)
// =====================================================================================

/**
 * @brief Complex 3D domain representing a drone arm with motor housing.
 * @details Implements IntegrationDomain interface for use with Monte Carlo estimators.
 *          Uses CSG union: isInside = arm ∪ motor_housing ∪ cabin
 */
class DroneArmDomain : public IntegrationDomain<DIM> {
public:
    // Geometric components
    std::unique_ptr<HyperRectangle<DIM>> arm;
    std::unique_ptr<HyperCylinder<DIM>> motor_housing;
    std::unique_ptr<PolyTope<DIM>> cabin;

    // Component offsets from origin
    std::array<double, DIM> motor_offset{};
    std::array<double, DIM> cabin_offset{};

    // Arm dimensions (public for stress surrogate calculations)
    double arm_length = 10.0;
    double arm_width  = 2.0;
    double arm_height = 1.0;

    DroneArmDomain() {
        // 1. Arm: rectangular beam centered at origin
        std::array<double, DIM> arm_dims = {arm_length, arm_width, arm_height};
        arm = std::make_unique<HyperRectangle<DIM>>(arm_dims);

        // 2. Motor housing: cylinder at arm tip
        motor_housing = std::make_unique<HyperCylinder<DIM>>(1.5, 1.2);
        motor_offset = {+arm_length / 2.0, 0.0, -0.6};

        // 3. Optional cabin from external files
        try {
            auto points = readPointsFromFile<DIM>("../drone_assets/cabin_points.txt");
            std::vector<std::array<double, DIM>> normals;
            std::vector<double> offsets;
            readNormalsAndOffsets<DIM>("../drone_assets/cabin_hull.txt", normals, offsets);

            cabin = std::make_unique<PolyTope<DIM>>(points, normals, offsets);
            cabin_offset = {-2.0, 0.0, 0.5};
        } catch (...) {
            cabin = nullptr;
        }
    }

    Bounds<DIM> getBounds() const override {
        Bounds<DIM> b;
        b[0] = {-6.0, 8.0};
        b[1] = {-3.0, 3.0};
        b[2] = {-2.0, 2.0};
        return b;
    }

    double getBoxVolume() const override {
        return 14.0 * 6.0 * 4.0;
    }

    bool isInside(const Point<DIM>& p) const override {
        // CSG Union: check arm first (most common)
        if (arm->isInside(p)) return true;

        // Check motor (with offset transformation)
        Point<DIM> p_motor;
        for (size_t i = 0; i < DIM; ++i) 
            p_motor[i] = p[i] - motor_offset[i];
        if (motor_housing->isInside(p_motor)) return true;

        // Check cabin (optional polytope)
        if (cabin) {
            Point<DIM> p_cabin;
            for (size_t i = 0; i < DIM; ++i) 
                p_cabin[i] = p[i] - cabin_offset[i];
            if (cabin->isInside(p_cabin)) return true;
        }

        return false;
    }

    /**
     * @brief Export domain geometry with ellipsoidal hole for visualization.
     */
    void exportGeometryWithHole(const std::string& out_dir,
                                 const Coordinates& params) const {
        try { 
            std::filesystem::create_directories(out_dir); 
        } catch (...) {}

        const std::string filename = out_dir + "/drone_domain_ellipsoid.txt";
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "Cannot open export file: " << filename << "\n";
            return;
        }

        const EllipsoidCache E = buildEllipsoidCache(params);

        out << "# Drone domain with ellipsoidal hole\n";
        out << "# Ellipsoid: center=(" << params[0] << "," << params[1] << "," << params[2] << ")\n";
        out << "#            axes=(" << params[3] << "," << params[4] << "," << params[5] << ")\n";
        out << "#            angles(yaw,pitch,roll)=(" << params[6] << "," << params[7] << "," << params[8] << ")\n";
        out << "# Format: type x y z\n";

        const double step = 0.15;

        // Sample arm region
        for (double x = -arm_length/2.0; x <= arm_length/2.0; x += step) {
            for (double y = -arm_width/2.0; y <= arm_width/2.0; y += step) {
                for (double z = -arm_height/2.0; z <= arm_height/2.0; z += step) {
                    Point<DIM> p; 
                    p[0] = x; p[1] = y; p[2] = z;
                    if (arm->isInside(p) && !isInsideEllipsoid(p, E))
                        out << "arm " << x << " " << y << " " << z << "\n";
                }
            }
        }

        // Sample motor region
        for (double x = motor_offset[0] - 2.0; x <= motor_offset[0] + 2.0; x += step) {
            for (double y = -2.0; y <= 2.0; y += step) {
                for (double z = -1.5; z <= 1.5; z += step) {
                    Point<DIM> p; 
                    p[0] = x; p[1] = y; p[2] = z;
                    Point<DIM> pl;
                    for (size_t i = 0; i < DIM; ++i) 
                        pl[i] = p[i] - motor_offset[i];
                    if (motor_housing->isInside(pl) && !isInsideEllipsoid(p, E))
                        out << "motor " << x << " " << y << " " << z << "\n";
                }
            }
        }

        std::cout << "Exported geometry to: " << filename << "\n";
    }
};

// =====================================================================================
// 3) Hole containment penalty (surface sampling)
// =====================================================================================

/**
 * @brief Check if ellipsoid hole is fully contained within the domain.
 * @details Samples points on the ellipsoid surface and returns the fraction
 *          that fall outside the domain (violation ratio).
 * @param domain The drone arm domain
 * @param params Ellipsoid parameters [hx,hy,hz,a,b,c,yaw,pitch,roll]
 * @param rng Random number generator for surface sampling
 * @param n_samples Number of surface samples to check
 * @return Fraction of surface points outside domain [0,1]
 */
static double computeContainmentPenalty(const DroneArmDomain& domain,
                                         const Coordinates& params,
                                         std::mt19937& rng,
                                         int n_samples) {
    // Quick check: center must be inside
    Point<DIM> center;
    center[0] = params[0];
    center[1] = params[1];
    center[2] = params[2];
    if (!domain.isInside(center)) 
        return 1.0;

    const double a = params[3], b = params[4], c = params[5];
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) 
        return 1.0;

    const Mat3 R = rotationZYX(params[6], params[7], params[8]);

    std::uniform_real_distribution<double> u11(-1.0, 1.0);
    int violations = 0;

    for (int i = 0; i < n_samples; ++i) {
        // Marsaglia method for uniform sphere sampling (no trig)
        double u, v, s;
        do {
            u = u11(rng);
            v = u11(rng);
            s = u*u + v*v;
        } while (s >= 1.0 || s < 1e-12);

        const double factor = std::sqrt(1.0 - s);
        const double sx = 2.0 * u * factor;
        const double sy = 2.0 * v * factor;
        const double sz = 1.0 - 2.0 * s;

        // Point on ellipsoid surface in local coords, then rotate to global
        std::array<double, 3> local = {a * sx, b * sy, c * sz};
        auto global = matVec(R, local);

        Point<DIM> p;
        p[0] = params[0] + global[0];
        p[1] = params[1] + global[1];
        p[2] = params[2] + global[2];

        if (!domain.isInside(p)) 
            violations++;
    }

    return static_cast<double>(violations) / static_cast<double>(n_samples);
}

// =====================================================================================
// 4) Evaluation metrics structure
// =====================================================================================

struct EvalMetrics {
    double objective   = 0.0;   ///< Combined objective value
    double cm_error    = 0.0;   ///< Distance from target CM
    double sigma_proxy = 0.0;   ///< Bending stress surrogate
    double mass_est    = 0.0;   ///< Estimated remaining mass
    Point<DIM> cm{};            ///< Computed center of mass
};

// =====================================================================================
// 5) Monte Carlo evaluation using library estimators
// =====================================================================================

/**
 * @brief Evaluate objective using Monte Carlo integration.
 * @details Uses MCMeanEstimator to compute CM and stress proxy.
 *          Penalties are applied for constraint violations.
 */
static EvalMetrics evaluateMonteCarlo(const DroneArmDomain& domain,
                                       const Coordinates& params,
                                       int n_samples,
                                       uint32_t seed,
                                       double F_vertical,
                                       double L_arm,
                                       double arm_height,
                                       double w_cm,
                                       double w_sigma) {
    EvalMetrics out;

    // Validate axes
    const double a = params[3], b = params[4], c = params[5];
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    // Check center feasibility
    Point<DIM> center;
    center[0] = params[0];
    center[1] = params[1];
    center[2] = params[2];
    if (!domain.isInside(center)) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    // Precompute ellipsoid cache
    const EllipsoidCache E = buildEllipsoidCache(params);

    // Containment penalty via surface sampling
    RngManager rng_mgr(seed ^ 0xA53A9u);
    auto containment_rng = rng_mgr.make_rng(0);
    const double contain_frac = computeContainmentPenalty(domain, params, containment_rng, 48);
    double penalty = 5e5 * contain_frac * contain_frac;

    // Monte Carlo sampling for mass and moments
    RngManager sampler(seed);
    auto rng = sampler.make_rng(0);

    const auto bounds = domain.getBounds();
    std::uniform_real_distribution<double> dist_x(bounds[0].first, bounds[0].second);
    std::uniform_real_distribution<double> dist_y(bounds[1].first, bounds[1].second);
    std::uniform_real_distribution<double> dist_z(bounds[2].first, bounds[2].second);

    double hits = 0.0;
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double sum_z2 = 0.0;

    for (int i = 0; i < n_samples; ++i) {
        Point<DIM> p;
        p[0] = dist_x(rng);
        p[1] = dist_y(rng);
        p[2] = dist_z(rng);

        if (!domain.isInside(p)) continue;
        if (isInsideEllipsoid(p, E)) continue;

        hits += 1.0;
        sum_x += p[0];
        sum_y += p[1];
        sum_z += p[2];
        sum_z2 += p[2] * p[2];
    }

    if (hits <= 1e-9) {
        out.objective = 1e6 + penalty;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    // Compute center of mass
    out.cm[0] = sum_x / hits;
    out.cm[1] = sum_y / hits;
    out.cm[2] = sum_z / hits;

    // CM error from target (1, 0, 0)
    const double tx = 1.0, ty = 0.0, tz = 0.0;
    out.cm_error = std::sqrt(
        (out.cm[0] - tx) * (out.cm[0] - tx) +
        (out.cm[1] - ty) * (out.cm[1] - ty) +
        (out.cm[2] - tz) * (out.cm[2] - tz)
    );

    // Estimate remaining volume/mass
    const double box_vol = domain.getBoxVolume();
    out.mass_est = (hits / static_cast<double>(n_samples)) * box_vol;

    // Second-moment proxy: I_proxy ≈ V * Var(z)
    const double Ez = sum_z / hits;
    const double Ez2 = sum_z2 / hits;
    const double varz = std::max(0.0, Ez2 - Ez * Ez);
    const double I_proxy = out.mass_est * varz;

    // Stress proxy (cantilever beam): σ ~ (F*L*c) / I
    const double eps = 1e-12;
    const double Mmax = F_vertical * L_arm;
    const double c_fiber = 0.5 * arm_height;
    out.sigma_proxy = (Mmax * c_fiber) / (I_proxy + eps);

    // Mass constraint: penalize if too much material removed
    if (out.mass_est < 0.2) 
        penalty += 1e4;

    // Combined objective
    out.objective = w_cm * out.cm_error + w_sigma * out.sigma_proxy + penalty;
    return out;
}

// =====================================================================================
// 6) High-precision verification (parallel)
// =====================================================================================

/**
 * @brief High-precision verification using OpenMP parallelization.
 * @details Uses RngManager for deterministic per-thread seeding.
 */
static EvalMetrics verifyMonteCarlo(const DroneArmDomain& domain,
                                     const Coordinates& params,
                                     int n_samples,
                                     uint32_t seed,
                                     double F_vertical,
                                     double L_arm,
                                     double arm_height,
                                     double w_cm,
                                     double w_sigma) {
    EvalMetrics out;

    const double a = params[3], b = params[4], c = params[5];
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    Point<DIM> center;
    center[0] = params[0];
    center[1] = params[1];
    center[2] = params[2];
    if (!domain.isInside(center)) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    const EllipsoidCache E = buildEllipsoidCache(params);
    const auto bounds = domain.getBounds();
    const double box_vol = domain.getBoxVolume();

    long double hits = 0.0L;
    long double sum_x = 0.0L, sum_y = 0.0L, sum_z = 0.0L, sum_z2 = 0.0L;

    RngManager rng_mgr(seed + 9999u);

    #pragma omp parallel reduction(+:hits, sum_x, sum_y, sum_z, sum_z2)
    {
        const int tid = omp_get_thread_num();
        auto rng = rng_mgr.make_rng(tid);

        #pragma omp for
        for (int i = 0; i < n_samples; ++i) {
            Point<DIM> p;
            const double rx = std::generate_canonical<double, 32>(rng);
            const double ry = std::generate_canonical<double, 32>(rng);
            const double rz = std::generate_canonical<double, 32>(rng);
            
            p[0] = bounds[0].first + rx * (bounds[0].second - bounds[0].first);
            p[1] = bounds[1].first + ry * (bounds[1].second - bounds[1].first);
            p[2] = bounds[2].first + rz * (bounds[2].second - bounds[2].first);

            if (!domain.isInside(p)) continue;
            if (isInsideEllipsoid(p, E)) continue;

            hits += 1.0L;
            sum_x += p[0];
            sum_y += p[1];
            sum_z += p[2];
            sum_z2 += p[2] * p[2];
        }
    }

    if (hits <= 1e-9L) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    out.cm[0] = static_cast<double>(sum_x / hits);
    out.cm[1] = static_cast<double>(sum_y / hits);
    out.cm[2] = static_cast<double>(sum_z / hits);

    const double tx = 1.0, ty = 0.0, tz = 0.0;
    out.cm_error = std::sqrt(
        (out.cm[0] - tx) * (out.cm[0] - tx) +
        (out.cm[1] - ty) * (out.cm[1] - ty) +
        (out.cm[2] - tz) * (out.cm[2] - tz)
    );

    out.mass_est = (static_cast<double>(hits) / n_samples) * box_vol;

    const double Ez = static_cast<double>(sum_z / hits);
    const double Ez2 = static_cast<double>(sum_z2 / hits);
    const double varz = std::max(0.0, Ez2 - Ez * Ez);
    const double I_proxy = out.mass_est * varz;

    const double eps = 1e-12;
    const double Mmax = F_vertical * L_arm;
    const double c_fiber = 0.5 * arm_height;
    out.sigma_proxy = (Mmax * c_fiber) / (I_proxy + eps);

    out.objective = w_cm * out.cm_error + w_sigma * out.sigma_proxy;
    return out;
}

// =====================================================================================
// 7) Gnuplot visualization helper
// =====================================================================================

static void writeGnuplotScript(const std::string& filename,
                                const Solution& best_pso,
                                const Solution& best_ga,
                                const std::string& png_output = "") {
    std::ofstream gp(filename);
    if (!gp.is_open()) {
        throw std::runtime_error("Cannot write gnuplot script: " + filename);
    }

    auto get = [](const Solution& s, size_t i) { return s.params.at(i); };

    gp << "### Auto-generated: PSO vs GA ellipsoid hole comparison\n";
    gp << "set view 62, 32\n";
    gp << "set grid\n";
    gp << "set key outside\n";
    gp << "set xlabel 'X'\n";
    gp << "set ylabel 'Y'\n";
    gp << "set zlabel 'Z'\n";
    gp << "set xrange [-6:8]\n";
    gp << "set yrange [-3:3]\n";
    gp << "set zrange [-2:2]\n";

    if (!png_output.empty()) {
        gp << "set term pngcairo size 1400,1000\n";
        gp << "set output '" << png_output << "'\n";
    } else {
        gp << "set term qt size 1100,800\n";
    }

    gp << "set parametric\n";
    gp << "set urange [0:2*pi]\n";
    gp << "set vrange [0:pi]\n";
    gp << "set isosamples 45, 30\n";

    // PSO parameters
    gp << std::setprecision(17);
    gp << "hx1=" << get(best_pso, 0) << "; hy1=" << get(best_pso, 1) << "; hz1=" << get(best_pso, 2) << "\n";
    gp << "a1=" << get(best_pso, 3) << "; b1=" << get(best_pso, 4) << "; c1=" << get(best_pso, 5) << "\n";
    gp << "yaw1=" << get(best_pso, 6) << "; pitch1=" << get(best_pso, 7) << "; roll1=" << get(best_pso, 8) << "\n";

    // GA parameters
    gp << "hx2=" << get(best_ga, 0) << "; hy2=" << get(best_ga, 1) << "; hz2=" << get(best_ga, 2) << "\n";
    gp << "a2=" << get(best_ga, 3) << "; b2=" << get(best_ga, 4) << "; c2=" << get(best_ga, 5) << "\n";
    gp << "yaw2=" << get(best_ga, 6) << "; pitch2=" << get(best_ga, 7) << "; roll2=" << get(best_ga, 8) << "\n";

    // Ellipsoid parametric equations
    gp << "lx(a,u,v)=a*cos(u)*sin(v)\n";
    gp << "ly(b,u,v)=b*sin(u)*sin(v)\n";
    gp << "lz(c,v)=c*cos(v)\n";
    gp << "cy(t)=cos(t); sy(t)=sin(t)\n";
    gp << "cp(t)=cos(t); sp(t)=sin(t)\n";
    gp << "cr(t)=cos(t); sr(t)=sin(t)\n";

    // Rotation transformation (Rz*Ry*Rx)
    gp << "gx(x,y,z,yaw,pitch,roll)=(cy(yaw)*cp(pitch))*x + (cy(yaw)*sp(pitch)*sr(roll) - sy(yaw)*cr(roll))*y + (cy(yaw)*sp(pitch)*cr(roll) + sy(yaw)*sr(roll))*z\n";
    gp << "gy(x,y,z,yaw,pitch,roll)=(sy(yaw)*cp(pitch))*x + (sy(yaw)*sp(pitch)*sr(roll) + cy(yaw)*cr(roll))*y + (sy(yaw)*sp(pitch)*cr(roll) - cy(yaw)*sr(roll))*z\n";
    gp << "gz(x,y,z,yaw,pitch,roll)=(-sp(pitch))*x + (cp(pitch)*sr(roll))*y + (cp(pitch)*cr(roll))*z\n";

    gp << "X(hx,a,b,c,u,v,yaw,pitch,roll)=hx + gx(lx(a,u,v), ly(b,u,v), lz(c,v), yaw,pitch,roll)\n";
    gp << "Y(hy,a,b,c,u,v,yaw,pitch,roll)=hy + gy(lx(a,u,v), ly(b,u,v), lz(c,v), yaw,pitch,roll)\n";
    gp << "Z(hz,a,b,c,u,v,yaw,pitch,roll)=hz + gz(lx(a,u,v), ly(b,u,v), lz(c,v), yaw,pitch,roll)\n";

    gp << "splot \\\n";
    gp << "  X(hx1,a1,b1,c1,u,v,yaw1,pitch1,roll1), Y(hy1,a1,b1,c1,u,v,yaw1,pitch1,roll1), Z(hz1,a1,b1,c1,u,v,yaw1,pitch1,roll1) w lines lc rgb '#3366cc' title 'PSO Hole', \\\n";
    gp << "  X(hx2,a2,b2,c2,u,v,yaw2,pitch2,roll2), Y(hy2,a2,b2,c2,u,v,yaw2,pitch2,roll2), Z(hz2,a2,b2,c2,u,v,yaw2,pitch2,roll2) w lines lc rgb '#dc3912' title 'GA Hole'\n";

    if (png_output.empty()) {
        gp << "pause -1\n";
    }
}

// =====================================================================================
// 8) Result printing helper
// =====================================================================================

static void printSolutionBlock(const std::string& title,
                                const Coordinates& params,
                                const EvalMetrics& m) {
    std::cout << "\n" << title << "\n";
    std::cout << "  objective:   " << m.objective << "\n";
    std::cout << "  cm_error:    " << m.cm_error << "\n";
    std::cout << "  sigma_proxy: " << m.sigma_proxy << "\n";
    std::cout << "  mass_est:    " << m.mass_est << "\n";
    std::cout << "  CM:          (" << m.cm[0] << ", " << m.cm[1] << ", " << m.cm[2] << ")\n";
    std::cout << "  params:\n";
    std::cout << "    center = (" << params[0] << ", " << params[1] << ", " << params[2] << ")\n";
    std::cout << "    axes   = (" << params[3] << ", " << params[4] << ", " << params[5] << ")\n";
    std::cout << "    angles = (yaw=" << params[6] << ", pitch=" << params[7] << ", roll=" << params[8] << ")\n";
}

// =====================================================================================
// 9) Main
// =====================================================================================

int main(int argc, char* argv[]) {
    // Parse command line arguments
    int num_threads = omp_get_max_threads();

    if (argc > 1) {
        std::string seed_arg = argv[1];
        if (seed_arg != "-") {
            try { 
                GLOBAL_SEED = static_cast<uint32_t>(std::stoul(seed_arg)); 
            } catch (...) {}
        }
    }
    if (argc > 2) {
        try {
            int requested = std::stoi(argv[2]);
            num_threads = (requested <= 0) ? 1 : requested;
        } catch (...) {}
    }
    omp_set_num_threads(num_threads);

    std::cout << "\n=== Drone CM + Strength Optimization (Ellipsoid Hole) ===\n";
    std::cout << "=== PSO vs GA Comparison ===\n";
    std::cout << "Seed: " << GLOBAL_SEED << " | Threads: " << num_threads << "\n";
    std::cout << "(Usage: ./drone_ellipsoid [seed|-] [num_threads])\n\n";

    // Construct domain once
    DroneArmDomain domain;

    // Problem constants
    const double F_vertical = 1.0;
    const double L_arm = domain.arm_length;
    const double arm_h = domain.arm_height;

    // Scalarization weights
    const double w_cm = 1.0;
    const double w_sigma = 0.03;

    // MC budgets
    const int n_fast = 25000;
    const int n_verify = 1000000;

    // Parameter bounds: [hx, hy, hz, a, b, c, yaw, pitch, roll]
    const double PI = std::acos(-1.0);
    
    Coordinates lower = {-4.0, -1.0, -0.4, 0.10, 0.10, 0.10, -PI, -PI, -PI};
    Coordinates upper = {+4.0, +1.0, +0.4, 0.90, 0.80, 0.45, +PI, +PI, +PI};

    // Objective function with deterministic MC noise
    auto objective = [&](const Coordinates& params) -> Real {
        Coordinates p = params;
        p[3] = std::max(p[3], 1e-6);
        p[4] = std::max(p[4], 1e-6);
        p[5] = std::max(p[5], 1e-6);

        const uint32_t seed = hashParams(p) + GLOBAL_SEED;
        const EvalMetrics m = evaluateMonteCarlo(domain, p, n_fast, seed, 
                                                  F_vertical, L_arm, arm_h, w_cm, w_sigma);
        return m.objective;
    };

    // =========================================================================
    // PSO Optimization
    // =========================================================================
    std::cout << "--- Running PSO ---\n";

    PSOConfig pso_cfg;
    pso_cfg.population_size = 40;
    pso_cfg.max_iterations = 70;
    pso_cfg.cognitive_coeff = 1.6;
    pso_cfg.social_coeff = 1.6;
    pso_cfg.inertia_weight = 0.70;

    PSO pso(pso_cfg);
    pso.setObjectiveFunction(objective);
    pso.setBounds(lower, upper);
    pso.setMode(OptimizationMode::MINIMIZE);

    // Progress callback
    pso.setCallback([](const Solution& best, size_t iter) {
        if (iter % 10 == 0) {
            std::cout << "  PSO iter " << iter << " | obj: " << best.value << "\n";
        }
    });

    Solution best_pso = pso.optimize();

    // =========================================================================
    // GA Optimization
    // =========================================================================
    std::cout << "\n--- Running GA ---\n";

    GAConfig ga_cfg;
    ga_cfg.population_size = 100;
    ga_cfg.max_generations = 140;
    ga_cfg.tournament_k = 3;
    ga_cfg.crossover_rate = 0.9;
    ga_cfg.mutation_rate = 0.12;
    ga_cfg.mutation_sigma = 0.10;
    ga_cfg.elitism_count = 2;

    GA ga(ga_cfg);
    ga.setObjectiveFunction(objective);
    ga.setBounds(lower, upper);
    ga.setMode(OptimizationMode::MINIMIZE);

    // Progress callback
    ga.setCallback([](const Solution& best, size_t gen) {
        if (gen % 20 == 0) {
            std::cout << "  GA gen " << gen << " | obj: " << best.value << "\n";
        }
    });

    Solution best_ga = ga.optimize();

    // =========================================================================
    // Fast metrics report
    // =========================================================================
    std::cout << std::fixed << std::setprecision(8);

    auto evalFast = [&](const Solution& sol) -> EvalMetrics {
        const uint32_t seed = hashParams(sol.params) + GLOBAL_SEED;
        return evaluateMonteCarlo(domain, sol.params, n_fast, seed, 
                                   F_vertical, L_arm, arm_h, w_cm, w_sigma);
    };

    const EvalMetrics pso_fast = evalFast(best_pso);
    const EvalMetrics ga_fast = evalFast(best_ga);

    printSolutionBlock("PSO (fast MC metrics)", best_pso.params, pso_fast);
    printSolutionBlock("GA  (fast MC metrics)", best_ga.params, ga_fast);

    // =========================================================================
    // High-precision verification
    // =========================================================================
    std::cout << "\n--- HIGH-PRECISION VERIFICATION (" << n_verify << " samples) ---\n";

    auto evalVerify = [&](const Solution& sol) -> EvalMetrics {
        const uint32_t seed = hashParams(sol.params) + GLOBAL_SEED;
        return verifyMonteCarlo(domain, sol.params, n_verify, seed,
                                 F_vertical, L_arm, arm_h, w_cm, w_sigma);
    };

    const EvalMetrics pso_ver = evalVerify(best_pso);
    const EvalMetrics ga_ver = evalVerify(best_ga);

    printSolutionBlock("PSO (verification)", best_pso.params, pso_ver);
    printSolutionBlock("GA  (verification)", best_ga.params, ga_ver);

    // Comparison summary
    std::cout << "\n--- COMPARISON SUMMARY ---\n";
    std::cout << "PSO: obj=" << pso_ver.objective << " | cm_err=" << pso_ver.cm_error 
              << " | sigma=" << pso_ver.sigma_proxy << "\n";
    std::cout << "GA:  obj=" << ga_ver.objective << " | cm_err=" << ga_ver.cm_error 
              << " | sigma=" << ga_ver.sigma_proxy << "\n";

    const bool pso_wins = (pso_ver.objective < ga_ver.objective);
    std::cout << "\nBest optimizer: " << (pso_wins ? "PSO" : "GA") << "\n";

    // =========================================================================
    // Export geometry and visualization
    // =========================================================================
    const Coordinates& best_params = pso_wins ? best_pso.params : best_ga.params;

    std::cout << "\nExporting geometry for best solution...\n";
    domain.exportGeometryWithHole("./drone_frames", best_params);

    try {
        std::filesystem::create_directories("./drone_frames");

        const std::string gp_file = "./drone_frames/visualize_pso_vs_ga.gp";
        const std::string png_file = "./drone_frames/holes_comparison.png";

        // Generate PNG
        writeGnuplotScript(gp_file, best_pso, best_ga, png_file);
        std::cout << "Generated gnuplot script: " << gp_file << "\n";

        // Generate interactive script
        writeGnuplotScript(gp_file, best_pso, best_ga, "");
        std::system(("gnuplot -persist " + gp_file).c_str());

    } catch (const std::exception& e) {
        std::cerr << "Visualization failed: " << e.what() << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
