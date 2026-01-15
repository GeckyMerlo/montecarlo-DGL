/**
 * @file drone_optimization_ellipsoid_strength.cpp
 * @brief Drone arm CM targeting + vertical bending-stress proxy with an ellipsoidal hole.
 *
 * @details
 * This example is deliberately designed to benefit from your library stack:
 * - Complex 3D domain with fast inside-tests (CSG union): arm + motor + optional cabin polytope
 * - Stochastic objective via Monte Carlo integration (mass, CM, second-moment proxies)
 * - Deterministic parameter-based seeding to make MC noise stationary (crucial for PSO/GA)
 * - Non-smooth constraints handled by penalties (hole must remain inside the body, etc.)
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
 * Notes:
 * - This is NOT FEM; sigma_proxy is a formal, lightweight surrogate:
 *     sigma_proxy ~ (M_max * c) / I_proxy
 *   where M_max = F * L, c ~ arm_height/2, and I_proxy is estimated from MC via Var(z).
 *
 * - "Hole inside body" is enforced by a surface sampling check of the ellipsoid; violations are penalized.
 *
 * Outputs:
 * - Best solution from PSO and GA (fast metrics + high-precision verification)
 * - Optional geometry export for visualization (points outside the hole)
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

#include <omp.h>
#include <filesystem>
#include <fstream>

#include "../montecarlo/geometry.hpp"
#include "../montecarlo/integrators/montecarlo_integrator.hpp"
#include "../montecarlo/domains/hyperrectangle.hpp"
#include "../montecarlo/domains/hypercylinder.hpp"
#include "../montecarlo/domains/polytope.hpp"

#include "../montecarlo/optimizers/PSO.hpp"
#include "../montecarlo/optimizers/GA.hpp"

using namespace geom;
using namespace optimizers;

/// Global seed for deterministic behavior (also used by optimizers via extern if needed)
uint32_t GLOBAL_SEED = 12345;

/// Problem geometry dimension
constexpr size_t DIM = 3;

// =====================================================================================
// 0) Helpers: deterministic hashing, basic linear algebra, safe clamps
// =====================================================================================

static uint32_t hash_params(const std::vector<double>& params) {
    uint32_t seed = 0;
    std::hash<double> hasher;
    for (double p : params) {
        seed ^= static_cast<uint32_t>(hasher(p)) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
    }
    return seed;
}

static double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

struct Mat3 {
    double m[3][3]{};
};

static Mat3 mul(const Mat3& A, const Mat3& B) {
    Mat3 C;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            C.m[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) C.m[i][j] += A.m[i][k] * B.m[k][j];
        }
    return C;
}

static Mat3 transpose(const Mat3& A) {
    Mat3 T;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            T.m[i][j] = A.m[j][i];
    return T;
}

static std::array<double,3> mat_vec(const Mat3& A, const std::array<double,3>& v) {
    std::array<double,3> r{};
    for (int i = 0; i < 3; ++i) {
        r[i] = A.m[i][0]*v[0] + A.m[i][1]*v[1] + A.m[i][2]*v[2];
    }
    return r;
}

/**
 * @brief Build rotation matrix from yaw/pitch/roll (Z-Y-X convention).
 * yaw   = rotation around Z
 * pitch = rotation around Y
 * roll  = rotation around X
 */
static Mat3 rot_zyx(double yaw, double pitch, double roll) {
    const double cy = std::cos(yaw),   sy = std::sin(yaw);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double cr = std::cos(roll),  sr = std::sin(roll);

    Mat3 Rz{{ {cy, -sy, 0.0}, {sy, cy, 0.0}, {0.0, 0.0, 1.0} }};
    Mat3 Ry{{ {cp, 0.0, sp},  {0.0, 1.0, 0.0}, {-sp, 0.0, cp} }};
    Mat3 Rx{{ {1.0, 0.0, 0.0}, {0.0, cr, -sr}, {0.0, sr, cr} }};

    return mul(mul(Rz, Ry), Rx);
}

struct EllipsoidFast {
    double hx, hy, hz;           // center
    double inv_a2, inv_b2, inv_c2; // 1/a^2, 1/b^2, 1/c^2
    Mat3 Rt;                     // R^T
};

// Build once per evaluation
static inline EllipsoidFast makeEllipsoidFast(const std::vector<double>& params) {
    EllipsoidFast E;
    E.hx = params[0]; E.hy = params[1]; E.hz = params[2];

    const double a = params[3], b = params[4], c = params[5];
    // (Assume checked > 0 elsewhere)
    E.inv_a2 = 1.0 / (a*a);
    E.inv_b2 = 1.0 / (b*b);
    E.inv_c2 = 1.0 / (c*c);

    const double yaw = params[6], pitch = params[7], roll = params[8];
    E.Rt = transpose(rot_zyx(yaw, pitch, roll));
    return E;
}

static inline bool isInsideEllipsoidFast(const Point<3>& p, const EllipsoidFast& E) {
    const double dx = p[0] - E.hx;
    const double dy = p[1] - E.hy;
    const double dz = p[2] - E.hz;

    // q = Rt * d
    const double qx = E.Rt.m[0][0]*dx + E.Rt.m[0][1]*dy + E.Rt.m[0][2]*dz;
    const double qy = E.Rt.m[1][0]*dx + E.Rt.m[1][1]*dy + E.Rt.m[1][2]*dz;
    const double qz = E.Rt.m[2][0]*dx + E.Rt.m[2][1]*dy + E.Rt.m[2][2]*dz;

    const double u = (qx*qx)*E.inv_a2 + (qy*qy)*E.inv_b2 + (qz*qz)*E.inv_c2;
    return u <= 1.0;
}

// =====================================================================================
// 1) Optional PolyTope readers (same structure as your original example)
// =====================================================================================

template <int dim>
static std::vector<geom::Point<dim>> read_points_from_file_drone(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open()) throw std::runtime_error("Cannot open file: " + filename);

    std::size_t num_points = 0;
    std::size_t file_dim   = 0;
    in >> num_points >> file_dim;
    if (!in.good()) throw std::runtime_error("Error reading header from file: " + filename);

    if (file_dim != static_cast<std::size_t>(dim)) {
        throw std::runtime_error(
            "Dimension mismatch: file has dim = " + std::to_string(file_dim) +
            " but template expects dim = " + std::to_string(dim));
    }

    std::vector<geom::Point<dim>> points;
    points.reserve(num_points);

    for (std::size_t i = 0; i < num_points; ++i) {
        geom::Point<dim> p;
        for (int k = 0; k < dim; ++k) {
            if (!(in >> p[k])) {
                throw std::runtime_error(
                    "Error reading coordinate " + std::to_string(k) +
                    " of point " + std::to_string(i) + " from file: " + filename);
            }
        }
        points.push_back(p);
    }
    return points;
}

template <int dim>
static void read_normals_and_offsets_from_file_drone(
    const std::string& filename,
    std::vector<std::array<double, dim>>& normals,
    std::vector<double>& offsets)
{
    std::ifstream in(filename);
    if (!in.is_open()) throw std::runtime_error("Cannot open normals file: " + filename);

    std::size_t file_dim = 0;
    std::size_t num_facets = 0;
    in >> file_dim >> num_facets;
    if (!in.good()) throw std::runtime_error("Error reading header (dim, num_facets) from: " + filename);

    if (file_dim != static_cast<std::size_t>(dim + 1)) {
        throw std::runtime_error(
            "Dimension mismatch in normals file: file has dim = " + std::to_string(file_dim) +
            " but template expects dim+1 = " + std::to_string(dim + 1));
    }

    normals.clear();
    offsets.clear();
    normals.reserve(num_facets);
    offsets.reserve(num_facets);

    for (std::size_t f = 0; f < num_facets; ++f) {
        std::array<double, dim> n{};
        double d = 0.0;

        for (std::size_t k = 0; k < static_cast<std::size_t>(dim); ++k) {
            if (!(in >> n[k])) {
                throw std::runtime_error("Error reading normal component from: " + filename);
            }
        }
        if (!(in >> d)) throw std::runtime_error("Error reading offset d from: " + filename);

        // If file stores plane as n·x + d <= 0, then b = -d for n·x <= b
        normals.push_back(n);
        offsets.push_back(-d);
    }
}

// =====================================================================================
// 2) Domain definition (CSG union): Arm + Motor + optional Cabin (PolyTope)
// =====================================================================================

class DroneArmDomain : public IntegrationDomain<DIM> {
public:
    std::unique_ptr<HyperRectangle<DIM>> arm;
    std::unique_ptr<HyperCylinder<DIM>> motor_housing;
    std::unique_ptr<PolyTope<DIM>> cabin;

    std::array<double, DIM> motor_offset{};
    std::array<double, DIM> cabin_offset{};

    // Arm dimensions (used by stress surrogate)
    double arm_length = 10.0; // along X
    double arm_width  = 2.0;  // along Y
    double arm_height = 1.0;  // along Z

    DroneArmDomain() {
        // Arm centered at origin (assumption consistent with your original comments)
        std::array<double, DIM> arm_dims = {arm_length, arm_width, arm_height};
        arm = std::make_unique<HyperRectangle<DIM>>(arm_dims);

        // Motor housing: cylinder (radius 1.5, height 1.2 along Z)
        motor_housing = std::make_unique<HyperCylinder<DIM>>(1.5, 1.2);
        motor_offset = {+arm_length/2.0, 0.0, -0.6};

        // Optional cabin
        try {
            auto points = read_points_from_file_drone<DIM>("../drone_assets/cabin_points.txt");
            std::vector<std::array<double, DIM>> normals;
            std::vector<double> offsets;
            read_normals_and_offsets_from_file_drone<DIM>("../drone_assets/cabin_hull.txt", normals, offsets);

            cabin = std::make_unique<PolyTope<DIM>>(points, normals, offsets);
            cabin_offset = {-2.0, 0.0, 0.5};
        } catch (...) {
            cabin = nullptr;
            cabin_offset = {0.0, 0.0, 0.0};
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
        if (arm->isInside(p)) return true;

        Point<DIM> p_motor;
        for (size_t i = 0; i < DIM; ++i) p_motor[i] = p[i] - motor_offset[i];
        if (motor_housing->isInside(p_motor)) return true;

        if (cabin) {
            Point<DIM> p_cabin;
            for (size_t i = 0; i < DIM; ++i) p_cabin[i] = p[i] - cabin_offset[i];
            if (cabin->isInside(p_cabin)) return true;
        }

        return false;
    }

    /**
     * @brief Export sampled geometry points excluding an ellipsoidal hole (for visualization).
     * This is intentionally simple and not used for physics, only for plotting.
     */
    void exportGeometryEllipsoid(const std::string& out_dir,
                                const std::vector<double>& params) const
    {
        try { std::filesystem::create_directories(out_dir); } catch (...) {}

        const std::string filename = out_dir + "/drone_domain_ellipsoid.txt";
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "Cannot open export file: " << filename << "\n";
            return;
        }

        const double hx = params[0], hy = params[1], hz = params[2];
        const double a = params[3],  b = params[4],  c = params[5];
        const double yaw = params[6], pitch = params[7], roll = params[8];

        const Mat3 R = rot_zyx(yaw, pitch, roll);
        const Mat3 Rt = transpose(R);

        auto is_in_hole = [&](double x, double y, double z) -> bool {
            std::array<double,3> d{ x - hx, y - hy, z - hz };
            const auto q = mat_vec(Rt, d);
            const double u = (q[0]/a)*(q[0]/a) + (q[1]/b)*(q[1]/b) + (q[2]/c)*(q[2]/c);
            return u <= 1.0;
        };

        out << "# type x y z\n";
        out << "# Ellipsoid hole: center=(" << hx << "," << hy << "," << hz << ") "
            << "axes=(" << a << "," << b << "," << c << ") "
            << "angles(yaw,pitch,roll)=(" << yaw << "," << pitch << "," << roll << ")\n";

        // Sample arm region (fine grid)
        const double step = 0.15;
        for (double x = -arm_length/2.0; x <= arm_length/2.0; x += step) {
            for (double y = -arm_width/2.0; y <= arm_width/2.0; y += step) {
                for (double z = -arm_height/2.0; z <= arm_height/2.0; z += step) {
                    Point<DIM> p; p[0]=x; p[1]=y; p[2]=z;
                    if (arm->isInside(p) && !is_in_hole(x,y,z)) out << "arm " << x << " " << y << " " << z << "\n";
                }
            }
        }

        // Sample motor region (coarse bounding)
        const double step_motor = 0.15;
        for (double x = (motor_offset[0] - 2.0); x <= (motor_offset[0] + 2.0); x += step_motor) {
            for (double y = -2.0; y <= 2.0; y += step_motor) {
                for (double z = -1.5; z <= 1.5; z += step_motor) {
                    Point<DIM> p; p[0]=x; p[1]=y; p[2]=z;
                    Point<DIM> pl; for (size_t i=0;i<DIM;++i) pl[i]=p[i]-motor_offset[i];
                    if (motor_housing->isInside(pl) && !is_in_hole(x,y,z)) out << "motor " << x << " " << y << " " << z << "\n";
                }
            }
        }

        std::cout << "Export geometry: " << filename << "\n";
    }
};

// =====================================================================================
// 3) Ellipsoid hole test + containment check (surface sampling)
// =====================================================================================

static bool isInsideEllipsoid(const Point<DIM>& p,
                              const std::vector<double>& params,
                              const Mat3& Rt)
{
    const double hx = params[0], hy = params[1], hz = params[2];
    const double a = params[3],  b = params[4],  c = params[5];

    std::array<double,3> d{ p[0]-hx, p[1]-hy, p[2]-hz };
    const auto q = mat_vec(Rt, d);

    const double u = (q[0]/a)*(q[0]/a) + (q[1]/b)*(q[1]/b) + (q[2]/c)*(q[2]/c);
    return u <= 1.0;
}

/**
 * @brief Surface containment check: sample points on ellipsoid surface; ensure all are inside the domain.
 * This enforces the "hole entirely inside body" constraint approximately but robustly.
 */
static double holeContainmentPenaltyFast(const DroneArmDomain& domain,
                                         const std::vector<double>& params,
                                         uint32_t seed,
                                         int n_surface_samples)
{
    // Quick fail if center not inside
    Point<3> hc; hc[0]=params[0]; hc[1]=params[1]; hc[2]=params[2];
    if (!domain.isInside(hc)) return 1.0;

    const double a = params[3], b = params[4], c = params[5];
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) return 1.0;

    // Build rotation once (global from local)
    const Mat3 R = rot_zyx(params[6], params[7], params[8]);

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u11(-1.0, 1.0);

    int violations = 0;

    for (int i = 0; i < n_surface_samples; ++i) {
        // Marsaglia method: sample point uniformly on sphere without trig
        double u, v, s;
        do {
            u = u11(rng);
            v = u11(rng);
            s = u*u + v*v;
        } while (s >= 1.0 || s < 1e-12);

        const double factor = std::sqrt(1.0 - s);

        // Unit direction on sphere
        const double sx = 2.0*u*factor;
        const double sy = 2.0*v*factor;
        const double sz = 1.0 - 2.0*s;

        // Point on ellipsoid surface in local coords
        std::array<double,3> pl{ a*sx, b*sy, c*sz };

        // Rotate to global coords: pg = R * pl
        const auto pg = mat_vec(R, pl);

        Point<3> p;
        p[0] = params[0] + pg[0];
        p[1] = params[1] + pg[1];
        p[2] = params[2] + pg[2];

        if (!domain.isInside(p)) violations++;
    }

    return static_cast<double>(violations) / static_cast<double>(n_surface_samples);
}

// =====================================================================================
// 4) Objective evaluation via Monte Carlo (CM + stiffness proxy + penalties)
// =====================================================================================

struct EvalMetrics {
    double objective   = 0.0;
    double cm_error    = 0.0;
    double sigma_proxy = 0.0;
    double mass_est    = 0.0;
    Point<DIM> cm{};
};

/**
 * @brief Evaluate candidate with MC:
 * - remaining material = domain AND NOT ellipsoid
 * - compute CM and a second-moment proxy along z (Var(z)*V) => I_proxy
 * - sigma_proxy ~ (F*L*c)/I_proxy
 * - penalties: invalid axes, hole not contained, too low remaining mass, etc.
 */
static EvalMetrics evaluate_mc(const DroneArmDomain& domain,
                              const std::vector<double>& params,
                              int n_samples,
                              uint32_t seed,
                              double F_vertical,
                              double L_arm,
                              double arm_height,
                              double w_cm,
                              double w_sigma)
{
    EvalMetrics out;

    // Unpack
    const double hx = params[0], hy = params[1], hz = params[2];
    const double a  = params[3], b  = params[4], c  = params[5];
    const double yaw = params[6], pitch = params[7], roll = params[8];

    // Axis feasibility
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    // Basic feasibility: center inside
    Point<DIM> hc; hc[0]=hx; hc[1]=hy; hc[2]=hz;
    if (!domain.isInside(hc)) {
        out.objective = 1e6;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    const EllipsoidFast E = makeEllipsoidFast(params);

    // Containment penalty (surface sampling)
    const double contain_frac = holeContainmentPenaltyFast(domain, params, seed ^ 0xA53A9u, 48);
    // Strong penalty if not contained
    double penalty = 0.0;
    penalty += 5e5 * contain_frac * contain_frac;

    // MC sampling
    std::mt19937 rng(seed);

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
        if (isInsideEllipsoidFast(p, E)) continue;

        hits += 1.0;
        sum_x += p[0];
        sum_y += p[1];
        sum_z += p[2];
        sum_z2 += p[2]*p[2];
    }

    if (hits <= 1e-9) {
        out.objective = 1e6 + penalty;
        out.cm_error = 1e6;
        out.sigma_proxy = 1e6;
        return out;
    }

    const double cmx = sum_x / hits;
    const double cmy = sum_y / hits;
    const double cmz = sum_z / hits;

    out.cm[0]=cmx; out.cm[1]=cmy; out.cm[2]=cmz;

    // CM target
    const double tx = 1.0, ty = 0.0, tz = 0.0;
    out.cm_error = std::sqrt((cmx-tx)*(cmx-tx) + (cmy-ty)*(cmy-ty) + (cmz-tz)*(cmz-tz));

    // Estimate remaining volume/mass by hit ratio
    const double box_vol = domain.getBoxVolume();
    out.mass_est = (hits / static_cast<double>(n_samples)) * box_vol;

    // Second-moment proxy along z: I_proxy ≈ ∫(z - zbar)^2 dV = V * Var(z)
    const double Ez  = sum_z / hits;
    const double Ez2 = sum_z2 / hits;
    const double varz = std::max(0.0, Ez2 - Ez*Ez);
    const double I_proxy = out.mass_est * varz;

    // Stress proxy (cantilever):
    // sigma ~ (F*L*c)/I, c ~ arm_height/2
    const double eps = 1e-12;
    const double Mmax = F_vertical * L_arm;
    const double c_fiber = 0.5 * arm_height;
    out.sigma_proxy = (Mmax * c_fiber) / (I_proxy + eps);

    // Additional simple constraints:
    // - avoid "almost deleting everything" (soft)
    if (out.mass_est < 0.2) penalty += 1e4;

    // Scalar objective
    out.objective = w_cm * out.cm_error + w_sigma * out.sigma_proxy + penalty;
    return out;
}

// High-precision verification (optionally OpenMP)
static EvalMetrics verify_mc(const DroneArmDomain& domain,
                            const std::vector<double>& params,
                            int n_samples,
                            uint32_t seed,
                            double F_vertical,
                            double L_arm,
                            double arm_height,
                            double w_cm,
                            double w_sigma,
                            bool use_openmp)
{
    EvalMetrics out;

    const double a = params[3], b = params[4], c = params[5];
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) {
        out.objective = 1e6; out.cm_error = 1e6; out.sigma_proxy = 1e6;
        return out;
    }

    // Quick center feasibility
    Point<DIM> hc; hc[0]=params[0]; hc[1]=params[1]; hc[2]=params[2];
    if (!domain.isInside(hc)) {
        out.objective = 1e6; out.cm_error = 1e6; out.sigma_proxy = 1e6;
        return out;
    }

    const EllipsoidFast E = makeEllipsoidFast(params);

    const auto bounds = domain.getBounds();
    const double box_vol = domain.getBoxVolume();

    long double hits = 0.0L;
    long double sum_x = 0.0L, sum_y = 0.0L, sum_z = 0.0L, sum_z2 = 0.0L;

    if (use_openmp) {
        #pragma omp parallel reduction(+:hits,sum_x,sum_y,sum_z,sum_z2)
        {
            const int tid = omp_get_thread_num();
            std::mt19937 rng(seed + 9999u + static_cast<uint32_t>(tid));

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
                if (isInsideEllipsoidFast(p, E)) continue;

                hits += 1.0L;
                sum_x += p[0];
                sum_y += p[1];
                sum_z += p[2];
                sum_z2 += p[2]*p[2];
            }
        }
    } else {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> dist_x(bounds[0].first, bounds[0].second);
        std::uniform_real_distribution<double> dist_y(bounds[1].first, bounds[1].second);
        std::uniform_real_distribution<double> dist_z(bounds[2].first, bounds[2].second);

        for (int i = 0; i < n_samples; ++i) {
            Point<DIM> p;
            p[0] = dist_x(rng);
            p[1] = dist_y(rng);
            p[2] = dist_z(rng);

            if (!domain.isInside(p)) continue;
            if (isInsideEllipsoidFast(p, E)) continue;

            hits += 1.0L;
            sum_x += p[0];
            sum_y += p[1];
            sum_z += p[2];
            sum_z2 += p[2]*p[2];
        }
    }

    if (hits <= 1e-9L) {
        out.objective = 1e6; out.cm_error = 1e6; out.sigma_proxy = 1e6;
        return out;
    }

    const double H = static_cast<double>(hits);
    const double cmx = static_cast<double>(sum_x / hits);
    const double cmy = static_cast<double>(sum_y / hits);
    const double cmz = static_cast<double>(sum_z / hits);
    out.cm[0]=cmx; out.cm[1]=cmy; out.cm[2]=cmz;

    const double tx=1.0, ty=0.0, tz=0.0;
    out.cm_error = std::sqrt((cmx-tx)*(cmx-tx) + (cmy-ty)*(cmy-ty) + (cmz-tz)*(cmz-tz));

    out.mass_est = (H / static_cast<double>(n_samples)) * box_vol;

    const double Ez  = static_cast<double>(sum_z / hits);
    const double Ez2 = static_cast<double>(sum_z2 / hits);
    const double varz = std::max(0.0, Ez2 - Ez*Ez);
    const double I_proxy = out.mass_est * varz;

    const double eps = 1e-12;
    const double Mmax = F_vertical * L_arm;
    const double c_fiber = 0.5 * arm_height;
    out.sigma_proxy = (Mmax * c_fiber) / (I_proxy + eps);

    out.objective = w_cm * out.cm_error + w_sigma * out.sigma_proxy;
    return out;
}

// =====================================================================================
//  5) PLOTTING UTILITIES (EXPORT GEOMETRY WITH ELLIPSOID HOLE)    
// =====================================================================================

static void writeGnuplotTwoEllipsoidsScript(
    const std::string& gp_filename,
    const Solution& best_pso,
    const Solution& best_ga,
    const std::string& out_png = ""  // se vuoto -> finestra interattiva
) {
    auto get = [](const Solution& s, size_t i) -> double { return s.params.at(i); };

    // PSO params layout: [hx,hy,hz, a,b,c, yaw,pitch,roll]
    const double hx1 = get(best_pso, 0), hy1 = get(best_pso, 1), hz1 = get(best_pso, 2);
    const double a1  = get(best_pso, 3), b1  = get(best_pso, 4), c1  = get(best_pso, 5);
    const double yaw1 = get(best_pso, 6), pitch1 = get(best_pso, 7), roll1 = get(best_pso, 8);

    const double hx2 = get(best_ga, 0), hy2 = get(best_ga, 1), hz2 = get(best_ga, 2);
    const double a2  = get(best_ga, 3), b2  = get(best_ga, 4), c2  = get(best_ga, 5);
    const double yaw2 = get(best_ga, 6), pitch2 = get(best_ga, 7), roll2 = get(best_ga, 8);

    std::ofstream gp(gp_filename);
    if (!gp.is_open()) {
        throw std::runtime_error("Cannot write gnuplot script: " + gp_filename);
    }

    gp << "### Auto-generated: visualize PSO vs GA ellipsoid holes\n";
    gp << "set view 62, 32\n";
    gp << "set ticslevel 0\n";
    gp << "set grid\n";
    gp << "set key outside\n";
    gp << "set xlabel 'X'\n";
    gp << "set ylabel 'Y'\n";
    gp << "set zlabel 'Z'\n";
    gp << "set size ratio -1\n";
    gp << "set xrange [-6:8]\n";
    gp << "set yrange [-3:3]\n";
    gp << "set zrange [-2:2]\n";

    if (!out_png.empty()) {
        gp << "set term pngcairo size 1400,1000\n";
        gp << "set output '" << out_png << "'\n";
    } else {
        gp << "set term qt size 1100,800\n";
    }

    gp << "set parametric\n";
    gp << "set urange [0:2*pi]\n";
    gp << "set vrange [0:pi]\n";
    gp << "set isosamples 45, 30\n";

    // Inject numeric params (fixed precision)
    gp << "hx1=" << std::setprecision(17) << hx1 << "\n";
    gp << "hy1=" << std::setprecision(17) << hy1 << "\n";
    gp << "hz1=" << std::setprecision(17) << hz1 << "\n";
    gp << "a1="  << std::setprecision(17) << a1  << "\n";
    gp << "b1="  << std::setprecision(17) << b1  << "\n";
    gp << "c1="  << std::setprecision(17) << c1  << "\n";
    gp << "yaw1="   << std::setprecision(17) << yaw1   << "\n";
    gp << "pitch1=" << std::setprecision(17) << pitch1 << "\n";
    gp << "roll1="  << std::setprecision(17) << roll1  << "\n";

    gp << "hx2=" << std::setprecision(17) << hx2 << "\n";
    gp << "hy2=" << std::setprecision(17) << hy2 << "\n";
    gp << "hz2=" << std::setprecision(17) << hz2 << "\n";
    gp << "a2="  << std::setprecision(17) << a2  << "\n";
    gp << "b2="  << std::setprecision(17) << b2  << "\n";
    gp << "c2="  << std::setprecision(17) << c2  << "\n";
    gp << "yaw2="   << std::setprecision(17) << yaw2   << "\n";
    gp << "pitch2=" << std::setprecision(17) << pitch2 << "\n";
    gp << "roll2="  << std::setprecision(17) << roll2  << "\n";

    // Math
    gp << "lx(a,u,v)=a*cos(u)*sin(v)\n";
    gp << "ly(b,u,v)=b*sin(u)*sin(v)\n";
    gp << "lz(c,v)=c*cos(v)\n";

    gp << "cy(t)=cos(t); sy(t)=sin(t)\n";
    gp << "cp(t)=cos(t); sp(t)=sin(t)\n";
    gp << "cr(t)=cos(t); sr(t)=sin(t)\n";

    // Rz*Ry*Rx applied to (x,y,z)
    gp << "gx(x,y,z,yaw,pitch,roll)=(cy(yaw)*cp(pitch))*x + (cy(yaw)*sp(pitch)*sr(roll) - sy(yaw)*cr(roll))*y + (cy(yaw)*sp(pitch)*cr(roll) + sy(yaw)*sr(roll))*z\n";
    gp << "gy(x,y,z,yaw,pitch,roll)=(sy(yaw)*cp(pitch))*x + (sy(yaw)*sp(pitch)*sr(roll) + cy(yaw)*cr(roll))*y + (sy(yaw)*sp(pitch)*cr(roll) - cy(yaw)*sr(roll))*z\n";
    gp << "gz(x,y,z,yaw,pitch,roll)=(-sp(pitch))*x + (cp(pitch)*sr(roll))*y + (cp(pitch)*cr(roll))*z\n";

    gp << "X(hx,a,b,c,u,v,yaw,pitch,roll)=hx + gx(lx(a,u,v), ly(b,u,v), lz(c,v), yaw,pitch,roll)\n";
    gp << "Y(hy,a,b,c,u,v,yaw,pitch,roll)=hy + gy(lx(a,u,v), ly(b,u,v), lz(c,v), yaw,pitch,roll)\n";
    gp << "Z(hz,a,b,c,u,v,yaw,pitch,roll)=hz + gz(lx(a,u,v), ly(b,u,v), lz(c,v), yaw,pitch,roll)\n";

    gp << "splot \\\n";
    gp << "  X(hx1,a1,b1,c1,u,v,yaw1,pitch1,roll1), Y(hy1,a1,b1,c1,u,v,yaw1,pitch1,roll1), Z(hz1,a1,b1,c1,u,v,yaw1,pitch1,roll1) w lines title 'Hole PSO', \\\n";
    gp << "  X(hx2,a2,b2,c2,u,v,yaw2,pitch2,roll2), Y(hy2,a2,b2,c2,u,v,yaw2,pitch2,roll2), Z(hz2,a2,b2,c2,u,v,yaw2,pitch2,roll2) w lines title 'Hole GA'\n";

    if (out_png.empty()) {
        gp << "pause -1\n";
    } else {
        gp << "unset output\n";
    }
}

// =====================================================================================
// 6) Main: baseline MC (optional), PSO vs GA, verification, export
// =====================================================================================

static void print_solution_block(const std::string& title,
                                 const std::vector<double>& params,
                                 const EvalMetrics& m)
{
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

int main(int argc, char* argv[]) {
    // Usage: ./drone_ellipsoid [seed|-] [num_threads]
    int num_threads = omp_get_max_threads();

    if (argc > 1) {
        std::string seed_arg = argv[1];
        if (seed_arg != "-") {
            try { GLOBAL_SEED = static_cast<uint32_t>(std::stoul(seed_arg)); } catch (...) {}
        }
    }
    if (argc > 2) {
        try {
            const int requested = std::stoi(argv[2]);
            num_threads = (requested <= 0) ? 1 : requested;
        } catch (...) {}
    }
    omp_set_num_threads(num_threads);

    std::cout << "\n--- Drone CM + Strength (Ellipsoid Hole) | PSO vs GA ---\n";
    std::cout << "Seed: " << GLOBAL_SEED << " | Threads: " << num_threads << "\n";
    std::cout << "(Usage: ./drone_ellipsoid [seed|-] [num_threads])\n\n";

    // Construct domain ONCE (important: avoid I/O in each evaluation)
    DroneArmDomain domain;

    // --- Problem constants (surrogate model)
    const double F_vertical = 1.0;                 // load proxy
    const double L_arm = domain.arm_length;        // cantilever length
    const double arm_h = domain.arm_height;        // for c = h/2

    // --- Scalarization weights
    // Make it "meaningfully multi-criteria": CM is primary; stiffness is relevant but secondary.
    const double w_cm = 1.0;
    const double w_sigma = 0.03;

    // --- MC budgets
    const int n_fast   = 25000;   // optimization budget
    const int n_verify = 1000000; // verification

    // --- Bounds for params = [hx,hy,hz,a,b,c,yaw,pitch,roll]
    // Center bounds: inside approximate arm extent
    // Axes bounds: enforce smaller than arm thickness/width to make containment feasible
    // Angles bounds: [-pi, pi]
    const double PI = std::acos(-1.0);

    std::vector<double> lower = {
        -4.0, -1.0, -0.4,
         0.10, 0.10, 0.10,
        -PI,  -PI,  -PI
    };
    std::vector<double> upper = {
        +4.0, +1.0, +0.4,
         0.90, 0.80, 0.45,
        +PI,  +PI,  +PI
    };

    // Objective function (stationary MC noise via hash(params)+GLOBAL_SEED)
    auto objective = [&](const std::vector<double>& params) -> double {
        // Additional safety: avoid extremely thin or enormous axes (repair-like soft)
        // (GA/PSO already keep within bounds, but this stabilizes near-boundary behavior.)
        std::vector<double> p = params;
        p[3] = std::max(p[3], 1e-6);
        p[4] = std::max(p[4], 1e-6);
        p[5] = std::max(p[5], 1e-6);

        const uint32_t seed = hash_params(p) + GLOBAL_SEED;
        const EvalMetrics m = evaluate_mc(domain, p, n_fast, seed, F_vertical, L_arm, arm_h, w_cm, w_sigma);
        return m.objective;
    };

    // =================================================================================
    // PSO
    // =================================================================================
    PSOConfig pso_cfg;
    pso_cfg.population_size = 40;
    pso_cfg.max_iterations  = 70;
    pso_cfg.cognitive_coeff = 1.6;
    pso_cfg.social_coeff    = 1.6;
    pso_cfg.inertia_weight  = 0.70;

    PSO pso(pso_cfg);
    pso.setObjectiveFunction(objective);
    pso.setBounds(lower, upper);
    pso.setMode(OptimizationMode::MINIMIZE);

    std::cout << "--- Running PSO ---\n";
    Solution best_pso = pso.optimize();

    // =================================================================================
    // GA
    // =================================================================================
    GAConfig ga_cfg;
    ga_cfg.population_size = 100;
    ga_cfg.max_generations = 140;
    ga_cfg.tournament_k = 3;
    ga_cfg.crossover_rate = 0.9;
    ga_cfg.mutation_rate  = 0.12;
    ga_cfg.mutation_sigma = 0.10;
    ga_cfg.elitism_count  = 2;

    GA ga(ga_cfg);
    ga.setObjectiveFunction(objective);
    ga.setBounds(lower, upper);
    ga.setMode(OptimizationMode::MINIMIZE);

    std::cout << "\n--- Running GA ---\n";
    Solution best_ga = ga.optimize();

    // =================================================================================
    // FAST metrics report (same MC estimator used in objective)
    // =================================================================================
    std::cout << std::fixed << std::setprecision(8);

    auto eval_fast = [&](const Solution& sol) -> EvalMetrics {
        std::vector<double> p = sol.params; // assumes your Solution stores "params"
        const uint32_t seed = hash_params(p) + GLOBAL_SEED;
        return evaluate_mc(domain, p, n_fast, seed, F_vertical, L_arm, arm_h, w_cm, w_sigma);
    };

    const EvalMetrics pso_fast = eval_fast(best_pso);
    const EvalMetrics ga_fast  = eval_fast(best_ga);

    print_solution_block("PSO (fast MC metrics)", best_pso.params, pso_fast);
    print_solution_block("GA  (fast MC metrics)", best_ga.params,  ga_fast);

    // =================================================================================
    // HIGH-PRECISION verification (1M samples), correlated per-thread RNG
    // =================================================================================
    std::cout << "\n--- HIGH-PRECISION VERIFICATION (" << n_verify << " samples) ---\n";

    auto eval_verify = [&](const Solution& sol) -> EvalMetrics {
        std::vector<double> p = sol.params;
        const uint32_t seed = hash_params(p) + GLOBAL_SEED;
        return verify_mc(domain, p, n_verify, seed, F_vertical, L_arm, arm_h, w_cm, w_sigma, /*use_openmp=*/true);
    };

    const EvalMetrics pso_ver = eval_verify(best_pso);
    const EvalMetrics ga_ver  = eval_verify(best_ga);

    print_solution_block("PSO (verification metrics)", best_pso.params, pso_ver);
    print_solution_block("GA  (verification metrics)", best_ga.params,  ga_ver);

    // Quick comparison summary
    std::cout << "\n--- COMPARISON SUMMARY ---\n";
    std::cout << "PSO verified objective: " << pso_ver.objective
              << " | cm_error: " << pso_ver.cm_error
              << " | sigma_proxy: " << pso_ver.sigma_proxy << "\n";
    std::cout << "GA  verified objective: " << ga_ver.objective
              << " | cm_error: " << ga_ver.cm_error
              << " | sigma_proxy: " << ga_ver.sigma_proxy << "\n";

    // =================================================================================
    // Export geometry for the best verified (lower objective)
    // =================================================================================
    const bool pso_better = (pso_ver.objective < ga_ver.objective);
    const std::vector<double>& best_params = pso_better ? best_pso.params : best_ga.params;

    std::cout << "\nExporting geometry for best solution (" << (pso_better ? "PSO" : "GA") << ")...\n";
    domain.exportGeometryEllipsoid("./drone_frames", best_params);

    std::cout << "\nDone.\n";
    try {
        std::filesystem::create_directories("./drone_frames");

        const std::string gp_file  = "./drone_frames/visualize_holes_pso_vs_ga.gp";
        const std::string png_file = "./drone_frames/holes_pso_vs_ga.png"; // se vuoi PNG automatico

        // 1) genera script + PNG
        writeGnuplotTwoEllipsoidsScript(gp_file, best_pso, best_ga, png_file);
        std::cout << "\nGenerated gnuplot script: " << gp_file << "\n";
        std::cout << "Generated PNG image:      " << png_file << "\n";

        writeGnuplotTwoEllipsoidsScript(gp_file, best_pso, best_ga, "");
        std::system(("gnuplot -persist " + gp_file).c_str());
        
    } catch (const std::exception& e) {
        std::cerr << "Plot generation failed: " << e.what() << "\n";
    }

    return 0;
}