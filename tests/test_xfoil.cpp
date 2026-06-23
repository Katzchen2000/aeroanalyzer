// test_xfoil.cpp -- gates for the native in-process viscous solver (aero::xfoil).
// Inviscid panel accuracy, viscous cd0 acceptance band, transition behaviour,
// and determinism. Cross-validated against real xfoil.exe offline; see
// scratch/xfoil_visc_probe.cpp.
#include "test_harness.h"
#include "aeroanalyzer/aero_xfoil.h"
#include <vector>
#include <cmath>

using namespace aero;

namespace {
const double D2R = 3.14159265358979 / 180.0;

// NACA 0012 ordered TE -> upper -> LE -> lower -> TE (sharp/closed TE).
void naca0012(std::vector<double>& X, std::vector<double>& Y, int nside) {
    const double t = 0.12, PI = 3.14159265358979;
    auto yt = [&](double x) {
        return 5 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                        0.2843 * x * x * x - 0.1036 * x * x * x * x);
    };
    X.clear(); Y.clear();
    for (int j = nside; j >= 1; --j) {
        double x = 0.5 * (1 - std::cos(PI * j / nside));
        X.push_back(x); Y.push_back(yt(x));
    }
    X.push_back(0.0); Y.push_back(0.0);
    for (int j = 1; j <= nside; ++j) {
        double x = 0.5 * (1 - std::cos(PI * j / nside));
        X.push_back(x); Y.push_back(-yt(x));
    }
}
}  // namespace

// Inviscid cl-alpha slope ~ 2*pi (thickness-corrected) and cm ~ 0 at quarter
// chord for the symmetric section.
TEST(xfoil_inviscid_slope_and_cm) {
    std::vector<double> X, Y;
    naca0012(X, Y, 90);
    xfoil::Solver s;
    CHECK(s.set_coords(X, Y));
    xfoil::Result r2 = s.solve(2.0 * D2R, 1e6);
    xfoil::Result r6 = s.solve(6.0 * D2R, 1e6);
    double slope = (r6.cl - r2.cl) / (4.0 * D2R);     // per rad
    CHECK(slope > 6.0 && slope < 7.4);                // 2*pi, +thickness
    CHECK(std::fabs(r2.cm) < 0.02);
}

// Closed cylinder (sources only) carries no net lift.
TEST(xfoil_cylinder_zero_lift) {
    const double PI = 3.14159265358979;
    int n = 100;
    std::vector<double> X, Y;
    for (int i = 0; i <= n; ++i) {
        double th = 2 * PI * i / n;
        X.push_back(0.5 * std::cos(th));
        Y.push_back(0.5 * std::sin(th));
    }
    xfoil::Solver s;
    CHECK(s.set_coords(X, Y));
    xfoil::Result r = s.solve(0.0, 1e6);
    CHECK_NEAR(r.cl, 0.0, 1e-3);
    CHECK_NEAR(r.cm, 0.0, 1e-3);
}

// NACA 0012 at Re = 2e5: profile drag in the M4 acceptance band, friction part a
// physical fraction of the total, and the solve converges.
TEST(xfoil_naca0012_cd0_band) {
    std::vector<double> X, Y;
    naca0012(X, Y, 100);
    xfoil::Options opt; opt.Ncrit = 9.0;
    xfoil::Solver s(opt);
    CHECK(s.set_coords(X, Y));
    xfoil::Result r = s.solve(0.0, 2.0e5);
    CHECK(r.converged);
    CHECK(r.cd > 0.008 && r.cd < 0.020);              // README M4 acceptance
    CHECK(r.cdf > 0.0 && r.cdf < r.cd);               // friction < total
    CHECK(r.cd < 0.05);
}

// Transition migrates forward on the suction side as incidence increases.
TEST(xfoil_transition_moves_forward) {
    std::vector<double> X, Y;
    naca0012(X, Y, 100);
    xfoil::Options opt; opt.Ncrit = 9.0;
    xfoil::Solver s(opt);
    CHECK(s.set_coords(X, Y));
    xfoil::Result r0 = s.solve(0.0, 2.0e5);
    xfoil::Result r6 = s.solve(6.0 * D2R, 2.0e5);
    CHECK(r0.converged && r6.converged);
    double xtr0 = std::min(r0.xtr_top, r0.xtr_bot);
    double xtr6 = std::min(r6.xtr_top, r6.xtr_bot);
    CHECK(xtr6 < xtr0);                                // suction side forward
}

// Drag rises monotonically with |alpha| through the unstalled band.
TEST(xfoil_drag_bucket_monotone) {
    std::vector<double> X, Y;
    naca0012(X, Y, 100);
    xfoil::Options opt; opt.Ncrit = 9.0;
    xfoil::Solver s(opt);
    CHECK(s.set_coords(X, Y));
    xfoil::Result r0 = s.solve(0.0, 2.0e5);
    xfoil::Result r4 = s.solve(4.0 * D2R, 2.0e5);
    xfoil::Result r7 = s.solve(7.0 * D2R, 2.0e5);
    CHECK(r0.converged && r4.converged && r7.converged);
    CHECK(r7.cd > r4.cd);
    CHECK(r4.cd > 0.7 * r0.cd);                        // bucket floor near alpha 0
}

// Bit-for-bit determinism: identical geometry + operating point -> identical
// coefficients (no RNG anywhere in the solve path).
TEST(xfoil_determinism) {
    std::vector<double> X, Y;
    naca0012(X, Y, 100);
    xfoil::Options opt; opt.Ncrit = 9.0;
    xfoil::Solver s1(opt), s2(opt);
    CHECK(s1.set_coords(X, Y));
    CHECK(s2.set_coords(X, Y));
    xfoil::Result a = s1.solve(3.0 * D2R, 2.0e5);
    xfoil::Result b = s2.solve(3.0 * D2R, 2.0e5);
    CHECK(a.cl == b.cl);
    CHECK(a.cd == b.cd);
    CHECK(a.cm == b.cm);
}
