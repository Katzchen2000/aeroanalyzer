#include "test_harness.h"
#include "aeroanalyzer/aero_neuralfoil.h"
#include <array>
#include <cmath>

using namespace aero;

// The weights live in the repo at data/Neurafoilbin; tests run from the project
// root. If they are absent (clean checkout without the .bin assets), skip rather
// than fail the suite.
static bool load_nf(nf::NeuralFoil& m) {
    if (m.load("data/Neurafoilbin")) return true;
    std::cerr << "[skip] NeuralFoil weights not found at data/Neurafoilbin\n";
    return false;
}

// Physical sanity: lift grows ~linearly with alpha, drag is a positive bucket,
// and confidence is a probability. Catches gross port errors (transposed
// weights, wrong activation, mis-decoded outputs) without a Python oracle.
TEST(neuralfoil_physics) {
    nf::NeuralFoil m;
    if (!load_nf(m)) return;

    std::array<double, 8> up = {0.18, 0.20, 0.18, 0.16, 0.15, 0.16, 0.14, 0.12};
    std::array<double, 8> lo = {-0.10, -0.12, -0.10, -0.08, -0.05, -0.02, 0.02, 0.05};
    const double te = 0.004, Re = 3.0e5, ncrit = 9.0, xtr = 1.0;

    double cl_prev = -1e9, cd_min = 1e9;
    for (double a = -4; a <= 10.0001; a += 2.0) {
        nf::Aero r = m.eval(up, lo, 0.0, te, a, Re, ncrit, xtr, xtr);
        CHECK(std::isfinite(r.cl) && std::isfinite(r.cd) && std::isfinite(r.cm));
        CHECK(r.cd > 0.0 && r.cd < 0.1);
        CHECK(r.confidence > 0.0 && r.confidence <= 1.0);
        CHECK(r.cl > cl_prev);            // monotonic lift over the linear range
        cl_prev = r.cl;
        cd_min = std::min(cd_min, r.cd);
    }
    // lift-curve slope ~ 2*pi/rad over -4..10 deg
    nf::Aero lo_a = m.eval(up, lo, 0.0, te, -4.0, Re, ncrit, xtr, xtr);
    nf::Aero hi_a = m.eval(up, lo, 0.0, te, 10.0, Re, ncrit, xtr, xtr);
    double slope_per_deg = (hi_a.cl - lo_a.cl) / 14.0;
    CHECK(slope_per_deg > 0.08 && slope_per_deg < 0.13);
    CHECK(cd_min > 0.004 && cd_min < 0.02);  // realistic low-Re drag bucket
}

// Up/down symmetry: a symmetric section has zero lift & moment at alpha=0, and
// CL(+a) = -CL(-a). This directly exercises the mirrored-pass averaging.
TEST(neuralfoil_symmetry) {
    nf::NeuralFoil m;
    if (!load_nf(m)) return;

    std::array<double, 8> up = {0.20, 0.20, 0.18, 0.16, 0.15, 0.14, 0.12, 0.10};
    std::array<double, 8> lo;
    for (int i = 0; i < 8; ++i) lo[i] = -up[i];
    const double te = 0.003, Re = 3.0e5, ncrit = 9.0, xtr = 1.0;

    nf::Aero r0 = m.eval(up, lo, 0.0, te, 0.0, Re, ncrit, xtr, xtr);
    CHECK_NEAR(r0.cl, 0.0, 0.02);
    CHECK_NEAR(r0.cm, 0.0, 0.02);

    nf::Aero rp = m.eval(up, lo, 0.0, te, 5.0, Re, ncrit, xtr, xtr);
    nf::Aero rn = m.eval(up, lo, 0.0, te, -5.0, Re, ncrit, xtr, xtr);
    CHECK_NEAR(rp.cl, -rn.cl, 1e-3);
    CHECK_NEAR(rp.cd, rn.cd, 1e-3);
}

// Golden cross-check against the reference NeuralFoil (Python, model_size="large").
// Inputs match tools/nf_golden.py; expected values are pasted from that script's
// output. Tolerance is 1% of |value| (or a small floor) to absorb float32 vs
// double accumulation. This is the gate that pins the transcribed I/O glue.
TEST(neuralfoil_golden) {
    nf::NeuralFoil m;
    if (!load_nf(m)) return;

    std::array<double, 8> up = {0.18, 0.20, 0.18, 0.16, 0.15, 0.16, 0.14, 0.12};
    std::array<double, 8> lo = {-0.10, -0.12, -0.10, -0.08, -0.05, -0.02, 0.02, 0.05};
    const double te = 0.004, Re = 3.0e5, ncrit = 9.0, xtr = 1.0;

    struct G { double alpha, cl, cd, cm; };
    // Reference values from tools/nf_golden.py (NeuralFoil model_size="large").
    // The C++ port matches these to ~5 digits once the two symmetry passes are
    // fused in raw output space before the nonlinear decode.
    const G ref[] = {
        {-4.0, -0.21794, 0.010114, -0.05356},
        { 0.0,  0.27227, 0.008529, -0.05775},
        { 4.0,  0.70667, 0.010340, -0.05329},
        { 8.0,  1.13872, 0.017293, -0.05297},
        {10.0,  1.30686, 0.024661, -0.04752},
    };
    for (const G& g : ref) {
        nf::Aero r = m.eval(up, lo, 0.0, te, g.alpha, Re, ncrit, xtr, xtr);
        CHECK_NEAR(r.cl, g.cl, std::max(0.01, 0.01 * std::fabs(g.cl)));
        CHECK_NEAR(r.cd, g.cd, std::max(2e-4, 0.01 * std::fabs(g.cd)));
        CHECK_NEAR(r.cm, g.cm, std::max(0.005, 0.01 * std::fabs(g.cm)));
    }
}
