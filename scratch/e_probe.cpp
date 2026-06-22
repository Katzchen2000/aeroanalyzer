// e_probe.cpp - diagnose panel span-efficiency (Glauert fit) collapse.
// Reproduces the exact 5-mode odd-series fit from span_efficiency() on the
// panel's spanwise loading, instrumenting the C_n so we can see WHY e pins to
// the 0.30 floor. Compares a clean rectangular case (expect e high, ->1) to the
// old-knee planform (AVL oracle e=0.4148 @ a=2deg) and a taper/sweep/washout
// sweep.
//
// Build: g++ -std=c++17 -O2 -fno-math-errno -fopenmp -DEIGEN_DONT_PARALLELIZE
//        -I include -I <eigen> src/*.cpp scratch/e_probe.cpp -o build/e_probe.exe
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/stability.h"
#include "aeroanalyzer/config.h"
#include <Eigen/Dense>
#include <cstdio>
#include <vector>
#include <cmath>

using namespace aero;

// Replicate span_efficiency()'s fit on a HALF-wing loading (y>=0, gamma per
// half strip), printing the C_n and the raw (pre-clamp) e.
static double fit_e(const std::vector<double>& yhalf,
                    const std::vector<double>& ghalf, double semi, bool verbose) {
    int K = (int)yhalf.size();
    if (K < 1 || semi <= 0) return 1.0;
    std::vector<double> Y, G;
    for (int k = K - 1; k >= 0; --k) { Y.push_back(-yhalf[k]); G.push_back(ghalf[k]); }
    for (int k = 0; k < K; ++k)      { Y.push_back( yhalf[k]); G.push_back(ghalf[k]); }
    const int modes[5] = {1, 3, 5, 7, 9};
    const int M = 5;
    int KK = (int)Y.size();
    Eigen::MatrixXd A(KK, M);
    Eigen::VectorXd b(KK);
    for (int k = 0; k < KK; ++k) {
        double yy = std::max(-0.999999, std::min(0.999999, Y[k] / semi));
        double th = std::acos(-yy);
        for (int j = 0; j < M; ++j) A(k, j) = std::sin(modes[j] * th);
        b(k) = G[k];
    }
    Eigen::VectorXd C = A.colPivHouseholderQr().solve(b);
    double num = C(0) * C(0), den = 0.0;
    for (int j = 0; j < M; ++j) den += modes[j] * C(j) * C(j);
    double e = (den > 1e-30) ? num / den : 1.0;
    if (verbose) {
        printf("    C_n = [");
        for (int j = 0; j < M; ++j) printf(" %+.4f", C(j));
        printf(" ]  raw_e=%.4f\n", e);
    }
    return e;
}

static void probe(const char* tag, WingGeometry w, double avl_e) {
    geom::loft(w, 20);
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "10"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate s; s.load("data/surrogates/polar_coeffs.csv", cfg);
    MassProps mp = massprops::compute(w, cfg);

    double a = 2.0 * DEG2RAD;
    panel::LoadingDump ld = panel::panel_loading_debug(w, cfg, a);
    AeroState st = potential::solve(w, mp, s, cfg, a, 0.0);
    AeroState tr = stability::trim(w, mp, s, cfg);   // operating point (drag obj)

    // semi-span estimate matching span_efficiency: y_back + 0.5*dy_back.
    double semi = w.semi_span;
    printf("[%s] sweep=%.1f taper=%.2f wash=%.2f  a2: CL=%.4f e=%.4f | "
           "TRIM: alpha=%.2f CL=%.4f e=%.4f CDi=%.5f trimmed=%d",
           tag, w.le_sweep * RAD2DEG, w.tip_chord / w.root_chord,
           w.washout * RAD2DEG, st.CL, st.e,
           tr.alpha * RAD2DEG, tr.CL, tr.e, tr.CDi, (int)tr.trimmed);
    if (avl_e > 0) printf("  (AVL e=%.4f @a2)", avl_e);
    printf("\n");
    // dump loading shape
    printf("    y/s, gamma:");
    for (size_t i = 0; i < ld.y.size(); ++i)
        printf(" (%.2f,%.4f)", ld.y[i] / semi, ld.gamma[i]);
    printf("\n");
    fit_e(ld.y, ld.gamma, semi, true);
}

static WingGeometry rect(double chord, double AR) {
    WingGeometry w;
    w.root_chord = chord; w.tip_chord = chord;
    w.semi_span = 0.5 * AR * chord;
    w.le_sweep = 0; w.washout = 0;
    w.section.wu = { 0.12, 0.12, 0.12, 0.12};
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.002; w.battery_x = 0.05;
    return w;
}

int main() {
    // 1. clean rectangular unswept -> e should be high (near 1).
    probe("rect_AR6", rect(0.20, 6.0), -1);

    // 2. old-knee planform (matches the baked AVL oracle e=0.4148 @ a=2deg).
    {
        WingGeometry w;
        w.root_chord = 0.222634; w.tip_chord = 0.102248;
        w.semi_span = 1.04751 / 2.0;
        w.le_sweep = 8.00683 * DEG2RAD; w.washout = -4.74184 * DEG2RAD;
        w.section.wu = {0.20, 0.17, 0.14, 0.11};
        w.section.wl = {-0.12, -0.09, -0.02, 0.06};
        w.section.te_thick = 0.004; w.battery_x = 0.06;
        probe("old_knee", w, 0.4148);
    }

    // 3. taper/sweep/washout sensitivity sweep.
    for (double tp : {0.4, 0.7, 1.0}) {
        WingGeometry w = rect(0.22, 6.0);
        w.tip_chord = tp * w.root_chord;
        probe("taper", w, -1);
    }
    for (double sw : {0.0, 15.0, 30.0}) {
        WingGeometry w = rect(0.22, 6.0);
        w.le_sweep = sw * DEG2RAD;
        probe("sweep", w, -1);
    }
    for (double wo : {0.0, -3.0, -6.0}) {
        WingGeometry w = rect(0.22, 6.0);
        w.washout = wo * DEG2RAD;
        probe("washout", w, -1);
    }
    return 0;
}
