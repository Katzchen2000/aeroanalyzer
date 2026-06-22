// vlm_probe.cpp — independent check of the VLM lift slope vs Prandtl.
// Builds a rectangular, unswept, symmetric-section wing and compares the
// panel-solver CL/alpha against 2*pi*AR/(AR+2).
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
#include <cmath>

using namespace aero;

static WingGeometry rect_wing(double chord, double AR, int nst) {
    WingGeometry w;
    double b_full = AR * chord;          // rectangular: AR = b/c
    w.root_chord = chord; w.tip_chord = chord;
    w.semi_span  = 0.5 * b_full;
    w.le_sweep = 0.0; w.washout = 0.0;
    // Thin symmetric section: zero camber -> alpha_L0 = 0, cm_ac = 0.
    w.section.wu = { 0.12, 0.12, 0.12, 0.12};
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.002;
    w.battery_x = 0.05;
    geom::loft(w, nst);
    return w;
}

static void probe(double AR) {
    Config cfg;
    // Heap-allocate so each wing has a UNIQUE address (defeats the &w cache key).
    WingGeometry* wp = new WingGeometry(rect_wing(0.20, AR, 40));
    WingGeometry& w = *wp;
    MassProps mp = massprops::compute(w, cfg);
    viscous::Surrogate surr; surr.load("", cfg);

    double a1 = 1.0 * DEG2RAD, a2 = 5.0 * DEG2RAD;
    AeroState s1 = potential::solve(w, mp, surr, cfg, a1, 0.0);
    AeroState s2 = potential::solve(w, mp, surr, cfg, a2, 0.0);
    double slope = (s2.CL - s1.CL) / (a2 - a1);
    double prandtl = 2.0 * PI * mp.AR / (mp.AR + 2.0);
    std::printf("AR=%.2f  CL@1deg=%.4f CL@5deg=%.4f  slope=%.4f  prandtl=%.4f  err=%+.1f%%  e=%.3f x_np=%.4f mac=%.4f\n",
                mp.AR, s1.CL, s2.CL, slope, prandtl,
                100.0 * (slope - prandtl) / prandtl, s1.e, s1.x_np, mp.mac);
}

int main() {
    std::printf("--- VLM lift-slope probe (rectangular, unswept, symmetric) ---\n");
    probe(4.0);
    probe(6.0);
    probe(8.0);
    probe(10.0);
    return 0;
}
