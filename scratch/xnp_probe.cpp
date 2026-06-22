// xnp_probe.cpp - compare the panel neutral point from the quarter-chord proxy
// vs the chordwise load-centre integration, against the AVL oracle.
//   old_knee AVL (a=2deg): Xnp=0.07306  (planform matches the baked oracle).
// A rect unswept wing must give x_np = 0.25c exactly for BOTH methods.
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

using namespace aero;

static void probe(const char* tag, WingGeometry w, double avl_xnp, int nc) {
    geom::loft(w, 20);
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", std::to_string(nc));
    cfg.set("panel_wake_chords", "20");
    viscous::Surrogate s; s.load("data/surrogates/polar_coeffs.csv", cfg);
    MassProps mp = massprops::compute(w, cfg);
    double a = 2.0 * DEG2RAD;
    panel::XnpDebug d = panel::panel_xnp_debug(w, cfg, a);
    double mac = mp.mac;
    printf("[%s nc=%d] mac=%.4f  proxy=%.5f load=%.5f  CL=%.4f", tag, nc, mac,
           d.proxy, d.load, d.cl);
    if (avl_xnp > 0) {
        printf("  AVL=%.5f | proxy-AVL=%+.1f%%MAC  load-AVL=%+.1f%%MAC",
               avl_xnp, 100.0 * (d.proxy - avl_xnp) / mac,
               100.0 * (d.load - avl_xnp) / mac);
    }
    printf("\n");
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

static WingGeometry old_knee() {
    WingGeometry w;
    w.root_chord = 0.222634; w.tip_chord = 0.102248;
    w.semi_span = 1.04751 / 2.0;
    w.le_sweep = 8.00683 * DEG2RAD; w.washout = -4.74184 * DEG2RAD;
    w.section.wu = {0.20, 0.17, 0.14, 0.11};
    w.section.wl = {-0.12, -0.09, -0.02, 0.06};
    w.section.te_thick = 0.004; w.battery_x = 0.06;
    return w;
}

int main() {
    // rect unswept: both methods must land on the quarter chord (0.25c = 0.05).
    for (int nc : {6, 10, 16}) probe("rect_AR6", rect(0.20, 6.0), 0.05, nc);
    // old-knee planform vs the baked AVL Xnp.
    for (int nc : {6, 10, 16, 24}) probe("old_knee", old_knee(), 0.07306, nc);
    // sweep sweep (no AVL number, just observe proxy vs load divergence).
    for (double sw : {0.0, 15.0, 30.0}) {
        WingGeometry w = rect(0.22, 6.0); w.le_sweep = sw * DEG2RAD;
        probe("sweep", w, -1, 10);
    }
    return 0;
}
