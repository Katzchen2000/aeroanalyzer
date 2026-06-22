// panel_solve_probe.cpp — validate the Morino solve: Trefftz lift slope vs AR
// (must rise toward 2*pi and track Prandtl) and the closed-body source flux.
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
#include <cmath>

using namespace aero;

static WingGeometry rect_wing(double chord, double AR, int nst) {
    WingGeometry w;
    double b = AR * chord;
    w.root_chord = chord; w.tip_chord = chord; w.semi_span = 0.5 * b;
    w.le_sweep = 0; w.washout = 0;
    w.section.wu = { 0.12, 0.12, 0.12, 0.12};
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.002; w.battery_x = 0.05;
    geom::loft(w, nst);
    return w;
}

int main() {
    Config cfg;
    const double PI = 3.14159265358979323846;
    std::printf("--- Morino panel solve: lift slope vs AR ---\n");
    for (double AR : {4.0, 6.0, 8.0, 10.0}) {
        WingGeometry* w = new WingGeometry(rect_wing(0.20, AR, 20));  // unique addr
        double a1 = 1.0 * PI / 180.0, a2 = 5.0 * PI / 180.0;
        panel::PanelSolveStats s1 = panel::panel_solve_debug(*w, cfg, a1);
        panel::PanelSolveStats s2 = panel::panel_solve_debug(*w, cfg, a2);
        double slope = (s2.cl - s1.cl) / (a2 - a1);
        double prandtl = 2.0 * PI * AR / (AR + 2.0);
        std::printf("AR=%5.1f N=%d  CL@1=%+.4f CL@5=%+.4f  slope=%.4f  prandtl=%.4f  err=%+.1f%%  flux=%.2e gmax=%.4f\n",
                    AR, s1.n_panels, s1.cl, s2.cl, slope, prandtl,
                    100.0 * (slope - prandtl) / prandtl, s2.sigma_flux, s2.gamma_max);
    }
    return 0;
}
