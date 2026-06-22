// mesh_probe.cpp — sanity-check the closed-body surface mesh.
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

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
    for (double AR : {6.0}) {
        WingGeometry w = rect_wing(0.20, AR, 20);
        double planform = 0.20 * (AR * 0.20);     // c * b
        panel::MeshStats s = panel::mesh_debug(w, cfg);
        std::printf("AR=%.1f  panels=%d strips=%d nc=%d  wetted=%.5f  planform=%.5f  wetted/2planform=%.3f  self_doublet=%.4f\n",
                    AR, s.n_panels, s.n_strips, s.nc, s.wetted_area, planform,
                    s.wetted_area / (2.0 * planform), s.self_doublet);
    }
    return 0;
}
