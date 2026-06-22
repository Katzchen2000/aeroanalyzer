// twist_cmp.cpp - compare panel vs VLM spanwise loading for a twisted wing,
// to localize the washout collapse (shape vs response).
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

using namespace aero;

static WingGeometry mk(double wash) {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.25; w.semi_span = 0.6;
    w.le_sweep = 0; w.washout = wash * DEG2RAD;
    w.section.wu = {0.12,0.12,0.12,0.12};
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.004; w.battery_x = 0.06;
    geom::loft(w, 20);
    return w;
}

int main() {
    Config cfg; cfg.set("panel_chordwise","8"); cfg.set("panel_wake_chords","20");
    viscous::Surrogate surr; surr.load("", cfg);

    for (double wash : {0.0, -3.0}) {
        WingGeometry w = mk(wash);
        printf("\n=== washout = %.1f deg ===\n", wash);

        // Panel loading at two alphas
        auto p0 = panel::panel_loading_debug(w, cfg, 0.0);
        auto p4 = panel::panel_loading_debug(w, cfg, 4.0*DEG2RAD);

        printf(" PANEL gamma(y):\n");
        printf("  %-8s %12s %12s %12s\n","y","g@0","g@4","dg");
        for (size_t i=0;i<p0.y.size();++i)
            printf("  %-8.4f %12.6f %12.6f %12.6f\n",
                   p0.y[i], p0.gamma[i], p4.gamma[i], p4.gamma[i]-p0.gamma[i]);

        // VLM cl_local at two alphas
        cfg.set("aero_model","vlm");
        MassProps mp = massprops::compute(w, cfg);
        AeroState v0 = potential::solve(w, mp, surr, cfg, 0.0, 0.0);
        AeroState v4 = potential::solve(w, mp, surr, cfg, 4.0*DEG2RAD, 0.0);
        cfg.set("aero_model","vlm");
        printf(" VLM cl_local(station):\n");
        printf("  %-8s %12s %12s %12s\n","y","cl@0","cl@4","dcl");
        for (size_t i=0;i<w.stations.size();++i)
            printf("  %-8.4f %12.6f %12.6f %12.6f\n",
                   w.stations[i].y, v0.cl_local[i], v4.cl_local[i],
                   v4.cl_local[i]-v0.cl_local[i]);
        printf(" VLM CL @0=%.4f @4=%.4f slope=%.3f\n",
               v0.CL, v4.CL, (v4.CL-v0.CL)/(4.0*DEG2RAD));
    }
    return 0;
}
