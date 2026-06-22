// uniform_vs_grad.cpp - isolate whether the collapse is driven by twist GRADIENT
// (warped panels) or by twist magnitude per se. Rectangular unswept wing; we
// overwrite station twists directly after lofting.
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

using namespace aero;

static WingGeometry base() {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.25; w.semi_span = 0.6;
    w.le_sweep = 0; w.washout = 0;
    w.section.wu = {0.12,0.12,0.12,0.12};
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.004; w.battery_x = 0.06;
    geom::loft(w, 20);
    return w;
}

static double slope(WingGeometry& w, Config& cfg) {
    auto s0 = panel::panel_solve_debug(w, cfg, 0.0);
    auto s4 = panel::panel_solve_debug(w, cfg, 4.0*DEG2RAD);
    return (s4.cl - s0.cl)/(4.0*DEG2RAD);
}

int main() {
    Config cfg; cfg.set("panel_chordwise","8"); cfg.set("panel_wake_chords","20");

    // 1) uniform twist (constant on every station) - should behave like alpha shift
    for (double T : {0.0, 2.0, 5.0}) {
        WingGeometry w = base();
        for (auto& s : w.stations) s.twist = T*DEG2RAD;
        printf("uniform twist=%+.1f deg : slope=%.4f\n", T, slope(w, cfg));
    }
    // 2) linear gradient of varying magnitude (tip twist), root=0
    for (double tip : {-0.5, -1.0, -3.0}) {
        WingGeometry w = base();
        for (auto& s : w.stations) {
            double t = (w.semi_span>0)? s.y/w.semi_span : 0.0;
            s.twist = tip*DEG2RAD * t;
        }
        printf("gradient tip=%+.1f deg : slope=%.4f\n", tip, slope(w, cfg));
    }
    // 3) tiny single-station perturbation (one interior station nudged)
    {
        WingGeometry w = base();
        w.stations[10].twist = 0.5*DEG2RAD;   // one station only
        printf("single-station nudge   : slope=%.4f\n", slope(w, cfg));
    }
    return 0;
}
