// rect_cond.cpp - diagnose the panel chordwise non-convergence. For the plain
// rectangular AR=6 wing, sweep nc and report CLa alongside cond(C) and the
// smallest collocation-point separation. Run with PANEL_DEBUG_COND set so the
// conditioning fields are populated.
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
using namespace aero;

static WingGeometry rect() {
    WingGeometry w; w.root_chord = 0.2; w.tip_chord = 0.2; w.semi_span = 0.6;
    w.le_sweep = 0; w.washout = 0;
    w.section.wu = {0.12, 0.12, 0.12, 0.12};
    w.section.wl = {-0.12, -0.12, -0.12, -0.12};
    w.section.te_thick = 0.002; w.battery_x = 0.05;
    geom::loft(w, 30);
    return w;
}

static void sweep(const char* spacing) {
    Config cfg; cfg.set("panel_wake_chords", "30");
    cfg.set("panel_chord_spacing", spacing);
    printf("=== panel_chord_spacing = %s ===\n", spacing);
    printf("%4s  %10s  %12s  %12s\n", "nc", "CLa", "cond(C)", "min_cp_gap");
    for (int nc : {4, 6, 8, 10, 12, 16, 20, 24}) {
        WingGeometry w = rect();
        char b[8]; snprintf(b, 8, "%d", nc); cfg.set("panel_chordwise", b);
        auto s1 = panel::panel_solve_debug(w, cfg, 1 * DEG2RAD);
        auto s2 = panel::panel_solve_debug(w, cfg, 5 * DEG2RAD);
        double cla = (s2.cl - s1.cl) / (4 * DEG2RAD);
        printf("%4d  %10.4f  %12.3e  %12.3e\n", nc, cla, s2.cond_C, s2.min_cp_gap);
    }
}

int main() {
    setvbuf(stdout, nullptr, _IONBF, 0);   // stream each row as it computes
    printf("Plain rectangular AR=6, no twist, symmetric section.  Prandtl CLa=4.71\n");
    sweep("cosine");
    sweep("halfcosine");
    return 0;
}
