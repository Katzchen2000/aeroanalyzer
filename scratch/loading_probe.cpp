// loading_probe.cpp — dump spanwise loading for twisted vs untwisted wings.
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

using namespace aero;

static WingGeometry mk(double washdeg) {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.25; w.semi_span = 0.6;
    w.le_sweep = 0; w.washout = washdeg * DEG2RAD;
    w.section.wu = {0.12,0.12,0.12,0.12};
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.004; w.battery_x = 0.06;
    geom::loft(w, 20);
    return w;
}

static void dump(const char* nm, double wash, double adeg) {
    Config cfg;
    WingGeometry w = mk(wash);
    panel::LoadingDump d = panel::panel_loading_debug(w, cfg, adeg * DEG2RAD);
    std::printf("%s (washout=%.0f, alpha=%.0f): ", nm, wash, adeg);
    for (size_t i = 0; i < d.y.size(); i += 4)
        std::printf("%+.4f ", d.gamma[i]);
    std::printf("\n");
}

int main() {
    dump("untwisted", 0, 0);
    dump("untwisted", 0, 8);
    dump("washout  ", -3, 0);
    dump("washout  ", -3, 8);
    return 0;
}
