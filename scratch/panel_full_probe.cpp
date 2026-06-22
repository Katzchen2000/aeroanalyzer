// panel_full_probe.cpp — exercise the full panel AeroState + trim.
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/stability.h"
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

static WingGeometry demo_wing() {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.13; w.semi_span = 0.6;
    w.le_sweep = 18.0 * DEG2RAD; w.washout = -3.0 * DEG2RAD;
    w.section.wu = {0.20, 0.17, 0.14, 0.11};
    w.section.wl = {-0.12, -0.09, -0.02, 0.06};
    w.section.te_thick = 0.004; w.battery_x = 0.06;
    geom::loft(w, 20);
    return w;
}

int main() {
    Config cfg; cfg.set("aero_model", "panel");
    viscous::Surrogate surr; surr.load("", cfg);

    std::printf("--- panel AeroState: rectangular lift slope & span eff ---\n");
    for (double AR : {4.0, 6.0, 8.0, 10.0}) {
        WingGeometry* w = new WingGeometry(rect_wing(0.20, AR, 20));
        MassProps mp = massprops::compute(*w, cfg);
        AeroState s1 = potential::solve(*w, mp, surr, cfg, 1.0 * DEG2RAD, 0.0);
        AeroState s2 = potential::solve(*w, mp, surr, cfg, 5.0 * DEG2RAD, 0.0);
        double slope = (s2.CL - s1.CL) / (4.0 * DEG2RAD);
        double prandtl = 2.0 * PI * mp.AR / (mp.AR + 2.0);
        std::printf("AR=%5.1f slope=%.4f prandtl=%.4f err=%+.1f%%  e=%.3f CDi=%.5f CDp=%.5f Cm=%+.4f\n",
                    mp.AR, slope, prandtl, 100*(slope-prandtl)/prandtl, s2.e, s2.CDi, s2.CDp, s2.CM);
    }

    std::printf("\n--- isolate factors (CL@4deg, slope, e) ---\n");
    auto test = [&](const char* name, double sweepdeg, double taper, double washdeg,
                    bool cambered) {
        WingGeometry* w = new WingGeometry();
        w->root_chord = 0.25; w->tip_chord = 0.25 * taper; w->semi_span = 0.6;
        w->le_sweep = sweepdeg * DEG2RAD; w->washout = washdeg * DEG2RAD;
        if (cambered) { w->section.wu = {0.20,0.17,0.14,0.11}; w->section.wl = {-0.12,-0.09,-0.02,0.06}; }
        else          { w->section.wu = {0.12,0.12,0.12,0.12}; w->section.wl = {-0.12,-0.12,-0.12,-0.12}; }
        w->section.te_thick = 0.004; w->battery_x = 0.06;
        geom::loft(*w, 20);
        MassProps mp = massprops::compute(*w, cfg);
        double cl0 = potential::solve(*w, mp, surr, cfg, 0.0, 0.0).CL;
        double cl2 = potential::solve(*w, mp, surr, cfg, 2.0*DEG2RAD, 0.0).CL;
        double cl4 = potential::solve(*w, mp, surr, cfg, 4.0*DEG2RAD, 0.0).CL;
        double cl8 = potential::solve(*w, mp, surr, cfg, 8.0*DEG2RAD, 0.0).CL;
        double slope = (cl4 - cl2)/(2.0*DEG2RAD);
        std::printf("%-22s CL@0=%+.4f @2=%+.4f @4=%+.4f @8=%+.4f slope=%.3f AR=%.2f\n",
                    name, cl0, cl2, cl4, cl8, slope, mp.AR);
    };
    test("rect symmetric",      0,  1.0,  0, false);
    test("rect cambered",       0,  1.0,  0, true);
    test("tapered symmetric",   0,  0.52, 0, false);
    test("swept symmetric",    18,  1.0,  0, false);
    test("washout symmetric",   0,  1.0, -3, false);
    test("swept+taper sym",    18,  0.52, 0, false);
    test("demo (all)",         18,  0.52,-3, true);

    std::printf("\n--- demo wing (swept, reflex) trim ---\n");
    WingGeometry w = demo_wing();
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = stability::trim(w, mp, surr, cfg);
    double q = 0.5 * RHO * V_CRUISE * V_CRUISE;
    double CL_req = mp.mass * GRAV / (q * mp.S_ref);
    std::printf("trimmed=%d  CL=%.5f CL_req=%.5f  CM=%+.2e  alpha=%.3fdeg delta=%.3fdeg SM=%.3f e=%.3f\n",
                (int)st.trimmed, st.CL, CL_req, st.CM, st.alpha*RAD2DEG, st.delta_e*RAD2DEG,
                st.static_margin, st.e);
    return 0;
}
