// panel_vs_vlm.cpp - final cross-check: panel vs validated VLM across
// representative wings (incl. sweep + washout + reflex).
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

using namespace aero;

static double cl(WingGeometry& w, MassProps& mp, viscous::Surrogate& s, Config cfg,
                const char* model, double a) {
    cfg.set("aero_model", model);
    cfg.set("panel_chordwise","8"); cfg.set("panel_wake_chords","20");
    return potential::solve(w, mp, s, cfg, a, 0.0).CL;
}

static void run(const char* lbl, WingGeometry w) {
    Config cfg; viscous::Surrogate s; s.load("", cfg);
    MassProps mp = massprops::compute(w, cfg);
    double a1=1.0*DEG2RAD, a2=5.0*DEG2RAD;
    double vp1=cl(w,mp,s,cfg,"panel",a1), vp2=cl(w,mp,s,cfg,"panel",a2);
    double vv1=cl(w,mp,s,cfg,"vlm",a1),   vv2=cl(w,mp,s,cfg,"vlm",a2);
    double prandtl = 2.0*PI*mp.AR/(mp.AR+2.0);
    printf("%-28s AR=%.2f  panel: CL@1=%+.4f slope=%.3f | vlm: CL@1=%+.4f slope=%.3f | prandtl=%.3f\n",
           lbl, mp.AR, vp1,(vp2-vp1)/(a2-a1), vv1,(vv2-vv1)/(a2-a1), prandtl);
}

int main() {
    { WingGeometry w; w.root_chord=0.20;w.tip_chord=0.20;w.semi_span=0.5;
      w.le_sweep=0;w.washout=0; w.section.wu={0.12,0.12,0.12,0.12};
      w.section.wl={-0.12,-0.12,-0.12,-0.12};w.section.te_thick=0.002;w.battery_x=0.05;
      geom::loft(w,24); run("rect, no twist", w); }
    { WingGeometry w; w.root_chord=0.20;w.tip_chord=0.20;w.semi_span=0.5;
      w.le_sweep=0;w.washout=-3*DEG2RAD; w.section.wu={0.12,0.12,0.12,0.12};
      w.section.wl={-0.12,-0.12,-0.12,-0.12};w.section.te_thick=0.002;w.battery_x=0.05;
      geom::loft(w,24); run("rect, -3 washout", w); }
    { WingGeometry w; w.root_chord=0.25;w.tip_chord=0.13;w.semi_span=0.6;
      w.le_sweep=18*DEG2RAD;w.washout=-3*DEG2RAD; w.section.wu={0.20,0.17,0.14,0.11};
      w.section.wl={-0.12,-0.09,-0.02,0.06};w.section.te_thick=0.004;w.battery_x=0.06;
      geom::loft(w,24); run("swept+washout+reflex (demo)", w); }
    return 0;
}
