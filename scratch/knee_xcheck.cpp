// knee_xcheck.cpp - reconstruct the knee design (pareto idx 7) from out/knee.dat
// + pareto params, run panel & VLM, compare to the AVL oracle.
//   AVL (knee.avl @ a=2deg): CLa=4.4557  e=0.4148  Xnp=0.07306  CG=0.14207
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/airfoil_io.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>

using namespace aero;

static void run(const char* model, WingGeometry& w, MassProps& mp,
                viscous::Surrogate& s, Config cfg) {
    cfg.set("aero_model", model);
    cfg.set("panel_chordwise","10"); cfg.set("panel_wake_chords","20");
    double a1=1.0*DEG2RAD, a2=5.0*DEG2RAD;
    AeroState s1 = potential::solve(w, mp, s, cfg, a1, 0.0);
    AeroState s2 = potential::solve(w, mp, s, cfg, a2, 0.0);
    AeroState s2d = potential::solve(w, mp, s, cfg, 2.0*DEG2RAD, 0.0);
    double CLa = (s2.CL - s1.CL)/(a2-a1);
    printf("  %-6s  CLa=%.4f  e=%.4f  x_np=%.5f  SM=%.4f  (CL@2=%.4f)\n",
           model, CLa, s2d.e, s2d.x_np, s2d.static_margin, s2d.CL);
}

int main() {
    bool ok=false;
    auto c = airfoil_io::load_dat("out/knee.dat", ok);
    if (!ok){ printf("failed to load knee.dat\n"); return 1; }
    double te = airfoil_io::estimate_te(c);
    Airfoil af = airfoil_io::to_airfoil(c, 3, te);   // 4 weights/surface

    WingGeometry w;
    w.root_chord = 0.222634;
    w.tip_chord  = 0.102248;
    w.semi_span  = 1.04751/2.0;
    w.le_sweep   = 8.00683*DEG2RAD;
    w.washout    = -4.74184*DEG2RAD;
    w.section    = af;
    geom::loft(w, 20);

    Config cfg; viscous::Surrogate s; s.load("data/surrogates/polar_coeffs.csv", cfg);
    MassProps mp = massprops::compute(w, cfg);
    printf("reconstructed: S_ref=%.5f mac=%.5f b=%.5f x_cg=%.5f AR=%.3f te=%.4f\n",
           mp.S_ref, mp.mac, mp.b_full, mp.x_cg, mp.AR, te);
    printf("AVL oracle  :  CLa=4.4557  e=0.4148  Xnp=0.07306  CG=0.14207\n");
    run("vlm",   w, mp, s, cfg);
    run("panel", w, mp, s, cfg);
    return 0;
}
