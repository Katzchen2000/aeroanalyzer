#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/airfoil_io.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
using namespace aero;
int main(){
  bool ok=false; auto c=airfoil_io::load_dat("out/knee.dat",ok);
  double te=airfoil_io::estimate_te(c); Airfoil af=airfoil_io::to_airfoil(c,3,te);
  Config cfg; viscous::Surrogate s; s.load("data/surrogates/polar_coeffs.csv",cfg);
  printf("AVL: CLa=4.4557 Xnp=0.07306 e=0.4148  (Nchord=12,Nspan=24)\n");
  for(int nst : {20,30,40}) for(int nc : {8,14,20}){
    WingGeometry w; w.root_chord=0.222634; w.tip_chord=0.102248; w.semi_span=1.04751/2.0;
    w.le_sweep=8.00683*DEG2RAD; w.washout=-4.74184*DEG2RAD; w.section=af; geom::loft(w,nst);
    MassProps mp=massprops::compute(w,cfg);
    cfg.set("aero_model","panel"); char b[8]; snprintf(b,8,"%d",nc); cfg.set("panel_chordwise",b);
    cfg.set("panel_wake_chords","30");
    double a1=1*DEG2RAD,a2=5*DEG2RAD;
    auto s1=potential::solve(w,mp,s,cfg,a1,0.0); auto s2=potential::solve(w,mp,s,cfg,a2,0.0);
    auto sd=potential::solve(w,mp,s,cfg,2*DEG2RAD,0.0);
    printf("nst=%2d nc=%2d N=%4d  CLa=%.4f  x_np=%.5f  e=%.4f\n",
           nst,nc,2*(nst-1)*nc, (s2.CL-s1.CL)/(a2-a1), sd.x_np, sd.e);
  }
  return 0;
}
