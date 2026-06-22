#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
using namespace aero;
int main(){
  Config cfg; cfg.set("panel_wake_chords","30");
  printf("Plain rectangular AR=6, no twist, symmetric section.  Prandtl CLa=4.71\n");
  for(int nc : {4,6,8,10,12,16,20,24}){
    WingGeometry w; w.root_chord=0.2; w.tip_chord=0.2; w.semi_span=0.6;
    w.le_sweep=0; w.washout=0; w.section.wu={0.12,0.12,0.12,0.12};
    w.section.wl={-0.12,-0.12,-0.12,-0.12}; w.section.te_thick=0.002; w.battery_x=0.05;
    geom::loft(w,30);
    char b[8]; snprintf(b,8,"%d",nc); cfg.set("panel_chordwise",b);
    auto s1=panel::panel_solve_debug(w,cfg,1*DEG2RAD);
    auto s2=panel::panel_solve_debug(w,cfg,5*DEG2RAD);
    printf("nc=%2d  CLa=%9.4f  flux=%.2e\n", nc, (s2.cl-s1.cl)/(4*DEG2RAD), s2.sigma_flux);
  }
  return 0;
}
