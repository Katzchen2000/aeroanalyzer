#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
using namespace aero;
static WingGeometry mk(double wash){ WingGeometry w; w.root_chord=0.25;w.tip_chord=0.25;w.semi_span=0.6;
  w.le_sweep=0;w.washout=wash*DEG2RAD; w.section.wu={0.12,0.12,0.12,0.12};w.section.wl={-0.12,-0.12,-0.12,-0.12};
  w.section.te_thick=0.004;w.battery_x=0.06; geom::loft(w,20); return w; }
int main(){ Config cfg;
  for(double wash : {0.0,-0.5,-3.0}){
    WingGeometry w=mk(wash);
    panel::MeshStats s=panel::mesh_debug(w,cfg);
    printf("wash=%.1f  panels=%d  inward_normals=%.0f\n", wash, s.n_panels, s.min_outward);
  }
  return 0; }
