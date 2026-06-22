#include "../include/geom.h"
#include "../include/panel.h"
#include "../include/config.h"
#include <cstdio>
#include <vector>

using namespace aero;

static WingGeometry mk(double wash, int n) {
    WingGeometry w;
    w.root_chord = 0.25;
    w.tip_chord = 0.25;
    w.semi_span = 0.6;
    w.le_sweep = 0;
    w.washout = wash * DEG2RAD;
    w.section.wu = {0.12, 0.12, 0.12, 0.12};
    w.section.wl = {-0.12, -0.12, -0.12, -0.12};
    w.section.te_thick = 0.004;
    w.battery_x = 0.06;
    geom::loft(w, n);
    return w;
}

int main() {
    Config cfg;
    
    for (double wash : {0.0, -0.5, -3.0}) {
        printf("\n=== wash=%.1f deg ===\n", wash);
        WingGeometry w = mk(wash, 20);
        
        // Build mesh
        auto mesh = panel::build_mesh(w, 20);  // You may need to expose this
        
        // Check panel geometry
        double areas[4] = {0};
        double normals[4] = {0};
        double z_range = 0;
        double z_min = 1e30, z_max = -1e30;
        
        for (const auto& p : mesh.panels) {
            // Check area distribution
            int idx = (p.surf ? 2 : 0) + (p.jc < 10 ? 0 : 1);
            areas[idx] += p.area;
            
            // Check normal consistency
            normals[p.surf]++;
            
            // Check z-coordinates
            for (int k = 0; k < 4; ++k) {
                z_min = std::min(z_min, p.q.c[k].z);
                z_max = std::max(z_max, p.q.c[k].z);
            }
        }
        
        printf("  Areas: upper-LE=%.4f upper-TE=%.4f lower-LE=%.4f lower-TE=%.4f\n",
               areas[2], areas[3], areas[0], areas[1]);
        printf("  Normal counts: upper=%d lower=%d\n", 
               (int)normals[1], (int)normals[0]);
        printf("  Z range: [%.4f, %.4f] (span=%.4f)\n", z_min, z_max, z_max - z_min);
    }
    return 0;
}