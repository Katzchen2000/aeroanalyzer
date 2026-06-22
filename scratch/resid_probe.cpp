#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
#include <iostream>
#include <exception>
#include <chrono>

using namespace aero;

static WingGeometry mk(double wash) {
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
    geom::loft(w, 20);
    return w;
}

int main() {
    // Use a small panel count for quick testing
    Config cfg;
    cfg.set("panel_chordwise", "2");   // small matrix (76 panels)
    cfg.set("panel_wake_chords", "2");
    
    printf("Testing with panel_chordwise=2 (small matrix)...\n");
    fflush(stdout);
    
    try {
        // Test only one case to see if it completes
        double wash = -3.0;
        printf("Building wing...");
        fflush(stdout);
        WingGeometry w = mk(wash);
        printf(" done.\n");
        fflush(stdout);
        
        printf("Calling panel_solve_debug at alpha=0...");
        fflush(stdout);
        auto start = std::chrono::high_resolution_clock::now();
        auto s0 = panel::panel_solve_debug(w, cfg, 0.0);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        printf(" done in %lld ms. CL@0=%+.4f\n", elapsed.count(), s0.cl);
        fflush(stdout);
        
        printf("Calling panel_solve_debug at alpha=4 deg...");
        fflush(stdout);
        start = std::chrono::high_resolution_clock::now();
        auto s4 = panel::panel_solve_debug(w, cfg, 4.0 * DEG2RAD);
        end = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        printf(" done in %lld ms. CL@4=%+.4f slope=%.3f\n", 
               elapsed.count(), s4.cl, (s4.cl - s0.cl) / (4.0 * DEG2RAD));
        fflush(stdout);
        
        // Now test the full wake-length scan with larger panels
        printf("\nNow testing with panel_chordwise=10 (default) for comparison...\n");
        fflush(stdout);
        
        for (double wl : {2.0, 10.0, 50.0}) {
            Config cfg2;
            cfg2.set("panel_chordwise", "10");
            char b[32]; snprintf(b, 32, "%g", wl);
            cfg2.set("panel_wake_chords", b);
            for (double wash2 : {0.0, -3.0}) {
                printf("  wake=%g wash=%.0f...", wl, wash2);
                fflush(stdout);
                WingGeometry w2 = mk(wash2);
                auto start2 = std::chrono::high_resolution_clock::now();
                auto s0_2 = panel::panel_solve_debug(w2, cfg2, 0.0);
                auto s4_2 = panel::panel_solve_debug(w2, cfg2, 4.0 * DEG2RAD);
                auto end2 = std::chrono::high_resolution_clock::now();
                auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
                printf(" done (%lld ms). CL@0=%+.4f CL@4=%+.4f slope=%.3f\n",
                       elapsed2.count(), s0_2.cl, s4_2.cl, (s4_2.cl - s0_2.cl) / (4.0 * DEG2RAD));
                fflush(stdout);
            }
        }
    } catch (const std::exception& e) {
        printf(" Exception: %s\n", e.what());
        fflush(stdout);
    } catch (...) {
        printf(" Unknown exception!\n");
        fflush(stdout);
    }
    return 0;
}