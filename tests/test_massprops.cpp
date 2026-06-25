#include "test_harness.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"

using namespace aero;

static WingGeometry rect_wing() {
    WingGeometry w;
    w.root_chord = 0.20; w.tip_chord = 0.20; w.semi_span = 0.60;
    w.le_sweep = 0.0; w.washout = 0.0;
    Airfoil af; af.wu = {0.18, 0.15, 0.12, 0.10};
    af.wl = {-0.18, -0.15, -0.12, -0.10};
    af.te_thick = 0.0;
    w.sections.assign(1, af);
    geom::loft(w, 20);
    return w;
}

// With infill=1 the structural volume must equal the enclosed section volume:
// 2 * A_hat * chord^2 * semi_span for a constant-chord wing. Validates the
// spanwise trapezoidal integration against a closed form.
TEST(volume_solid_matches_closed_form) {
    WingGeometry w = rect_wing();
    Config cfg;
    cfg.set("infill_root", "1.0");
    cfg.set("infill_tip", "1.0");
    cfg.set("shell_thickness", "0.0001");
    cfg.set("material_density", "1000");
    MassProps mp = massprops::compute(w, cfg);
    double A_hat = massprops::section_area_hat(w.sections[0]);
    double expect = 2.0 * A_hat * 0.20 * 0.20 * 0.60;
    CHECK_NEAR(mp.volume, expect, 1e-6);
}

// Structure-only CG sits at the area centroid (~0.42c) for an unswept wing.
TEST(cg_structure_only) {
    WingGeometry w = rect_wing();
    Config cfg;
    cfg.set("mass_motor", "0");
    cfg.set("mass_battery", "0");
    cfg.set("target_mass_aux", "0");
    cfg.set("mass_servo_each", "0");
    MassProps mp = massprops::compute(w, cfg);
    CHECK_NEAR(mp.x_cg, 0.42 * 0.20, 1e-6);
    CHECK(mp.mass > 0.0);
}

// Battery shifting aft must move the CG aft (the SM trim handle, plan §4).
TEST(battery_shift_moves_cg) {
    WingGeometry w = rect_wing();
    Config cfg;
    w.battery_x = 0.02;
    double cg_fwd = massprops::compute(w, cfg).x_cg;
    w.battery_x = 0.18;
    double cg_aft = massprops::compute(w, cfg).x_cg;
    CHECK(cg_aft > cg_fwd);
}

// Spar clearance is finite and bounded by the local half-thickness.
TEST(spar_clearance_bounded) {
    WingGeometry w = rect_wing();
    Config cfg;
    MassProps mp = massprops::compute(w, cfg);
    double half_t_max = 0.5 * (geom::cst_upper(w.sections[0], 0.15) -
                               geom::cst_lower(w.sections[0], 0.15)) * w.root_chord;
    CHECK(mp.spar_clearance < half_t_max + 1e-9);
}
