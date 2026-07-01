#include "test_harness.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"

using namespace aero;

static WingGeometry rect_wing() {
    WingGeometry w;
    w.semi_span = 0.60;
    geom::set_linear_planform(w, 0.20, 0.20, 0.0, 0.0);
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

// Prop keep-out: heavily-swept wing whose near-axis TE reaches aft of the fixed
// disk face must produce a negative prop_clearance (gate fires).
TEST(prop_clearance_fires_on_aft_te) {
    WingGeometry w;
    // Very long root chord so TE = root_chord = 0.50 m >> prop face at 0.203/2+0.01=~0.22 m
    w.semi_span = 0.60;
    geom::set_linear_planform(w, 0.50, 0.10, 0.0, 0.0);
    Airfoil af; af.wu = {0.18,0.15,0.12,0.10}; af.wl = {-0.12,-0.10,-0.08,-0.05};
    af.te_thick = 0.0;
    w.sections.assign(1, af);
    geom::loft(w, 20);
    Config cfg;
    cfg.set("prop_diameter", "0.203");
    cfg.set("prop_hub_gap",  "0.010");
    cfg.set("prop_blade_clear", "0.005");
    MassProps mp = massprops::compute(w, cfg);
    // face = 0.50 + 0.010 - 0.005 = 0.505; root TE at 0.50 → cl = 0.505-0.50 = 0.005 > 0
    // Use a more extreme case: zero sweep, root_chord=0.50, prop face = 0.505 → barely clear
    // Make root_chord > face to breach: set root_chord=0.60 via re-loft
    WingGeometry w2;
    w2.semi_span = 0.60;
    geom::set_linear_planform(w2, 0.60, 0.10, 0.0, 0.0);
    w2.sections.assign(1, af);
    geom::loft(w2, 20);
    MassProps mp2 = massprops::compute(w2, cfg);
    // face = 0.60 + 0.010 - 0.005 = 0.605; root TE = 0.60 → cl = 0.005 (just clear)
    // Increase prop_hub_gap to zero → face = 0.60 - 0.005 = 0.595 < 0.60 → negative
    cfg.set("prop_hub_gap", "0.0");
    cfg.set("prop_blade_clear", "0.010");
    MassProps mp3 = massprops::compute(w2, cfg);
    // face = 0.60 + 0.0 - 0.010 = 0.590; root TE at x=0.60 → cl = 0.590-0.60 = -0.010
    CHECK(mp3.prop_clearance < 0.0);
}

// Prop keep-out: short unswept wing stays clear of the disk.
TEST(prop_clearance_ok_for_short_wing) {
    WingGeometry w = rect_wing();  // root_chord=0.20
    Config cfg;
    cfg.set("prop_diameter", "0.203");
    cfg.set("prop_hub_gap",  "0.010");
    cfg.set("prop_blade_clear", "0.005");
    // face = 0.20 + 0.010 - 0.005 = 0.205; root TE at x=0.20 → cl = 0.005 > 0
    MassProps mp = massprops::compute(w, cfg);
    CHECK(mp.prop_clearance > 0.0);
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
