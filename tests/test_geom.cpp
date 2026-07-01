#include "test_harness.h"
#include "aeroanalyzer/geom.h"

using namespace aero;

// CST round-trip: fit must recover the surfaces it was sampled from (plan §3
// "ancestry injection"). If this fails, low-order CST cannot represent the
// shapes the GA explores.
TEST(cst_round_trip) {
    Airfoil f;
    f.wu = {0.20, 0.18, 0.15, 0.12};
    f.wl = {-0.15, -0.12, -0.05, 0.08};
    f.te_thick = 0.003;

    std::vector<std::pair<double, double>> up, lo;
    for (int i = 0; i <= 40; ++i) {
        double x = 0.5 * (1.0 - std::cos(PI * i / 40));
        up.push_back({x, geom::cst_upper(f, x)});
        lo.push_back({x, geom::cst_lower(f, x)});
    }
    Airfoil g = geom::fit_cst(up, lo, 3, f.te_thick);
    for (int i = 1; i <= 10; ++i) {
        double x = i / 11.0;
        CHECK_NEAR(geom::cst_upper(g, x), geom::cst_upper(f, x), 1e-4);
        CHECK_NEAR(geom::cst_lower(g, x), geom::cst_lower(f, x), 1e-4);
    }
}

// Symmetric section -> zero camber -> zero-lift angle and cm_ac ~ 0.
TEST(symmetric_airfoil_zero_camber) {
    Airfoil f;
    f.wu = {0.20, 0.15, 0.10, 0.08};
    f.wl = {-0.20, -0.15, -0.10, -0.08};
    f.te_thick = 0.0;
    geom::ThinAirfoil t = geom::thin_airfoil(f);
    CHECK_NEAR(t.alpha_L0, 0.0, 1e-3);
    CHECK_NEAR(t.cm_ac, 0.0, 1e-3);
}

// Positive camber -> negative zero-lift angle (textbook sign).
TEST(cambered_airfoil_negative_alpha0) {
    Airfoil f;
    f.wu = {0.26, 0.23, 0.20, 0.18};
    f.wl = {-0.06, -0.03, 0.00, 0.03};   // mean camber line clearly positive
    f.te_thick = 0.002;
    geom::ThinAirfoil t = geom::thin_airfoil(f);
    CHECK(t.alpha_L0 < 0.0);
}

// Aft-camber change must move cm_ac (reflex authority exists).
TEST(reflex_changes_cm) {
    Airfoil a; a.wu = {0.20, 0.16, 0.12, 0.10}; a.wl = {-0.16, -0.12, -0.06, -0.04};
    Airfoil b = a; b.wl[3] = 0.12;   // raise aft lower surface -> reflex
    double cma = geom::thin_airfoil(a).cm_ac;
    double cmb = geom::thin_airfoil(b).cm_ac;
    CHECK(std::fabs(cmb - cma) > 1e-3);
}

// Rectangular planform via set_linear_planform: known area / aspect ratio.
TEST(planform_rectangular) {
    WingGeometry w;
    w.semi_span = 0.60;
    geom::set_linear_planform(w, 0.20, 0.20, 0.0, 0.0);
    w.sections.resize(1);   // single uniform section, loft handles it
    geom::loft(w, 20);
    double S, mac, xle, b, AR;
    geom::planform(w, S, mac, xle, b, AR);
    CHECK_NEAR(S, 2.0 * 0.60 * 0.20, 1e-6);   // 0.24 m^2
    CHECK_NEAR(b, 1.20, 1e-9);
    CHECK_NEAR(AR, 1.44 / 0.24, 1e-6);        // 6.0
    CHECK_NEAR(mac, 0.20, 1e-6);
}

// Lofting respects cosine spacing endpoints and monotonic span.
TEST(loft_spanwise) {
    WingGeometry w;
    w.semi_span = 0.6;
    geom::set_linear_planform(w, 0.25, 0.12, 20.0 * DEG2RAD, -4.0 * DEG2RAD);
    w.sections.resize(1);
    geom::loft(w, 20);
    CHECK_NEAR(w.stations.front().y, 0.0, 1e-9);
    CHECK_NEAR(w.stations.back().y, 0.6, 1e-9);
    CHECK_NEAR(w.stations.back().twist, -4.0 * DEG2RAD, 1e-9);
    CHECK(w.stations.back().x_le > w.stations.front().x_le);  // swept aft
}

// New genome layout: N_GENES == 67, G_SEC addressing correct.
TEST(genome_layout_67_genes) {
    CHECK(geom::N_GENES == 67);
    geom::GenomeSpec spec = geom::default_genome();
    CHECK((int)spec.size() == 67);
    // all sections share the same CST bounds
    for (int k = 1; k < geom::N_SECTIONS; ++k) {
        CHECK_NEAR(spec.lo[geom::G_SEC(k,0,0)], spec.lo[geom::G_SEC(0,0,0)], 1e-12);
        CHECK_NEAR(spec.hi[geom::G_SEC(k,1,3)], spec.hi[geom::G_SEC(0,1,3)], 1e-12);
    }
}

// loft() piecewise CST blend: knot points reproduce control sections exactly;
// inner segment blends correctly. Unchanged behavior from the old 5-section loft.
TEST(loft_5section_blend) {
    WingGeometry w;
    w.semi_span = 0.6;
    geom::set_linear_planform(w, 0.25, 0.13, 0.0, 0.0);
    // 5 control sections with distinct wu[0] values
    w.sections.resize(5);
    for (int k = 0; k < 5; ++k) {
        w.sections[k].wu = {0.10 + 0.04*k, 0.15, 0.12, 0.09};
        w.sections[k].wl = {-0.10, -0.08, -0.04, 0.02};
        w.sections[k].te_thick = (k == 4) ? 0.0 : 0.005;  // closed tip
    }
    // n=100: cos(pi*66/99) = cos(2pi/3) = -0.5 → y=0.45 = 0.75*b exactly
    geom::loft(w, 100);

    // η=0 root station == s0
    for (int j = 0; j < 4; ++j)
        CHECK_NEAR(w.stations.front().af.wu[j], w.sections[0].wu[j], 1e-9);
    // η=1 tip station == s4
    for (int j = 0; j < 4; ++j)
        CHECK_NEAR(w.stations.back().af.wu[j], w.sections[4].wu[j], 1e-9);
    CHECK_NEAR(w.stations.back().af.te_thick, 0.0, 1e-9);

    // station 66 lands at η=0.75 exactly; must equal s2
    double target = 0.75 * w.semi_span;
    int idx = 0;
    double best = 1e9;
    for (int i = 0; i < (int)w.stations.size(); ++i) {
        double d = std::fabs(w.stations[i].y - target);
        if (d < best) { best = d; idx = i; }
    }
    CHECK(best < 1e-6);
    for (int j = 0; j < 4; ++j)
        CHECK_NEAR(w.stations[idx].af.wu[j], w.sections[2].wu[j], 1e-9);
}

// decode() forces all section te_thick = 0 (sharp TE everywhere).
TEST(decode_sharp_te) {
    geom::GenomeSpec spec = geom::default_genome();
    // use lo-bound genome (all values at lower bound = valid, minimal wing)
    std::vector<double> genes(spec.size());
    for (int i = 0; i < (int)spec.size(); ++i) genes[i] = spec.lo[i];
    WingGeometry w = geom::decode(genes, spec);
    for (int k = 0; k < (int)w.sections.size(); ++k)
        CHECK_NEAR(w.sections[k].te_thick, 0.0, 1e-12);
    // CST at x=1: class_fn(x=1) = 0 so cst_upper/lower both zero regardless of weights
    for (int k = 0; k < (int)w.sections.size(); ++k) {
        double gap = geom::cst_upper(w.sections[k], 1.0) - geom::cst_lower(w.sections[k], 1.0);
        CHECK_NEAR(gap, 0.0, 1e-9);
    }
}

// ---- new organic-Bezier evaluator tests ---------------------------------

// Endpoint exactness: Bezier curves interpolate their first/last control point.
TEST(bezier_endpoint_exactness) {
    WingGeometry w;
    w.semi_span = 0.6;
    w.chord_cp = {0.28, 0.24, 0.22, 0.18, 0.15, 0.11};
    w.sweep_cp = {0.0, 0.10, 0.20, 0.30, 0.40, 0.55};
    w.twist_cp = {1.0*DEG2RAD, 0.5*DEG2RAD, -0.5*DEG2RAD, -2.0*DEG2RAD, -3.5*DEG2RAD, -5.0*DEG2RAD};
    w.dih_cp   = {0.0, 2.0*DEG2RAD, 5.0*DEG2RAD, 10.0*DEG2RAD, 20.0*DEG2RAD, 35.0*DEG2RAD, 50.0*DEG2RAD};

    CHECK_NEAR(geom::chord_at(w, 0.0), w.chord_cp.front(), 1e-12);
    CHECK_NEAR(geom::chord_at(w, 1.0), w.chord_cp.back(), 1e-12);
    CHECK_NEAR(geom::twist_at(w, 0.0), w.twist_cp.front(), 1e-12);
    CHECK_NEAR(geom::twist_at(w, 1.0), w.twist_cp.back(), 1e-12);
    CHECK_NEAR(geom::xle_at(w, 0.0), 0.0, 1e-12);
    CHECK_NEAR(geom::xle_at(w, 1.0), w.semi_span * w.sweep_cp.back(), 1e-9);
    CHECK_NEAR(geom::dihedral_at(w, 0.0), 0.0, 1e-12);
    CHECK_NEAR(geom::dihedral_at(w, 1.0), w.dih_cp.back(), 1e-12);
}

// Smoothness / no crease: a rising dihedral curve near the tip must be C1 —
// finite-difference slope has no jump anywhere, unlike the old winglet fold.
TEST(dihedral_curve_no_crease) {
    WingGeometry w;
    w.semi_span = 0.6;
    geom::set_linear_planform(w, 0.25, 0.13, 0.0, 0.0);
    w.dih_cp = {0.0, 0.0, 2.0*DEG2RAD, 8.0*DEG2RAD, 20.0*DEG2RAD, 40.0*DEG2RAD, 60.0*DEG2RAD};

    const int N = 400;
    const double h = 1.0 / N;
    double max_slope_jump = 0.0;
    double prev_slope = 0.0;
    for (int i = 1; i <= N; ++i) {
        double eta = i * h;
        double slope = (geom::dihedral_at(w, eta) - geom::dihedral_at(w, eta - h)) / h;
        if (i > 1) {
            double jump = std::fabs(slope - prev_slope);
            if (jump > max_slope_jump) max_slope_jump = jump;
        }
        prev_slope = slope;
    }
    // A hard crease (like the old winglet fold) produces a slope jump of tens of
    // rad/unit-eta at the fold station; a smooth Bezier stays tiny between samples.
    CHECK(max_slope_jump < 1.0);
}

// Chord positivity: convex-hull guarantee — all-positive control points keep
// chord_at(eta) > 0 across the whole span, even for an aggressive random curve.
TEST(chord_positivity_convex_hull) {
    WingGeometry w;
    w.semi_span = 0.6;
    w.chord_cp = {0.30, 0.05, 0.28, 0.04, 0.20, 0.06};  // wild but all > 0
    for (int i = 0; i <= 50; ++i) {
        double eta = i / 50.0;
        CHECK(geom::chord_at(w, eta) > 0.0);
    }
}

// Twist non-linearity: a curved twist_cp set must NOT reduce to a straight
// line between root and tip (guards against silently falling back to the
// old purely-linear washout model).
TEST(twist_curve_is_nonlinear) {
    WingGeometry w;
    w.semi_span = 0.6;
    w.twist_cp = {0.0, 15.0*DEG2RAD, -20.0*DEG2RAD, 5.0*DEG2RAD, -15.0*DEG2RAD, -5.0*DEG2RAD};
    double t0 = geom::twist_at(w, 0.0);
    double t1 = geom::twist_at(w, 1.0);
    double mid_linear = 0.5 * (t0 + t1);
    double mid_actual = geom::twist_at(w, 0.5);
    CHECK(std::fabs(mid_actual - mid_linear) > 1.0 * DEG2RAD);
}

// Single source of truth: chord_at/xle_at/twist_at/dihedral_at must match the
// values loft() actually wrote into w.stations at that station's eta (guards
// massprops/avl_export drift back to independently re-derived geometry).
TEST(loft_matches_evaluators) {
    WingGeometry w;
    w.semi_span = 0.6;
    w.chord_cp = {0.28, 0.24, 0.22, 0.18, 0.15, 0.11};
    w.sweep_cp = {0.0, 0.10, 0.20, 0.30, 0.40, 0.55};
    w.twist_cp = {1.0*DEG2RAD, 0.5*DEG2RAD, -0.5*DEG2RAD, -2.0*DEG2RAD, -3.5*DEG2RAD, -5.0*DEG2RAD};
    w.dih_cp   = {0.0, 2.0*DEG2RAD, 5.0*DEG2RAD, 10.0*DEG2RAD, 20.0*DEG2RAD, 35.0*DEG2RAD, 50.0*DEG2RAD};
    w.sections.resize(1);
    geom::loft(w, 30);

    for (const auto& s : w.stations) {
        CHECK_NEAR(s.chord, geom::chord_at(w, s.eta), 1e-9);
        CHECK_NEAR(s.x_le, geom::xle_at(w, s.eta), 1e-9);
        CHECK_NEAR(s.twist, geom::twist_at(w, s.eta), 1e-9);
    }
}

// Non-planar height scales with tip dihedral: a steeper rising tip curve
// produces a larger nonplanar_h (drives the induced-drag credit); a flat
// wing (all dih_cp = 0) produces nonplanar_h ~ 0.
TEST(nonplanar_height_scales_with_tip_dihedral) {
    WingGeometry flat;
    flat.semi_span = 0.6;
    geom::set_linear_planform(flat, 0.25, 0.13, 0.0, 0.0);
    flat.sections.resize(1);
    geom::loft(flat, 20);
    CHECK_NEAR(flat.nonplanar_h, 0.0, 1e-9);

    WingGeometry gentle = flat;
    gentle.dih_cp = {0.0, 2.0*DEG2RAD, 4.0*DEG2RAD, 6.0*DEG2RAD, 8.0*DEG2RAD, 10.0*DEG2RAD, 12.0*DEG2RAD};
    geom::loft(gentle, 20);

    WingGeometry steep = flat;
    steep.dih_cp = {0.0, 5.0*DEG2RAD, 15.0*DEG2RAD, 30.0*DEG2RAD, 45.0*DEG2RAD, 60.0*DEG2RAD, 75.0*DEG2RAD};
    geom::loft(steep, 20);

    CHECK(gentle.nonplanar_h > flat.nonplanar_h);
    CHECK(steep.nonplanar_h > gentle.nonplanar_h);
}
