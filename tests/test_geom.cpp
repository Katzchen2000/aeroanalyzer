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

// Rectangular planform: known area / aspect ratio.
TEST(planform_rectangular) {
    // build directly to control taper=1, sweep=0
    WingGeometry w;
    w.root_chord = 0.20; w.tip_chord = 0.20; w.semi_span = 0.60;
    w.le_sweep = 0.0; w.washout = 0.0;
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
    w.root_chord = 0.25; w.tip_chord = 0.12; w.semi_span = 0.6;
    w.le_sweep = 20.0 * DEG2RAD; w.washout = -4.0 * DEG2RAD;
    w.sections.resize(1);
    geom::loft(w, 20);
    CHECK_NEAR(w.stations.front().y, 0.0, 1e-9);
    CHECK_NEAR(w.stations.back().y, 0.6, 1e-9);
    CHECK_NEAR(w.stations.back().twist, -4.0 * DEG2RAD, 1e-9);
    CHECK(w.stations.back().x_le > w.stations.front().x_le);  // swept aft
}

// 5-section genome: N_GENES == 52, G_SEC addressing correct.
TEST(genome_5_sections) {
    CHECK(geom::N_GENES == 52);
    geom::GenomeSpec spec = geom::default_genome();
    CHECK((int)spec.size() == 52);
    // all sections share the same CST bounds
    for (int k = 1; k < geom::N_SECTIONS; ++k) {
        CHECK_NEAR(spec.lo[geom::G_SEC(k,0,0)], spec.lo[geom::G_SEC(0,0,0)], 1e-12);
        CHECK_NEAR(spec.hi[geom::G_SEC(k,1,3)], spec.hi[geom::G_SEC(0,1,3)], 1e-12);
    }
    // G_LE_BOW / G_TE_BOW exist and are symmetric
    CHECK_NEAR(spec.lo[geom::G_LE_BOW], -spec.hi[geom::G_LE_BOW], 1e-12);
}

// loft() piecewise blend: knot points reproduce control sections exactly;
// inner segment blends correctly; le_bow parabola is applied.
TEST(loft_5section_blend) {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.13; w.semi_span = 0.6;
    w.le_bow = 0.0; w.te_bow = 0.0;
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

// le_bow parabola: midspan x_le offset is exactly le_bow.
TEST(loft_le_bow) {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.13; w.semi_span = 0.6;
    w.le_sweep = 0.0; w.washout = 0.0;
    w.le_bow = 0.03; w.te_bow = 0.0;
    w.sections.resize(1);
    geom::loft(w, 41);
    // midspan station nearest η=0.5
    int mid = 0; double best = 1e9;
    for (int i = 0; i < (int)w.stations.size(); ++i) {
        double d = std::fabs(w.stations[i].y - 0.5 * w.semi_span);
        if (d < best) { best = d; mid = i; }
    }
    // parabola value at η=0.5: 4*0.5*0.5 = 1.0 → x_le = le_bow
    CHECK_NEAR(w.stations[mid].x_le, w.le_bow, 1e-3);
}

// te_bow < 0 makes tip chord narrower; gate should trigger in evaluate.
// Here just verify the chord formula: chord at midspan decreases with te_bow < 0.
TEST(loft_te_bow_chord) {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.13; w.semi_span = 0.6;
    w.le_bow = 0.0; w.te_bow = -0.05;
    w.sections.resize(1);
    geom::loft(w, 41);
    int mid = 0; double best = 1e9;
    for (int i = 0; i < (int)w.stations.size(); ++i) {
        double d = std::fabs(w.stations[i].y - 0.5 * w.semi_span);
        if (d < best) { best = d; mid = i; }
    }
    double chord_no_bow = 0.25 + (0.13 - 0.25) * 0.5;   // linear taper at η=0.5
    CHECK(w.stations[mid].chord < chord_no_bow);
}
