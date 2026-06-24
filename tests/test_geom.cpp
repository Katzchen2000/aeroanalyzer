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
    w.section.wu = {0.18, 0.15, 0.12, 0.10};
    w.section.wl = {-0.18, -0.15, -0.12, -0.10};
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
    geom::loft(w, 20);
    CHECK_NEAR(w.stations.front().y, 0.0, 1e-9);
    CHECK_NEAR(w.stations.back().y, 0.6, 1e-9);
    CHECK_NEAR(w.stations.back().twist, -4.0 * DEG2RAD, 1e-9);
    CHECK(w.stations.back().x_le > w.stations.front().x_le);  // swept aft
}

// Variable-section genome: root + tip CST (4th order each) -> 24 genes, tip
// bounds mirror the root bounds.
TEST(genome_has_tip_section) {
    CHECK(geom::N_GENES == 26);   // 18 base genes + 8 tip CST weights
    geom::GenomeSpec spec = geom::default_genome();
    CHECK((int)spec.size() == 26);
    CHECK_NEAR(spec.lo[geom::G_TIP_WU0], spec.lo[geom::G_WU0], 1e-12);
    CHECK_NEAR(spec.hi[geom::G_TIP_WL3], spec.hi[geom::G_WL3], 1e-12);
}

// loft() blends root->tip CST linearly in eta; endpoints reproduce the sections
// exactly and the tip TE closes when section_tip.te_thick = 0.
TEST(loft_blends_root_to_tip) {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.13; w.semi_span = 0.6;
    w.section.wu = {0.20, 0.17, 0.14, 0.11};
    w.section.wl = {-0.12, -0.09, -0.02, 0.06};
    w.section.te_thick = 0.006;
    w.section_tip.wu = {0.10, 0.09, 0.08, 0.05};
    w.section_tip.wl = {-0.06, -0.04, 0.00, 0.03};
    w.section_tip.te_thick = 0.0;                 // closed tip
    geom::loft(w, 21);
    for (int j = 0; j < 4; ++j) {
        CHECK_NEAR(w.stations.front().af.wu[j], w.section.wu[j], 1e-9);
        CHECK_NEAR(w.stations.back().af.wu[j],  w.section_tip.wu[j], 1e-9);
    }
    CHECK_NEAR(w.stations.back().af.te_thick, 0.0, 1e-9);   // tip TE closed
    const Station& s = w.stations[10];
    double t = s.y / w.semi_span;
    for (int j = 0; j < 4; ++j)
        CHECK_NEAR(s.af.wu[j],
                   (1.0 - t) * w.section.wu[j] + t * w.section_tip.wu[j], 1e-9);
}
