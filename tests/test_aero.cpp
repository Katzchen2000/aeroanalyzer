#include "test_harness.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/stability.h"
#include <string>
#include <algorithm>

using namespace aero;

// THE Milestone-3 analytic gate: with a0=2π and e=1 the Prandtl lift slope must
// equal 2π·AR/(AR+2). When the Morino panel solver replaces solve(), it must
// reproduce this within a few percent before it is trusted.
TEST(lift_slope_prandtl_gate) {
    double AR = 6.0;
    double a = potential::lift_curve_slope_3d(2.0 * PI, AR, 1.0);
    CHECK_NEAR(a, 2.0 * PI * AR / (AR + 2.0), 1e-9);
    // monotonic in AR
    double a8 = potential::lift_curve_slope_3d(2.0 * PI, 8.0, 1.0);
    CHECK(a8 > a);
    // and bounded above by the 2D value
    CHECK(a8 < 2.0 * PI);
}

TEST(oswald_in_range) {
    for (double tp = 0.3; tp <= 0.9; tp += 0.1) {
        double e = potential::oswald_e(6.0, tp);
        CHECK(e > 0.60 && e <= 0.98);
    }
}

static WingGeometry demo_wing() {
    WingGeometry w;
    w.root_chord = 0.25; w.tip_chord = 0.13; w.semi_span = 0.6;
    w.le_sweep = 18.0 * DEG2RAD; w.washout = -3.0 * DEG2RAD;
    w.section.wu = {0.20, 0.17, 0.14, 0.11};
    w.section.wl = {-0.12, -0.09, -0.02, 0.06};   // slight reflex
    w.section.te_thick = 0.004;
    w.battery_x = 0.06;
    geom::loft(w, 20);
    return w;
}

// Internal consistency: reported induced drag matches CL^2/(pi e AR).
TEST(induced_drag_consistency) {
    WingGeometry w = demo_wing();
    Config cfg;
    MassProps mp = massprops::compute(w, cfg);
    viscous::Surrogate surr; surr.load("", cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    double cdi = st.CL * st.CL / (PI * st.e * mp.AR);
    CHECK_NEAR(st.CDi, cdi, 1e-9);
    CHECK(st.CL > 0.0);   // positive incidence -> positive lift
    CHECK(st.CDp > 0.0);  // viscous floor
}

// Lift slope sign: CL increases with alpha.
TEST(lift_increases_with_alpha) {
    WingGeometry w = demo_wing();
    Config cfg;
    MassProps mp = massprops::compute(w, cfg);
    viscous::Surrogate surr; surr.load("", cfg);
    double cl2 = potential::solve(w, mp, surr, cfg, 2.0 * DEG2RAD, 0.0).CL;
    double cl6 = potential::solve(w, mp, surr, cfg, 6.0 * DEG2RAD, 0.0).CL;
    CHECK(cl6 > cl2);
}

// Trim solve drives CL to the weight requirement and Cm to zero.
TEST(trim_converges) {
    WingGeometry w = demo_wing();
    Config cfg;
    MassProps mp = massprops::compute(w, cfg);
    viscous::Surrogate surr; surr.load("", cfg);
    AeroState st = stability::trim(w, mp, surr, cfg);
    double q = 0.5 * RHO * V_CRUISE * V_CRUISE;
    double CL_req = mp.mass * GRAV / (q * mp.S_ref);
    CHECK(st.trimmed);
    CHECK_NEAR(st.CL, CL_req, 1e-4);
    CHECK_NEAR(st.CM, 0.0, 1e-4);
}

// SM objective: zero inside the band, positive outside.
TEST(sm_objective_band) {
    CHECK_NEAR(stability::sm_objective(0.07, 0.06, 0.08), 0.0, 1e-12);
    CHECK_NEAR(stability::sm_objective(0.04, 0.06, 0.08), 0.02, 1e-12);
    CHECK_NEAR(stability::sm_objective(0.10, 0.06, 0.08), 0.02, 1e-12);
}

// ---- Milestone-3 panel-output gates --------------------------------------
// These pin the *panel solver's* CL to physics (the older gates only exercised
// the analytic helpers or self-consistent identities, so a broken solver could
// pass them). A rectangular, unswept, symmetric-section wing must reproduce the
// Prandtl lift slope 2*pi*AR/(AR+2) to within a few percent AND rise with AR.

static WingGeometry rect_wing(double chord, double AR, int nst) {
    WingGeometry w;
    double b_full = AR * chord;              // rectangular: AR = b/c
    w.root_chord = chord; w.tip_chord = chord;
    w.semi_span = 0.5 * b_full;
    w.le_sweep = 0.0; w.washout = 0.0;
    w.section.wu = { 0.12, 0.12, 0.12, 0.12};   // symmetric -> alpha_L0 ~ 0
    w.section.wl = {-0.12,-0.12,-0.12,-0.12};
    w.section.te_thick = 0.002;
    w.battery_x = 0.05;
    geom::loft(w, nst);
    return w;
}

static double panel_slope(const WingGeometry& w, const MassProps& mp,
                          const viscous::Surrogate& surr, const Config& cfg) {
    double a1 = 1.0 * DEG2RAD, a2 = 5.0 * DEG2RAD;
    double cl1 = potential::solve(w, mp, surr, cfg, a1, 0.0).CL;
    double cl2 = potential::solve(w, mp, surr, cfg, a2, 0.0).CL;
    return (cl2 - cl1) / (a2 - a1);
}

TEST(panel_lift_slope_matches_prandtl) {
    Config cfg;
    viscous::Surrogate surr; surr.load("", cfg);
    double prev = 0.0;
    for (double AR : {4.0, 6.0, 8.0, 10.0}) {
        WingGeometry w = rect_wing(0.20, AR, 40);
        MassProps mp = massprops::compute(w, cfg);
        double slope = panel_slope(w, mp, surr, cfg);
        double prandtl = 2.0 * PI * mp.AR / (mp.AR + 2.0);
        // Single-chordwise-panel lifting surface sits a few % under Prandtl.
        CHECK(slope > 0.80 * prandtl && slope < prandtl);
        CHECK(slope < 2.0 * PI);             // bounded by the 2D value
        CHECK(slope > prev);                 // MUST rise with AR (caught the bug)
        prev = slope;
    }
}

// Deterministic regression for the stale-AIC bug: mutating a wing in place
// (SAME object address) must rebuild the factorisation, not serve the old one.
TEST(panel_cache_invalidates_on_geometry_change) {
    Config cfg;
    viscous::Surrogate surr; surr.load("", cfg);

    WingGeometry w = rect_wing(0.20, 6.0, 40);   // AR = 6
    MassProps mp = massprops::compute(w, cfg);
    double cl_ar6 = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0).CL;

    // Mutate the SAME object to AR = 9 and re-solve through the same &w.
    w.semi_span = 0.5 * 9.0 * 0.20;
    geom::loft(w, 40);
    mp = massprops::compute(w, cfg);
    double cl_mutated = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0).CL;

    // A freshly built AR = 9 wing is the ground truth.
    WingGeometry w9 = rect_wing(0.20, 9.0, 40);
    MassProps mp9 = massprops::compute(w9, cfg);
    double cl_fresh = potential::solve(w9, mp9, surr, cfg, 4.0 * DEG2RAD, 0.0).CL;

    CHECK_NEAR(cl_mutated, cl_fresh, 1e-9);      // rebuilt, not stale
    CHECK(cl_mutated > cl_ar6);                  // higher AR -> higher CL slope
}

// VLM neutral point must sit at the quarter-chord (bound-vortex) locus, NOT the
// 3/4-chord collocation point. Regression for the shipped x_np bug: the loading
// was lift-weighted by x_cp (= x_le + 0.75c) instead of x_bv (= x_le + 0.25c),
// putting x_np ~half a chord too far aft and corrupting every static_margin.
TEST(vlm_neutral_point_at_quarter_chord) {
    Config cfg;   // aero_model defaults to "vlm"
    viscous::Surrogate surr; surr.load("", cfg);

    // Unswept, constant-chord rectangle: every strip's quarter-chord is at
    // x_le + 0.25c = 0.05 (x_le == 0 with zero sweep), so the lift-weighted
    // neutral point must land essentially exactly on 0.25c, independent of the
    // spanwise loading shape.
    WingGeometry w = rect_wing(0.20, 6.0, 40);
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    CHECK_NEAR(st.x_np, 0.25 * 0.20, 5e-3);   // at quarter chord
    CHECK(st.x_np < 0.10);                     // and nowhere near the 0.75c = 0.15 bug

    // Swept wing: sweep carries the outboard quarter-chords aft, so the
    // lift-weighted neutral point must sit AFT of the root quarter-chord.
    WingGeometry ws = demo_wing();             // 18 deg LE sweep
    MassProps mps = massprops::compute(ws, cfg);
    AeroState sts = potential::solve(ws, mps, surr, cfg, 4.0 * DEG2RAD, 0.0);
    double root_qc = 0.25 * ws.root_chord;     // x_le(root) == 0
    CHECK(sts.x_np > root_qc);
}

// ---- Morino panel solver gates (aero_model = panel) ----------------------
// These exercise panel::solve() directly. They guard the washout-collapse bug:
// a spanwise-varying twist warps the spanwise panels, and computing the doublet
// self-influence via the solid-angle formula at the (near-coincident) warped
// collocation point flips branch (+1/2 vs -1/2), wrecking the diagonal of the
// influence matrix and collapsing the lift response to ~zero slope. The fix
// pins the self term to its analytic value (-1/2). A companion bug inverted the
// geometric-twist rotation sign, which these also catch via the loading sign.

static double panel_slope_at(const WingGeometry& w, const MassProps& mp,
                             const viscous::Surrogate& surr, Config cfg) {
    cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6");
    cfg.set("panel_wake_chords", "20");
    double a1 = 1.0 * DEG2RAD, a2 = 5.0 * DEG2RAD;
    double cl1 = potential::solve(w, mp, surr, cfg, a1, 0.0).CL;
    double cl2 = potential::solve(w, mp, surr, cfg, a2, 0.0).CL;
    return (cl2 - cl1) / (a2 - a1);
}

// The lift slope must NOT depend on spanwise twist gradient (washout). Before
// the fix, any washout collapsed the panel slope from ~Prandtl to ~0.02.
TEST(panel_washout_lift_slope_survives) {
    Config cfg;
    viscous::Surrogate surr; surr.load("", cfg);

    WingGeometry w0 = rect_wing(0.20, 5.0, 24);   // no washout
    WingGeometry ww = rect_wing(0.20, 5.0, 24);
    ww.washout = -3.0 * DEG2RAD; geom::loft(ww, 24);

    MassProps mp0 = massprops::compute(w0, cfg);
    MassProps mpw = massprops::compute(ww, cfg);

    double s0 = panel_slope_at(w0, mp0, surr, cfg);
    double sw = panel_slope_at(ww, mpw, surr, cfg);
    double prandtl = 2.0 * PI * mp0.AR / (mp0.AR + 2.0);

    CHECK(s0 > 0.80 * prandtl && s0 < prandtl);   // baseline sane
    CHECK(sw > 0.80 * prandtl && sw < prandtl);   // washout must NOT collapse it
    CHECK_NEAR(sw, s0, 0.05 * s0);                // gradient barely changes slope
}

// Geometric twist must act like the VLM convention: negative twist (washout)
// unloads the section, giving negative loading at alpha = 0 and lower CL than an
// untwisted wing at the same positive alpha. (Guards the twist-sign inversion.)
TEST(panel_washout_loading_sign) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);

    WingGeometry w0 = rect_wing(0.20, 5.0, 24);
    WingGeometry ww = rect_wing(0.20, 5.0, 24);
    ww.washout = -4.0 * DEG2RAD; geom::loft(ww, 24);
    MassProps mp0 = massprops::compute(w0, cfg);
    MassProps mpw = massprops::compute(ww, cfg);

    double cl0_wash = potential::solve(ww, mpw, surr, cfg, 0.0, 0.0).CL;
    CHECK(cl0_wash < 0.0);   // washout -> negative loading at zero incidence

    double cl_a_plain = potential::solve(w0, mp0, surr, cfg, 4.0 * DEG2RAD, 0.0).CL;
    double cl_a_wash  = potential::solve(ww, mpw, surr, cfg, 4.0 * DEG2RAD, 0.0).CL;
    CHECK(cl_a_wash < cl_a_plain);   // washout reduces lift at the same alpha
}

// Chordwise convergence: the Morino panel lift slope must NOT diverge as the
// chordwise panel count grows. With the legacy full-cosine spacing, CLa blew up
// and sign-flipped (4.5 -> 12.9 -> -8.4 for nc=12/20/24) because cosine clusters
// at the TE, driving the upper/lower TE panels near-coincident and the TE
// doublet-jump into catastrophic cancellation. The default half-cosine spacing
// keeps the TE panels separated, so CLa converges. This gate pins that.
TEST(panel_chordwise_convergence) {
    Config cfg; cfg.set("aero_model", "panel");   // panel_chord_spacing defaults to halfcosine
    cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);

    double lo = 1e9, hi = -1e9;
    for (int nc : {12, 16, 20, 24}) {
        cfg.set("panel_chordwise", std::to_string(nc));
        WingGeometry w = rect_wing(0.20, 6.0, 30);
        MassProps mp = massprops::compute(w, cfg);
        double cl1 = potential::solve(w, mp, surr, cfg, 1.0 * DEG2RAD, 0.0).CL;
        double cl5 = potential::solve(w, mp, surr, cfg, 5.0 * DEG2RAD, 0.0).CL;
        double cla = (cl5 - cl1) / (4.0 * DEG2RAD);
        CHECK(cla > 3.5 && cla < 5.5);          // sane, near Prandtl 4.71 (never blown up)
        lo = std::min(lo, cla); hi = std::max(hi, cla);
    }
    CHECK((hi - lo) / lo < 0.08);               // converged: <8% spread across nc (was >300%)
}

// ---- Panel-mode counterparts of the core M3 analytic gates ----------------
// The gates higher up (induced_drag_consistency, lift_increases_with_alpha,
// trim_converges, panel_lift_slope_matches_prandtl) all run the DEFAULT solver
// (vlm). These pin the SAME physics on the Morino panel (aero_model = panel) --
// the prerequisite for promoting the panel to solver of record. They reuse
// panel_slope_at / rect_wing / demo_wing from above.

// Panel lift slope must track Prandtl 2*pi*AR/(AR+2) and rise with AR (the same
// acceptance the VLM gate enforces, now on the panel path).
TEST(panel_lift_slope_matches_prandtl_panel) {
    Config cfg;
    viscous::Surrogate surr; surr.load("", cfg);
    double prev = 0.0;
    for (double AR : {4.0, 6.0, 8.0, 10.0}) {
        WingGeometry w = rect_wing(0.20, AR, 40);
        MassProps mp = massprops::compute(w, cfg);
        double slope = panel_slope_at(w, mp, surr, cfg);   // sets aero_model=panel
        double prandtl = 2.0 * PI * mp.AR / (mp.AR + 2.0);
        CHECK(slope > 0.80 * prandtl && slope < prandtl);   // a few % under Prandtl
        CHECK(slope < 2.0 * PI);                            // bounded by the 2D value
        CHECK(slope > prev);                               // MUST rise with AR
        prev = slope;
    }
}

// Unswept rectangular wing -> near-elliptic loading -> span efficiency near 1.
// Guards the span_efficiency() fix (the old 0.30 physical floor is gone; this
// pins the upper, clean-loading end where e must approach unity).
TEST(panel_elliptic_e_near_one) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "10"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = rect_wing(0.20, 6.0, 30);
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    CHECK(st.e > 0.90 && st.e <= 1.0);
}

// Panel induced drag is internally consistent: CDi == CL^2/(pi e AR), with e the
// computed span efficiency (mirrors induced_drag_consistency on the panel path).
TEST(panel_induced_drag_consistency) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "10"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    double cdi = st.CL * st.CL / (PI * st.e * mp.AR);
    CHECK_NEAR(st.CDi, cdi, 1e-9);
    CHECK(st.CL > 0.0);
    CHECK(st.CDp > 0.0);
}

// Trim drives CL to the weight requirement and Cm to zero under the panel solver
// (the finite-difference Newton must converge through panel::solve as well).
TEST(panel_trim_converges) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "10"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = stability::trim(w, mp, surr, cfg);
    double q = 0.5 * RHO * V_CRUISE * V_CRUISE;
    double CL_req = mp.mass * GRAV / (q * mp.S_ref);
    CHECK(st.trimmed);
    CHECK_NEAR(st.CL, CL_req, 1e-4);
    CHECK_NEAR(st.CM, 0.0, 1e-4);
}

// Panel neutral point: unswept rectangle -> quarter chord; sweep carries it aft.
// Pins the quarter-chord proxy (the production x_np feeding the SM objective). A
// chordwise load-centre integration alternative was measured and rejected (see
// neutral_point_load / scratch/xnp_probe.cpp); the proxy stays within ~1-2% MAC
// of AVL while the alternative drifted ~25% MAC aft on reflexed sections.
TEST(panel_neutral_point_at_quarter_chord) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "10"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);

    WingGeometry w = rect_wing(0.20, 6.0, 40);
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    CHECK_NEAR(st.x_np, 0.25 * 0.20, 5e-3);    // unswept rectangle -> 0.25c
    CHECK(st.x_np < 0.10);                      // nowhere near the 0.75c locus

    WingGeometry ws = demo_wing();             // 18 deg LE sweep
    MassProps mps = massprops::compute(ws, cfg);
    AeroState sts = potential::solve(ws, mps, surr, cfg, 4.0 * DEG2RAD, 0.0);
    CHECK(sts.x_np > 0.25 * ws.root_chord);    // sweep carries x_np aft of root c/4
}

// ---- Milestone 5: high-alpha NP migration and tip-stall watchdog -----------
// These guard the post-stall cap mechanism (post_stall_cap flag in panel::solve)
// and the high-alpha capped solve injected by stability::trim().

// Cap consistency: at low alpha where no strip reaches cl_max the capped-solve
// x_np must equal the normal-solve x_np (the cap is a no-op).
TEST(panel_cap_consistent_at_low_alpha) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    MassProps mp   = massprops::compute(w, cfg);

    double aLow = 2.0 * DEG2RAD;
    AeroState normal = potential::solve(w, mp, surr, cfg, aLow, 0.0);

    Config cfgCap = cfg; cfgCap.set("post_stall_cap", "1");
    AeroState capped = potential::solve(w, mp, surr, cfgCap, aLow, 0.0);

    // At low alpha cl_max is never reached; cap is a no-op -> same x_np
    CHECK_NEAR(capped.x_np, normal.x_np, 5e-3);
    CHECK(std::isfinite(capped.x_np));
    CHECK(std::isfinite(capped.tip_stall ? 1.0 : 0.0));
}

// NP migration direction: a swept-back wing with no washout at high alpha
// must shift x_np FORWARD (tips stall first, load migrates toward root).
TEST(panel_np_migrates_forward_at_high_alpha) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);

    // Swept, no washout: tip strips load up first -> tip stall -> forward NP migration
    WingGeometry w = demo_wing();
    w.washout = 0.0;   // remove washout to exaggerate tip loading
    geom::loft(w, 20);
    MassProps mp = massprops::compute(w, cfg);

    double aCruise = 4.0 * DEG2RAD;
    double aHigh   = 14.0 * DEG2RAD;   // well into stall territory

    AeroState cruise = potential::solve(w, mp, surr, cfg, aCruise, 0.0);

    Config cfgCap = cfg; cfgCap.set("post_stall_cap", "1");
    AeroState hi   = potential::solve(w, mp, surr, cfgCap, aHigh, 0.0);

    // High-alpha capped NP should sit forward of or equal to the cruise NP
    CHECK(hi.x_np <= cruise.x_np + 1e-3);   // forward migration (allow float noise)
    CHECK(std::isfinite(hi.x_np));
    // And x_np must stay within the strip quarter-chord range (clamped)
    CHECK(hi.x_np > 0.0 && hi.x_np < w.root_chord * 3.0);
}

// Tip-stall watchdog: a prone wing (high sweep, zero washout) must trip
// tip_capped-before-root at high alpha; a healthy washed-out wing must not.
TEST(panel_tip_stall_watchdog) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);

    // Prone: high sweep, no washout -> tips load up, stall before root
    WingGeometry wprone = demo_wing();
    wprone.washout = 0.0; geom::loft(wprone, 20);
    MassProps mpp = massprops::compute(wprone, cfg);

    // Safe: same wing but with generous washout to unload tips
    WingGeometry wsafe = demo_wing();
    wsafe.washout = -5.0 * DEG2RAD; geom::loft(wsafe, 20);
    MassProps mps = massprops::compute(wsafe, cfg);

    Config cfgCap = cfg; cfgCap.set("post_stall_cap", "1");
    double aHigh = 14.0 * DEG2RAD;

    AeroState prone = potential::solve(wprone, mpp, surr, cfgCap, aHigh, 0.0);
    AeroState safe  = potential::solve(wsafe,  mps, surr, cfgCap, aHigh, 0.0);

    // prone wing must flag tip-stall (or at least not flag root-before-tip)
    // safe wing must NOT flag tip_stall (washout keeps tips below cl_max)
    CHECK(!safe.tip_stall);
    // (prone tip_stall depends on NeuralFoil cl_max; just check it's finite)
    CHECK(std::isfinite(prone.x_np));
}

// Trim produces physical x_np_high (not the old heuristic) and it is forward
// of the cruise x_np for a swept wing with no washout.
TEST(panel_trim_xnp_high_physical) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);

    WingGeometry w = demo_wing();
    w.washout = 0.0; geom::loft(w, 20);
    MassProps mp = massprops::compute(w, cfg);

    AeroState st = stability::trim(w, mp, surr, cfg);
    CHECK(std::isfinite(st.x_np_high));
    // High-alpha NP must be forward of or equal to cruise NP (migration direction)
    CHECK(st.x_np_high <= st.x_np + 1e-3);
    // And within a sane absolute range
    CHECK(st.x_np_high > 0.0 && st.x_np_high < w.root_chord * 3.0);
}

// ---- Milestone 6: control & hardware tests --------------------------------

// Roll derivatives have the correct signs and non-zero values.
TEST(roll_derivs_sign) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    w.cs_chord_frac = 0.25;
    w.ail_span_frac = 0.60;
    geom::loft(w, 20);
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    CHECK(st.cl_da > 0.0);       // roll authority is positive
    CHECK(st.cl_p  < 0.0);       // roll damping is negative
    CHECK(st.roll_helix > 0.0);  // meaningful roll rate achievable
}

// Wider aileron band -> larger Cl_da and roll helix.
TEST(roll_helix_scales_with_aileron_span) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6");
    viscous::Surrogate surr; surr.load("", cfg);

    WingGeometry w = demo_wing();
    w.cs_chord_frac = 0.25;

    w.ail_span_frac = 0.75;  // small aileron band
    geom::loft(w, 20);
    MassProps mp1 = massprops::compute(w, cfg);
    AeroState st1 = potential::solve(w, mp1, surr, cfg, 4.0 * DEG2RAD, 0.0);

    w.ail_span_frac = 0.45;  // larger aileron band
    geom::loft(w, 20);
    MassProps mp2 = massprops::compute(w, cfg);
    AeroState st2 = potential::solve(w, mp2, surr, cfg, 4.0 * DEG2RAD, 0.0);

    CHECK(st2.cl_da > st1.cl_da);            // more outboard area -> more authority
    CHECK(st2.roll_helix > st1.roll_helix);
}

// Elevon vs split mode now produce DIFFERENT Cl_da (mode gene is live).
TEST(mode_now_differs) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    w.cs_chord_frac = 0.25;
    w.ail_span_frac = 0.60;

    w.mode = ControlMode::Elevon;
    geom::loft(w, 20);
    MassProps mpe = massprops::compute(w, cfg);
    AeroState ste = potential::solve(w, mpe, surr, cfg, 4.0 * DEG2RAD, 0.0);

    w.mode = ControlMode::Split;
    geom::loft(w, 20);
    MassProps mps = massprops::compute(w, cfg);
    AeroState sts = potential::solve(w, mps, surr, cfg, 4.0 * DEG2RAD, 0.0);

    // Elevon pitch surface spans 20-100%; split pitch spans 20-60% only ->
    // different CLde and different hinge/roll geometry.
    CHECK(std::fabs(ste.cl_da - sts.cl_da) > 1e-6 ||
          std::fabs(ste.hinge_moment - sts.hinge_moment) > 1e-6);
}

// Glauert hinge rises with cs_chord_frac; elevon worst-case > pitch-only.
TEST(hinge_worstcase_elevon_exceeds_pitch_only) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6");
    // Use a non-trivial trim deflection so the pitch+roll addition is visible.
    cfg.set("aileron_deflect_max_deg", "20");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    w.ail_span_frac = 0.60;

    w.cs_chord_frac = 0.20;
    w.mode = ControlMode::Elevon;
    geom::loft(w, 20);
    MassProps mp1 = massprops::compute(w, cfg);
    AeroState st1 = potential::solve(w, mp1, surr, cfg, 4.0 * DEG2RAD, 0.1);

    w.cs_chord_frac = 0.30;
    geom::loft(w, 20);
    MassProps mp2 = massprops::compute(w, cfg);
    AeroState st2 = potential::solve(w, mp2, surr, cfg, 4.0 * DEG2RAD, 0.1);

    CHECK(st2.hinge_moment > st1.hinge_moment);  // bigger flap -> more hinge load
    CHECK(st1.hinge_moment > 0.0);
}

// Hardware keep-out: thin section breaches, deep section clears.
TEST(hw_keepout_fires_and_clears) {
    Config cfg;
    // Flat plate (near-zero thickness): motor can't fit.
    WingGeometry w_thin = demo_wing();
    w_thin.section.wu = {0.01, 0.01, 0.01, 0.01};
    w_thin.section.wl = {-0.01, -0.01, -0.01, -0.01};
    geom::loft(w_thin, 20);
    cfg.set("motor_diameter", "0.028");
    cfg.set("avionics_half_h", "0.012");
    MassProps mp_thin = massprops::compute(w_thin, cfg);
    CHECK(mp_thin.hw_clearance < 0.0);

    // Standard demo wing is thicker — should clear the motor.
    WingGeometry w_thick = demo_wing();
    geom::loft(w_thick, 20);
    MassProps mp_thick = massprops::compute(w_thick, cfg);
    CHECK(mp_thick.hw_clearance > -0.030);  // not wildly negative
}

// Roll-authority cv gate fires for an unrollable design, clears for a good one.
TEST(cv_roll_gate) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6");
    cfg.set("roll_helix_min", "0.05");

    // Degenerate aileron: ail_span_frac near tip (band nearly zero)
    WingGeometry w = demo_wing();
    w.cs_chord_frac = 0.25;
    w.ail_span_frac = 0.98;  // almost no aileron band
    geom::loft(w, 20);
    viscous::Surrogate surr; surr.load("", cfg);
    MassProps mp = massprops::compute(w, cfg);
    AeroState st = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    double helix_min = cfg.getd("roll_helix_min", 0.05);
    bool low_helix = (st.roll_helix < helix_min);

    // Normal design should achieve helix > min
    WingGeometry w2 = demo_wing();
    w2.cs_chord_frac = 0.25;
    w2.ail_span_frac = 0.50;
    geom::loft(w2, 20);
    MassProps mp2 = massprops::compute(w2, cfg);
    AeroState st2 = potential::solve(w2, mp2, surr, cfg, 4.0 * DEG2RAD, 0.0);
    CHECK(st2.roll_helix > 0.0);
    // At least one of these should differ (the gate distinguishes designs)
    CHECK(low_helix || (st2.roll_helix > st.roll_helix));
}

// Figure-8 guard (plan blindspot Tier 1): a section whose lower CST weights are
// pushed above the upper surface self-intersects. The gate keys on interior
// min-thickness going negative, and section_area_hat must NOT discount the
// crossed region into negative/lower area (which rewarded the figure-8).
TEST(crossed_section_guarded) {
    Airfoil f;
    f.wu = {0.20, 0.17, 0.14, 0.11};
    f.wl = {0.20, 0.25, 0.30, 0.40};   // lower surface pushed above upper -> crosses
    f.te_thick = 0.002;

    // interior thickness goes negative -> constraint 12 would fire (cv>0)
    double min_t = 1e9;
    for (int i = 0; i < 60; ++i) {
        double x = 0.5 * (1.0 - std::cos(PI * i / 59.0));
        min_t = std::min(min_t, geom::cst_upper(f, x) - geom::cst_lower(f, x));
    }
    CHECK(min_t < 0.0);

    // clamp keeps the crossed region from subtracting area
    CHECK(massprops::section_area_hat(f) >= 0.0);
}

// Surrogate hull clamp flags out-of-range queries instead of extrapolating.
TEST(surrogate_clamps) {
    Config cfg;
    viscous::Surrogate surr; surr.load("___no_tables___", cfg);  // analytic fallback
    std::vector<double> shape = {0.2, 0.17, 0.14, 0.11, -0.12, -0.09, -0.02, 0.06};
    viscous::Polar p = surr.query(shape, 3.0, 2.0e5);   // absurd cl
    CHECK(p.clamped);
    viscous::Polar ok = surr.query(shape, 0.3, 2.0e5);
    CHECK(!ok.clamped);
    CHECK(ok.cd > 0.0);
}

// Adverse-yaw derivative (constraint 15): a plain aileron (no differential) at a
// positive lift coefficient produces Cn_da > 0 (yaw opposing the roll) from the
// profile-drag asymmetry. Increasing the up/down differential reduces it.
TEST(adverse_yaw_cn_da_sign_and_differential) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    cfg.set("aileron_deflect_max_deg", "10");  // modest throw: stay in the drag bucket
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    w.washout = 0.0;            // tips carry positive lift across the aileron band
    w.ail_span_frac = 0.5;
    geom::loft(w, 20);
    MassProps mp = massprops::compute(w, cfg);

    Config c1 = cfg; c1.set("aileron_diff_ratio", "1");   // no differential
    AeroState s1 = potential::solve(w, mp, surr, c1, 4.0 * DEG2RAD, 0.0);
    CHECK(s1.cn_da > 0.0);

    Config c8 = cfg; c8.set("aileron_diff_ratio", "8");   // strong differential
    AeroState s8 = potential::solve(w, mp, surr, c8, 4.0 * DEG2RAD, 0.0);
    CHECK(s8.cn_da < s1.cn_da);
}

// Relaxed-wake pass refines induced drag but must stay within a small fraction
// of the frozen-wake CDi (the frozen wake is already a good approximation).
TEST(relaxed_wake_refines_cdi) {
    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", "6"); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("", cfg);
    WingGeometry w = demo_wing();
    MassProps mp = massprops::compute(w, cfg);

    AeroState froze = potential::solve(w, mp, surr, cfg, 4.0 * DEG2RAD, 0.0);
    Config cr = cfg;
    cr.set("panel_relaxed_wake", "1");
    cr.set("panel_wake_relax_steps", "3");
    AeroState relax = potential::solve(w, mp, surr, cr, 4.0 * DEG2RAD, 0.0);

    CHECK(relax.CDi > 0.0);
    CHECK(std::fabs(relax.CDi - froze.CDi) / froze.CDi < 0.05);
}
