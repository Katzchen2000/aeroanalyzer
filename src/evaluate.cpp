#include "aeroanalyzer/evaluate.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/stability.h"
#include <cmath>
#include <algorithm>

namespace aero {

Evaluator::Evaluator(const Config& cfg) : cfg_(cfg) {
    spec_ = geom::default_genome();
    seeds_ = seeds::load_seeds(cfg_);
    seeds::widen_cst_bounds(spec_, seeds_, cfg_.getd("cst_bound_margin", 0.04));
    surr_.load("data/surrogates", cfg_);
    n_stations_ = cfg_.geti("n_stations", 20);
}

EvalResult Evaluator::run(const std::vector<double>& genes, bool relaxed_wake) const {
    EvalResult r;
    r.geom = geom::decode(genes, spec_);
    r.geom.winglet_taper = cfg_.getd("winglet_taper", 1.0);
    r.geom.winglet_blend = cfg_.getd("winglet_blend", 0.06);
    geom::loft(r.geom, n_stations_);
    r.mp = massprops::compute(r.geom, cfg_);
    // reporting re-eval (detail) turns on the relaxed-wake pass for exact CDi;
    // the GA hot path keeps the frozen wake.
    if (relaxed_wake) {
        Config c = cfg_;
        c.set("panel_relaxed_wake", "1");
        r.aero = stability::trim(r.geom, r.mp, surr_, c);
    } else {
        r.aero = stability::trim(r.geom, r.mp, surr_, cfg_);
    }

    const double V = cfg_.getd("v_cruise", V_CRUISE);
    const double q = 0.5 * RHO * V * V;

    // ---- objectives (all minimized) ----
    // OBJ_DRAG = CD/CL = inverse L/D (pure aero efficiency, weight-agnostic).
    // Mass owns the weight axis; drag-force (mass·g/(L/D)) would double-count it.
    {
        double cl_eff = std::max(0.05, r.aero.CL);
        r.objectives[OBJ_DRAG] = r.aero.CD / cl_eff;
    }
    r.objectives[OBJ_MASS] = r.mp.mass;                      // kg
    double sm_lo = cfg_.getd("sm_band_lo", 0.06);
    double sm_hi = cfg_.getd("sm_band_hi", 0.08);
    r.objectives[OBJ_SM] = stability::sm_objective(r.aero.static_margin,
                                                   sm_lo, sm_hi);

    // ---- constraints -> total violation cv (0 = feasible) ----
    double cv = 0.0;

    // (1) TE is always sharp (te_thick=0 from decode); no te gate needed.

    // (2) tip Reynolds floor (plan §3)
    double re_min = cfg_.getd("re_tip_min", 100000.0);
    double Re_tip = RHO * V * r.geom.tip_chord / MU;
    if (Re_tip < re_min) cv += 5.0 * (re_min - Re_tip) / re_min;
    // (3) spar-to-OML clearance: continuous penalty under 1mm (plan §4)
    double clear_target = 0.001;
    if (r.mp.spar_clearance < clear_target) {
        double d = (clear_target - r.mp.spar_clearance) / clear_target;
        cv += (r.mp.spar_clearance < 0.0 ? 30.0 : 4.0) * d;  // breach is worse
    }

    // (4) hinge-moment gate, fatal (plan §7)
    double hinge_max = cfg_.getd("hinge_moment_max", 1.2);
    if (r.aero.hinge_moment > hinge_max)
        cv += 40.0 * (r.aero.hinge_moment - hinge_max) / hinge_max;

    // (5) high-alpha NP must stay aft of CG, else delete (plan §7)
    if (r.aero.x_np_high < r.mp.x_cg) {
        double d = (r.mp.x_cg - r.aero.x_np_high) / (r.mp.mac > 0 ? r.mp.mac : 1);
        cv += 60.0 * d;
    }

    // (6) tip-stall watchdog at cruise (plan §6)
    if (r.aero.tip_stall) cv += 20.0;

    // (7) trim must converge
    if (!r.aero.trimmed) cv += 100.0;

    // (8) basic static stability floor (NP aft of CG at cruise)
    if (r.aero.static_margin < 0.0) cv += 25.0 * (-r.aero.static_margin);

    // (9) roll authority: steady helix pb/2V must meet minimum (M6)
    double helix_min = cfg_.getd("roll_helix_min", 0.05);
    if (helix_min > 0.0 && r.aero.roll_helix < helix_min)
        cv += 30.0 * (helix_min - r.aero.roll_helix) / helix_min;
    // (10) hardware keep-out: motor + avionics must fit in section (M6)
    if (r.mp.hw_clearance < 0.0) {
        double hw_ref = cfg_.getd("motor_diameter", 0.028) * 0.5;
        if (hw_ref <= 0.0) hw_ref = 0.014;
        cv += 30.0 * (-r.mp.hw_clearance) / hw_ref;
    }
    // (18) pusher propeller keep-out: wing TE must stay clear of prop disk
    if (r.mp.prop_clearance < 0.0) {
        double prop_r = 0.5 * cfg_.getd("prop_diameter", 0.203);
        if (prop_r > 0.0) cv += 30.0 * (-r.mp.prop_clearance) / prop_r;
    }

    // (11) static-margin band floor: the OBJ_SM band is only an objective, so
    // the drag-minimizing corner of the front sits below sm_lo (less washout =>
    // less trim drag). This hard cv gate forces every feasible design into the
    // >= sm_band_lo region. sm_floor_penalty = 0 recovers the soft-only behavior.
    double sm_floor_pen = cfg_.getd("sm_floor_penalty", 30.0);
    if (sm_floor_pen > 0.0 && r.aero.static_margin < sm_lo)
        cv += sm_floor_pen * (sm_lo - r.aero.static_margin) / sm_lo;

    // (12) interior section validity: lower surface must stay below upper across
    // all K control sections. The piecewise loft spans between them so checking
    // every section bounds the full wing.
    double t_floor = cfg_.getd("min_thickness_frac", 0.005);  // 0.5% chord
    double min_t = 1e9;
    for (const auto& af : r.geom.sections) {
        for (int i = 0; i < 60; ++i) {
            double x = 0.5 * (1.0 - std::cos(PI * i / 59.0));
            // Skip LE/TE: CST thickness ~0 by construction; TE policed by #1.
            if (x < 0.05 || x > 0.95) continue;
            min_t = std::min(min_t, geom::cst_upper(af, x) - geom::cst_lower(af, x));
        }
    }
    if (min_t < t_floor) cv += 40.0 * (t_floor - min_t) / t_floor;

    // (14) chord-collapse: catches extreme exp/gull combos that make tip chord negative.
    // Winglet stations get a smaller floor (they're lightly loaded and can be narrow).
    double chord_min_m    = cfg_.getd("chord_min_m", 0.03);
    double chord_min_wl_m = cfg_.getd("chord_min_winglet_m", 0.015);
    for (const auto& s : r.geom.stations) {
        double floor_m = s.in_winglet ? chord_min_wl_m : chord_min_m;
        if (s.chord < floor_m) cv += 50.0 * (floor_m - s.chord) / floor_m;
    }

    // (13) surrogate confidence: keep OOD polars from silently passing.
    double conf_threshold = cfg_.getd("confidence_threshold", 0.5);
    if (conf_threshold > 0.0 && r.aero.polar_confidence < conf_threshold)
        cv += 10.0 * (conf_threshold - r.aero.polar_confidence) / conf_threshold;

    // (15) adverse-yaw protection: penalize positive Cn_da (aileron yaw opposing
    // the roll). Differential throw (aileron_diff_ratio) is the design lever.
    // Configurable like sm_floor_penalty: adverse yaw is a matter of degree (most
    // plain ailerons are mildly adverse), so the weight tunes how hard the GA is
    // pushed toward differential-friendly designs; 0 disables the gate.
    double adv_pen = cfg_.getd("adverse_yaw_penalty", 50.0);
    if (adv_pen > 0.0 && r.aero.cn_da > 0.0) cv += adv_pen * r.aero.cn_da;

    // (16) Dutch-roll damping floor (hard constraint, same pattern as SM floor).
    // dutch_roll_zeta_min default 0.05 — "good enough" per user requirement.
    // Set to 0 to disable (reverts to reporting-only behaviour).
    double zdr_min = cfg_.getd("dutch_roll_zeta_min", 0.05);
    if (zdr_min > 0.0 && r.aero.dutch_roll_zeta < zdr_min) {
        double dyn_pen = cfg_.getd("dynamic_stab_penalty", 5.0);
        cv += dyn_pen * (zdr_min - r.aero.dutch_roll_zeta);
    }
    // phugoid floor remains opt-in (typically stable by Lanchester)
    double dyn_pen = cfg_.getd("dynamic_stab_penalty", 5.0);
    double zph_min = cfg_.getd("phugoid_zeta_min", 0.0);
    if (dyn_pen > 0.0 && zph_min > 0.0 && r.aero.phugoid_zeta < zph_min)
        cv += dyn_pen * (zph_min - r.aero.phugoid_zeta);

    // (17) Banked-turn gates (n·CL; no extra panel solve — same deriv model as cruise)
    // Directional divergence in the turn
    if (r.aero.cn_beta_turn <= 0.0)
        cv += 30.0 * (-r.aero.cn_beta_turn + 1e-4);
    // Dutch-roll in the turn
    if (zdr_min > 0.0 && r.aero.dutch_roll_zeta_turn < zdr_min)
        cv += dyn_pen * (zdr_min - r.aero.dutch_roll_zeta_turn);
    // Tip-stall in the turn
    if (r.aero.tip_stall_turn) cv += 20.0;

    r.cv = cv;
    return r;
}

void Evaluator::evaluate(Candidate& c) const {
    EvalResult r = run(c.genes);
    c.objectives = r.objectives;
    c.cv = r.cv;
    c.evaluated = true;
}

EvalResult Evaluator::detail(const std::vector<double>& genes) const {
    return run(genes, /*relaxed_wake=*/true);
}

}  // namespace aero
