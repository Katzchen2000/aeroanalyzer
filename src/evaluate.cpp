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

EvalResult Evaluator::run(const std::vector<double>& genes) const {
    EvalResult r;
    r.geom = geom::decode(genes, spec_);
    geom::loft(r.geom, n_stations_);
    r.mp = massprops::compute(r.geom, cfg_);
    r.aero = stability::trim(r.geom, r.mp, surr_, cfg_);

    const double V = cfg_.getd("v_cruise", V_CRUISE);
    const double q = 0.5 * RHO * V * V;

    // ---- objectives (all minimized) ----
    r.objectives[OBJ_DRAG] = r.aero.CD * q * r.mp.S_ref;     // drag force, N
    r.objectives[OBJ_MASS] = r.mp.mass;                      // kg
    double sm_lo = cfg_.getd("sm_band_lo", 0.06);
    double sm_hi = cfg_.getd("sm_band_hi", 0.08);
    r.objectives[OBJ_SM] = stability::sm_objective(r.aero.static_margin,
                                                   sm_lo, sm_hi);

    // ---- constraints -> total violation cv (0 = feasible) ----
    double cv = 0.0;

    // (1) TE thickness clamp (plan §3)
    double te_min = cfg_.getd("te_thick_min_mm", 0.8) / 1000.0;
    double te_actual = r.geom.section.te_thick * r.geom.root_chord;
    if (te_actual < te_min) cv += 50.0 * (te_min - te_actual) / te_min;

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

    // (11) static-margin band floor: the OBJ_SM band is only an objective, so
    // the drag-minimizing corner of the front sits below sm_lo (less washout =>
    // less trim drag). This hard cv gate forces every feasible design into the
    // >= sm_band_lo region. sm_floor_penalty = 0 recovers the soft-only behavior.
    double sm_floor_pen = cfg_.getd("sm_floor_penalty", 30.0);
    if (sm_floor_pen > 0.0 && r.aero.static_margin < sm_lo)
        cv += sm_floor_pen * (sm_lo - r.aero.static_margin) / sm_lo;

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
    return run(genes);
}

}  // namespace aero
