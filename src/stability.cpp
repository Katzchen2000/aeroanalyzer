#include "aeroanalyzer/stability.h"
#include "aeroanalyzer/aero_potential.h"
#include <cmath>
#include <algorithm>

namespace aero {
namespace stability {

double sm_objective(double sm, double lo, double hi) {
    if (sm < lo) return lo - sm;
    if (sm > hi) return sm - hi;
    return 0.0;
}

AeroState trim(const WingGeometry& w, const MassProps& mp,
               const viscous::Surrogate& surr, const Config& cfg) {
    const double V = cfg.getd("v_cruise", V_CRUISE);
    const double q = 0.5 * RHO * V * V;
    const double W = mp.mass * GRAV;
    const double CL_req = (q * mp.S_ref > 0) ? W / (q * mp.S_ref) : 0.0;

    double alpha = 2.0 * DEG2RAD;   // initial guess
    double delta = 0.0;

    // Wake-freeze strategy for the dense panel AIC (the GA hot-path cost driver:
    // a cold AIC build is ~hundreds of ms, and the cache key includes wake_alpha
    // so any alpha change forces a full rebuild). Two independent knobs:
    //   panel_trim_freeze_jac  (default 1): the FD Jacobian probes (alpha+h /
    //       delta+h) reuse the residual factorization -- a cleaner derivative
    //       (no wake-rebuild discontinuity in dCL/dalpha) and 1 fewer build/iter.
    //   panel_trim_freeze_wake (default 1): the residual solves freeze the wake
    //       too, so the WHOLE Newton solve uses ONE factorization (built at the
    //       initial alpha) instead of rebuilding every iteration -- the dominant
    //       saving. The wake inclination's effect on attached loading is second
    //       order, so the trim point is within a fraction of a percent of the
    //       live-wake result (validated in scratch/panel_timing.cpp). Both knobs
    //       are inert for the VLM model, which has no wake cache.
    // Setting either to 0 restores the legacy rebuild-every-solve behaviour.
    bool freezeWake = cfg.geti("panel_trim_freeze_wake", 1) != 0;
    bool freezeJac  = cfg.geti("panel_trim_freeze_jac", 1) != 0;
    Config cfgRes = cfg, cfgJac = cfg;
    if (freezeWake) { cfgRes.set("panel_freeze_wake", "1"); cfgJac.set("panel_freeze_wake", "1"); }
    else if (freezeJac) { cfgJac.set("panel_freeze_wake", "1"); }

    AeroState st = potential::solve(w, mp, surr, cfgRes, alpha, delta);
    const int max_iter = 50;
    const double tol = 1e-8;
    const double h = 1e-6;          // FD step (rad)

    for (int it = 0; it < max_iter; ++it) {
        st = potential::solve(w, mp, surr, cfgRes, alpha, delta);
        double r0 = st.CL - CL_req;
        double r1 = st.CM;
        if (std::fabs(r0) < tol && std::fabs(r1) < tol) {
            st.trimmed = true;
            break;
        }
        // finite-difference Jacobian J = d[r]/d[alpha,delta] (frozen-wake probes)
        AeroState sa = potential::solve(w, mp, surr, cfgJac, alpha + h, delta);
        AeroState sd = potential::solve(w, mp, surr, cfgJac, alpha, delta + h);
        double j00 = (sa.CL - st.CL) / h;   // dCL/dalpha
        double j01 = (sd.CL - st.CL) / h;   // dCL/ddelta
        double j10 = (sa.CM - st.CM) / h;   // dCm/dalpha
        double j11 = (sd.CM - st.CM) / h;   // dCm/ddelta
        double det = j00 * j11 - j01 * j10;
        if (std::fabs(det) < 1e-12) break;  // singular -> leave untrimmed
        double da = (-r0 * j11 + r1 * j01) / det;
        double dd = (-j00 * r1 + j10 * r0) / det;
        // step clamp for robustness against the nonlinear panel model later
        double maxstep = 5.0 * DEG2RAD;
        da = std::max(-maxstep, std::min(maxstep, da));
        dd = std::max(-maxstep, std::min(maxstep, dd));
        alpha += da;
        delta += dd;
    }
    st = potential::solve(w, mp, surr, cfgRes, alpha, delta);
    st.trimmed = st.trimmed ||
        (std::fabs(st.CL - CL_req) < 1e-5 && std::fabs(st.CM) < 1e-5);

    // ---- High-alpha capped solve: physical NP migration and tip-stall (M5) -
    // Reuses the warm frozen-wake AIC (same geometry, one back-sub + strip loop).
    // Only for the panel model; VLM keeps the heuristic from solve() above.
    if (cfg.gets("aero_model", "panel") == "panel") {
        double aHi = cfg.getd("high_alpha_deg", 12.0) * DEG2RAD;
        Config cfgHi = cfgRes;          // inherits panel_freeze_wake=1 when set
        cfgHi.set("post_stall_cap", "1");
        AeroState hi = potential::solve(w, mp, surr, cfgHi, aHi, delta);
        st.x_np_high = hi.x_np;        // capped-load centroid at high alpha
        st.tip_stall = hi.tip_stall;   // tips-before-root stall
    }
    return st;
}

}  // namespace stability
}  // namespace aero
