#include "aeroanalyzer/stability.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/control.h"
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
    // The relaxed-wake pass tilts (and rebuilds) the wake; running it inside the
    // frozen-wake Newton iterations would corrupt the cached factorization. Strip
    // it here and apply exactly one pass after convergence (below).
    cfgRes.set("panel_relaxed_wake", "0");
    cfgJac.set("panel_relaxed_wake", "0");

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

    // ---- One relaxed-wake pass at the trimmed point (opt-in; reporting only) --
    // Refines induced drag with a free-aligned wake. Panel model only.
    if (cfg.geti("panel_relaxed_wake", 0) != 0 &&
        cfg.gets("aero_model", "panel") == "panel") {
        Config cfgRx = cfg;
        cfgRx.set("panel_relaxed_wake", "1");
        cfgRx.set("panel_freeze_wake", "0");   // live wake for the realignment
        AeroState rx = potential::solve(w, mp, surr, cfgRx, alpha, delta);
        st.e = rx.e; st.CDi = rx.CDi; st.CD = rx.CD;
        st.polar_confidence = std::min(st.polar_confidence, rx.polar_confidence);
    }

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
        st.polar_confidence = std::min(st.polar_confidence, hi.polar_confidence);
    }

    // ---- Dutch-roll helper (used twice: cruise + banked-turn) ---------------
    struct DrResult { double omega = 0.0, zeta = 0.0; };
    auto dutch_roll = [&](const control::LateralDerivs& ld, double cl_p_val,
                          double CL_val, double CD_val) -> DrResult {
        if (mp.Izz <= 0.0 || mp.Ixx <= 0.0 || ld.cn_beta <= 0.0)
            return {0.0, (ld.cn_beta <= 0.0) ? -1.0 : 0.0};
        const double qSb  = 0.5 * RHO * V * V * mp.S_ref * mp.b_full;
        const double qSb2 = qSb * mp.b_full / (2.0 * V);
        double Lb = qSb  * ld.cl_beta  / mp.Ixx;
        double Lp = qSb2 * cl_p_val    / mp.Ixx;
        double Lr = qSb2 * ld.cl_r     / mp.Ixx;
        double Nb = qSb  * ld.cn_beta  / mp.Izz;
        double Nr = qSb2 * ld.cn_r     / mp.Izz;
        double Np = qSb2 * ld.cn_p     / mp.Izz;
        double omega2      = (std::fabs(Lp) > 1e-9) ? Nb - (Lb * Np) / Lp : Nb;
        double two_zeta_om = (std::fabs(Lp) > 1e-9) ? -Nr + (Lr * Np) / Lp : -Nr;
        if (omega2 <= 0.0) return {0.0, -1.0};
        DrResult dr;
        dr.omega = std::sqrt(omega2);
        dr.zeta  = two_zeta_om / (2.0 * dr.omega);
        return dr;
    };

    // ---- Dynamic stability: cruise point ------------------------------------
    {
        auto ld = control::lateral_derivs(w, mp, st.CL, st.CD, cfg);
        st.cn_beta = ld.cn_beta;
        st.cn_r    = ld.cn_r;

        auto dr = dutch_roll(ld, st.cl_p, st.CL, st.CD);
        st.dutch_roll_omega = dr.omega;
        st.dutch_roll_zeta  = dr.zeta;

        // Phugoid: Lanchester approximation  zeta_ph = CD / (sqrt(2) * CL)
        st.phugoid_zeta = (st.CL > 1e-9) ? st.CD / (std::sqrt(2.0) * st.CL) : 0.0;
    }

    // ---- Banked-turn prediction (n·CL; zero extra panel solve) --------------
    // ponytail: CDi ∝ CL² at fixed span efficiency; CDp unchanged.
    {
        double n      = cfg.getd("load_factor_turn", 2.0);
        double CL_t   = n * st.CL;
        double CDi_t  = (st.CL > 1e-9) ? st.CDi * (n * n) : st.CDi;
        double CD_t   = st.CDp + CDi_t;
        auto ld_t = control::lateral_derivs(w, mp, CL_t, CD_t, cfg);
        st.cn_beta_turn = ld_t.cn_beta;
        auto dr_t = dutch_roll(ld_t, st.cl_p, CL_t, CD_t);
        st.dutch_roll_zeta_turn = dr_t.zeta;

        // Tip-stall proxy: outboard station cl scales with n
        double cl_max_t = cfg.getd("cl_max_turn", 1.2);
        if (!st.cl_local.empty()) {
            // check outermost 25% of span
            int n_sta = (int)st.cl_local.size();
            int i_ail = (int)(0.75 * n_sta);
            for (int i = i_ail; i < n_sta; ++i) {
                if (st.cl_local[i] * n > cl_max_t) { st.tip_stall_turn = true; break; }
            }
        }
    }

    return st;
}

}  // namespace stability
}  // namespace aero
