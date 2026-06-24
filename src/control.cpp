// control.cpp — Shared flight-control derivatives (M6).
//
// Replaces the verbatim-duplicated pitch_control / flap_tau / hinge blocks in
// aero_potential.cpp and aero_panel.cpp with one implementation.
//
// Physics summary:
//   Pitch : Glauert flap theory — CLde = a*tau*(Sf/S), Cmde = -CLde*arm
//   Roll  : strip theory — Cl_da = 2*a*tau*ail_ym/(S*b), Cl_p = -4a/(S*b^2)*Idamp
//   Hinge : thin-airfoil Ch_alpha, Ch_delta (exact Glauert, replaces heuristic)

#include "aeroanalyzer/control.h"
#include <cmath>
#include <algorithm>

namespace aero {
namespace control {

namespace {

// Glauert flap effectiveness tau = 1 - (tf - sin tf)/pi,
// tf = arccos(2*cf - 1).  Moved here from both aero backends.
static double flap_tau(double cf) {
    if (cf <= 0.0) return 0.0;
    if (cf >= 1.0) return 1.0;
    double tf = std::acos(2.0 * cf - 1.0);
    return 1.0 - (tf - std::sin(tf)) / PI;
}

// Thin-airfoil hinge-moment magnitudes (both positive; caller applies sign).
// Using tf = arccos(2*cf - 1) (the same "TE-origin" angle as flap_tau):
//   ch_alpha = (sin tf - tf*cos tf) / pi      [derived: theta_f = pi-tf substituted]
//   ch_delta = ((pi-tf) + sin tf*cos tf) / pi
// Numerically verified vs textbook values at cf=0.30: ch_alpha~0.544, ch_delta~0.252.
static void flap_hinge_coeffs(double cf, double& ch_alpha, double& ch_delta) {
    if (cf <= 0.0) { ch_alpha = 0.0; ch_delta = 0.0; return; }
    if (cf >= 1.0) { ch_alpha = 0.0; ch_delta = 0.0; return; }
    double tf = std::acos(2.0 * cf - 1.0);
    double s  = std::sin(tf);
    double c  = std::cos(tf);
    ch_alpha = (s - tf * c) / PI;           // > 0
    ch_delta = ((PI - tf) + s * c) / PI;    // > 0
}

}  // namespace

Derivs compute(const WingGeometry& w, const MassProps& mp,
               double a, double alpha, double delta_e, const Config& cfg) {
    Derivs d{};
    const double b = mp.b_full > 0 ? mp.b_full : 2.0 * w.semi_span;
    const double S = mp.S_ref  > 0 ? mp.S_ref  : 1.0;

    const double cf      = w.cs_chord_frac;  // control-surface chord fraction
    const double ail_t   = w.ail_span_frac;  // aileron inboard edge (semi-span frac)
    const double tau     = flap_tau(cf);

    // ---- pitch surface extent (depends on mode) --------------------------
    //   Elevon: full-span surface [0.20, 1.0], symmetric for pitch.
    //   Split : inboard elevator [0.20, ail_t], pitch only (disjoint from ailerons).
    const double pitch_y0 = 0.20;
    const double pitch_y1 = (w.mode == ControlMode::Elevon) ? 1.0 : ail_t;

    // Accumulate integrals over stations (right half-wing)
    double Sf_half  = 0.0;  // pitch surface area (one side)
    double arm_num  = 0.0;  // area-weighted CP x, for x_cp
    double cfc_sum  = 0.0;  // mean control-surface chord numerator
    double n_pitch  = 0.0;
    double ail_ym   = 0.0;  // Σ cf*c_i*y_i*Δy_i over aileron band (right side)
    double roll_damp = 0.0; // Σ c_i*y_i²*Δy_i over all strips (right side)

    for (const auto& s : w.stations) {
        const double t = (w.semi_span > 0) ? s.y / w.semi_span : 0.0;

        // Pitch surface
        if (t >= pitch_y0 && t <= pitch_y1) {
            Sf_half += cf * s.chord * s.width;
            double xcp = s.x_le + (1.0 - 0.5 * cf) * s.chord;
            arm_num   += xcp * (cf * s.chord * s.width);
            cfc_sum   += cf * s.chord;
            n_pitch   += 1.0;
        }

        // Aileron band: [ail_t, 1.0] in both modes
        if (t >= ail_t) {
            ail_ym += cf * s.chord * s.y * s.width;
        }

        // Roll damping: all strips (no flap factor — this is full-chord damp)
        roll_damp += s.chord * s.y * s.y * s.width;
    }

    // ---- pitch derivatives ----------------------------------------------
    const double Sf = 2.0 * Sf_half;
    d.Sf      = Sf;
    d.cf_chord = (n_pitch > 0) ? cfc_sum / n_pitch : cf * w.root_chord;
    d.x_cp    = (Sf_half > 0) ? arm_num / Sf_half : w.root_chord;

    const double CLde = a * tau * (S > 0 ? Sf / S : 0.0);
    d.CLde = CLde;
    const double arm  = (mp.mac > 0) ? (d.x_cp - mp.x_cg) / mp.mac : 0.0;
    d.Cmde = -CLde * arm;

    // ---- roll control derivative ----------------------------------------
    // Cl_da = 2*a*tau * Σ(cf*c_i*y_i*Δy_i) / (S*b)  (both halves, factor 2)
    d.Cl_da = (S > 0 && b > 0) ? 2.0 * a * tau * ail_ym / (S * b) : 0.0;

    // ---- roll damping ---------------------------------------------------
    // Cl_p = -4*a / (S*b^2) * Σ(c_i*y_i^2*Δy_i)  (always negative)
    d.Cl_p = (S > 0 && b > 0) ? -4.0 * a * roll_damp / (S * b * b) : 0.0;
    if (d.Cl_p >= 0.0) d.Cl_p = -1.0e-6;  // guard degenerate geometry

    // ---- roll helix with 4:1 differential (both halves antisymmetric) ---
    // pb/2V = -Cl_da * da_eff / Cl_p  (positive when well-designed)
    // ponytail: 4:1 differential reduces adverse yaw (not modeled here);
    // it sets the deflection budget. Upgrade path: add Cn_da/Cn_p yaw model.
    const double da_max    = cfg.getd("aileron_deflect_max_deg", 20.0) * DEG2RAD;
    const double diff_ratio = cfg.getd("aileron_diff_ratio", 4.0);
    const double da_eff    = da_max * (1.0 + 1.0 / diff_ratio);
    d.roll_helix = std::max(0.0, -d.Cl_da * da_eff / d.Cl_p);

    // ---- hinge moment — Glauert thin-airfoil (replaces 0.45|δ|+0.05|α|) -
    double ch_a, ch_d;
    flap_hinge_coeffs(cf, ch_a, ch_d);  // both positive magnitudes

    // Worst-case combined deflection:
    //   Elevon: one surface does pitch AND full roll simultaneously.
    //   Split : elevator sees pitch only, ailerons see roll only -> take max.
    const double d_worst = (w.mode == ControlMode::Elevon)
        ? std::fabs(delta_e) + da_max
        : std::max(std::fabs(delta_e), da_max);

    const double V    = cfg.getd("v_cruise", V_CRUISE);
    const double q    = 0.5 * RHO * V * V;
    const double Ch   = ch_a * std::fabs(alpha) + ch_d * d_worst;
    const double H_Nm = Ch * q * Sf * d.cf_chord;
    d.hinge_moment    = H_Nm * (100.0 / GRAV);  // N·m -> kg·cm

    return d;
}

double adverse_yaw_cn_da(const WingGeometry& w, const MassProps& mp,
                         const viscous::Surrogate& surr,
                         const std::vector<double>& cl_local,
                         double a, const Config& cfg) {
    const double S = mp.S_ref > 0 ? mp.S_ref : 1.0;
    const double b = mp.b_full > 0 ? mp.b_full : 2.0 * w.semi_span;
    const double V = cfg.getd("v_cruise", V_CRUISE);
    const double cosL = std::cos(w.le_sweep);
    const double tau  = flap_tau(w.cs_chord_frac);
    const double da_max = cfg.getd("aileron_deflect_max_deg", 20.0) * DEG2RAD;
    const double diff   = std::max(1.0, cfg.getd("aileron_diff_ratio", 4.0));
    const double d_down = da_max / diff;   // down-going aileron (differential: smaller)
    const double d_up   = da_max;          // up-going aileron (larger throw)
    const double ail_t  = w.ail_span_frac;

    // num = integral of (cd_down - cd_up) * chord * y * dy over the aileron band.
    double num = 0.0;
    std::vector<double> shape;
    for (std::size_t i = 0; i < w.stations.size(); ++i) {
        const Station& s = w.stations[i];
        double t = (w.semi_span > 0) ? s.y / w.semi_span : 0.0;
        if (t < ail_t) continue;
        double cl_base = (i < cl_local.size()) ? cl_local[i] : 0.0;
        double dcl_down = a * tau * d_down;
        double dcl_up   = a * tau * d_up;
        shape.clear();
        shape.insert(shape.end(), s.af.wu.begin(), s.af.wu.end());
        shape.insert(shape.end(), s.af.wl.begin(), s.af.wl.end());
        double Re = RHO * V * cosL * s.chord / MU;
        double cd_down = surr.query(shape, (cl_base + dcl_down) / (cosL * cosL),
                                    Re, s.af.te_thick).cd;
        double cd_up   = surr.query(shape, (cl_base - dcl_up) / (cosL * cosL),
                                    Re, s.af.te_thick).cd;
        num += (cd_down - cd_up) * s.chord * s.y * s.width;
    }
    return (2.0 / (S * b)) * num;
}

}  // namespace control
}  // namespace aero
