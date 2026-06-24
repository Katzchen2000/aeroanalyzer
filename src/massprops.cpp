#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/geom.h"
#include <cmath>
#include <algorithm>

namespace aero {
namespace massprops {

double section_area_hat(const Airfoil& f) {
    const int N = 120;
    double area = 0.0, prev_t = 0.0, prev_x = 0.0;
    for (int i = 0; i < N; ++i) {
        double th = PI * i / (N - 1);
        double x = 0.5 * (1.0 - std::cos(th));
        // clamp: a crossed (figure-8) section must not subtract area, else the
        // mass objective rewards self-intersection (plan blindspot Tier 1B).
        double t = std::max(0.0, geom::cst_upper(f, x) - geom::cst_lower(f, x));
        if (i > 0) area += 0.5 * (t + prev_t) * (x - prev_x);
        prev_t = t; prev_x = x;
    }
    return area;
}

double section_perimeter_hat(const Airfoil& f) {
    const int N = 120;
    double per = 0.0;
    double pux = 0, puy = 0, plx = 0, ply = 0;
    for (int i = 0; i < N; ++i) {
        double th = PI * i / (N - 1);
        double x = 0.5 * (1.0 - std::cos(th));
        double yu = geom::cst_upper(f, x), yl = geom::cst_lower(f, x);
        if (i > 0) {
            per += std::hypot(x - pux, yu - puy);  // upper surface
            per += std::hypot(x - plx, yl - ply);  // lower surface
        }
        pux = x; puy = yu; plx = x; ply = yl;
    }
    return per;
}

MassProps compute(const WingGeometry& w, const Config& cfg) {
    MassProps mp;
    geom::planform(w, mp.S_ref, mp.mac, mp.x_mac_le, mp.b_full, mp.AR);

    const double rho_mat  = cfg.getd("material_density", 1070.0);
    const double t_shell  = cfg.getd("shell_thickness", 0.0012);
    const double infill_r = cfg.getd("infill_root", 0.10);
    const double infill_t = cfg.getd("infill_tip", 0.03);

    // ---- structural mass + CG via spanwise trapezoidal accumulation ----
    double vol_half = 0.0;     // m^3 (one half)
    double struct_m = 0.0;     // kg (one half)
    double struct_mx = 0.0;    // kg*m
    for (const auto& s : w.stations) {
        double t = (w.semi_span > 0) ? s.y / w.semi_span : 0.0;
        double A_hat = section_area_hat(s.af);       // chord=1 area, lofted section
        double P_hat = section_perimeter_hat(s.af);  // chord=1 perimeter
        double infill = infill_r + (infill_t - infill_r) * t;     // gradient
        double enclosed = A_hat * s.chord * s.chord;              // m^2
        double shell = P_hat * s.chord * t_shell;                // m^2 (thin wall)
        shell = std::min(shell, enclosed);
        double solid = shell + infill * (enclosed - shell);      // equiv solid area
        double dvol = solid * s.width;                           // m^3
        double dm = dvol * rho_mat;                              // kg
        double x_centroid = s.x_le + 0.42 * s.chord;            // ~area centroid
        vol_half  += dvol;
        struct_m  += dm;
        struct_mx += dm * x_centroid;
    }
    mp.volume = 2.0 * vol_half;

    // ---- point masses (plan §2, §4) ----
    auto chord_at = [&](double t) {
        return w.root_chord + (w.tip_chord - w.root_chord) * t;
    };
    auto xle_at = [&](double t) { return t * w.semi_span * std::tan(w.le_sweep); };

    double m_motor = cfg.getd("mass_motor", 0.060);
    double x_motor = w.root_chord;                       // Y=0 trailing edge
    double m_batt  = cfg.getd("mass_battery", 0.210);
    double x_batt  = w.battery_x + 0.035;                // box CG (70mm long)
    double m_avi   = cfg.getd("target_mass_aux", 0.045); // 20%-50% span block
    double t_avi   = 0.35;
    double x_avi   = xle_at(t_avi) + 0.30 * chord_at(t_avi);
    double m_srv   = cfg.getd("mass_servo_each", 0.012);
    double t_srv   = 0.30;
    double x_srv   = xle_at(t_srv) + 0.72 * chord_at(t_srv);

    double m_struct_full = 2.0 * struct_m;
    double mx_struct_full = 2.0 * struct_mx;

    double m_pts = m_motor + m_batt + m_avi + 2.0 * m_srv;
    double mx_pts = m_motor * x_motor + m_batt * x_batt + m_avi * x_avi
                  + 2.0 * m_srv * x_srv;

    mp.mass = m_struct_full + m_pts;
    mp.x_cg = (mp.mass > 0) ? (mx_struct_full + mx_pts) / mp.mass : 0.0;

    // ---- analytic spar clearance (plan §4) ----
    // clearance = local half-thickness at the spar's chordwise station - radius.
    // Negative => OML breach. No B-Rep booleans.
    //
    // NOTE (modeling decision surfaced by the optimizer): a *straight* spanwise
    // spar held at a constant x = 15% root chord exits a swept wing almost
    // immediately (outboard the LE has marched aft past that x), so it breaches
    // the OML on every swept candidate. The buildable default is therefore a
    // spar that follows the swept 15%-local-chord line. Set spar_straight=1 in
    // the config to model a literal straight carbon tube instead.
    double root_frac = cfg.getd("spar_root_frac", 0.15);
    double r_spar = 0.5 * cfg.getd("spar_diameter", 0.006);
    bool straight = cfg.getb("spar_straight", false);
    double x_spar_const = root_frac * w.root_chord;
    double min_clear = 1.0;
    for (const auto& s : w.stations) {
        double frac;
        if (straight)
            frac = (s.chord > 0) ? (x_spar_const - s.x_le) / s.chord : -1.0;
        else
            frac = root_frac;                 // spar follows the sweep
        double clear;
        if (frac < 0.02 || frac > 0.98) {
            clear = -(r_spar);   // spar exits the section near LE/TE
        } else {
            double half_t = 0.5 * (geom::cst_upper(s.af, frac)
                                   - geom::cst_lower(s.af, frac)) * s.chord;
            clear = half_t - r_spar;
        }
        min_clear = std::min(min_clear, clear);
    }
    mp.spar_clearance = min_clear;

    // ---- Hardware keep-out clearances (M6) --------------------------------
    // Mirrors the spar clearance pattern: reuse cst_upper/cst_lower for local
    // half-thickness, compare against the hardware's space requirement.
    // ponytail: battery-fit check omitted (plan names motor + avionics); add when
    //           battery box dimensions become a design variable.

    // ponytail: motor (root) + avionics keep-outs use the root section; the
    // lofted spanwise variation is second-order for these fixed-spot clearance
    // checks. Use the lofted station section if avionics placement gets tuned.
    // Motor at root TE: check half-thickness at 80% chord >= motor radius.
    double motor_r   = 0.5 * cfg.getd("motor_diameter", 0.028);
    double half_t_mot = 0.5 * (geom::cst_upper(w.section, 0.80)
                               - geom::cst_lower(w.section, 0.80)) * w.root_chord;
    double mot_clear = half_t_mot - motor_r;

    // Avionics block at t=0.35, 30% chord: check half-thickness >= avionics_half_h.
    double avi_hh    = cfg.getd("avionics_half_h", 0.012);
    double t_avi_hw  = 0.35;
    double chord_avi = w.root_chord + (w.tip_chord - w.root_chord) * t_avi_hw;
    double half_t_avi = 0.5 * (geom::cst_upper(w.section, 0.30)
                                - geom::cst_lower(w.section, 0.30)) * chord_avi;
    double avi_clear = half_t_avi - avi_hh;

    mp.hw_clearance = std::min(mot_clear, avi_clear);

    return mp;
}

}  // namespace massprops
}  // namespace aero
