#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/geom.h"
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>

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
    // strip cache for Izz (computed after x_cg is known)
    struct StripMass { double dm, x_cen, y, z; };
    std::vector<StripMass> strips;
    strips.reserve(w.stations.size());
    for (const auto& s : w.stations) {
        double t = (w.semi_span > 0) ? s.y / w.semi_span : 0.0;
        double A_hat = section_area_hat(s.af);       // chord=1 area, lofted section
        double P_hat = section_perimeter_hat(s.af);  // chord=1 perimeter
        double infill = infill_r + (infill_t - infill_r) * t;     // gradient
        double enclosed = A_hat * s.chord * s.chord;              // m^2
        double shell = P_hat * s.chord * t_shell;                // m^2 (thin wall)
        shell = std::min(shell, enclosed);
        double solid = shell + infill * (enclosed - shell);      // equiv solid area
        double dvol = solid * s.ds;                              // m^3 (true arc length)
        double dm = dvol * rho_mat;                              // kg
        double x_centroid = s.x_le + 0.42 * s.chord;            // ~area centroid
        vol_half  += dvol;
        struct_m  += dm;
        struct_mx += dm * x_centroid;
        strips.push_back({dm, x_centroid, s.y, s.z});
    }
    mp.volume = 2.0 * vol_half;

    // ---- point masses (plan §2, §4) ----
    // Single source of truth: read the smooth spanwise curves via geom::*_at,
    // never re-derive chord/x_le locally (that duplication used to silently
    // drift from the wing the aero solver actually lofted).
    auto chord_at = [&](double t) { return geom::chord_at(w, t); };
    auto xle_at   = [&](double t) { return geom::xle_at(w, t); };

    double m_motor = cfg.getd("mass_motor", 0.060);
    double x_motor = w.root_chord;                       // Y=0 trailing edge
    double m_batt   = cfg.getd("mass_battery",    0.210);
    double batt_len = cfg.getd("battery_len_m",   0.070);
    double batt_wid = cfg.getd("battery_width_m", 0.035);
    double batt_hgt = cfg.getd("battery_height_m",0.015);
    double batt_x0  = std::min(w.battery_x, std::max(0.0, w.root_chord - batt_len));
    double x_batt   = batt_x0 + 0.5 * batt_len;          // box CG
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

    // ---- bending-aware spar mass (plan Tier 3) --------------------------------
    // Chord-proportional ultimate load -> cumulative shear/moment (tip->root,
    // projected y for moment arms) -> spar cap area M/(sigma_allow*h) -> mass via
    // arc length s.ds. Appended to `strips` so the Izz/Ixx loop below picks it up
    // for free. spar_enable=0 reproduces pre-Tier-3 masses exactly (regression).
    if (cfg.getb("spar_enable", true) && w.stations.size() >= 2 && mp.mass > 0.0) {
        double n_ult       = cfg.getd("n_ult", 6.0);
        double sigma_allow = cfg.getd("spar_sigma_allow", 150e6);  // Pa; ponytail: calibration
                                                                    // knob, not literal material
                                                                    // strength -- tune so spar
                                                                    // mass lands ~10-15% of
                                                                    // structural mass.
        double rho_spar    = cfg.getd("spar_density", 1600.0);     // kg/m^3
        double spar_frac   = cfg.getd("spar_root_frac", 0.15);     // chordwise station (shared
                                                                    // with the clearance check)
        double W  = mp.mass * GRAV;      // pre-spar weight estimate; spar is a small fraction
        double Sr = mp.S_ref;

        int n = (int)w.stations.size();
        std::vector<double> V(n, 0.0), M(n, 0.0);   // shear, moment; tip (n-1) -> root (0)
        for (int i = n - 2; i >= 0; --i) {
            double y0 = w.stations[i].y,     y1 = w.stations[i + 1].y;
            double Lp0 = n_ult * W * w.stations[i].chord     / Sr;
            double Lp1 = n_ult * W * w.stations[i + 1].chord / Sr;
            double dy = y1 - y0;
            V[i] = V[i + 1] + 0.5 * (Lp0 + Lp1) * dy;
            M[i] = M[i + 1] + 0.5 * (V[i] + V[i + 1]) * dy;
        }

        double spar_m_half = 0.0, spar_mx_half = 0.0;
        for (int i = 0; i < n; ++i) {
            const auto& s = w.stations[i];
            double h = std::max((geom::cst_upper(s.af, spar_frac)
                                 - geom::cst_lower(s.af, spar_frac)) * s.chord, 1e-4);
            double A = M[i] / (sigma_allow * h);       // spar cap area, m^2
            double dm = A * rho_spar * s.ds;            // material follows arc length
            double x_cen = s.x_le + spar_frac * s.chord;
            spar_m_half  += dm;
            spar_mx_half += dm * x_cen;
            strips.push_back({dm, x_cen, s.y, s.z});
        }
        double spar_m_full  = 2.0 * spar_m_half;
        double spar_mx_full = 2.0 * spar_mx_half;
        double mx_total = mp.mass * mp.x_cg + spar_mx_full;
        mp.mass += spar_m_full;
        mp.x_cg = mx_total / mp.mass;
        mp.spar_mass = spar_m_full;
    }

    // ---- yaw/roll inertia Izz, Ixx about CG -----------------------------------
    // Structural strips: each half-station at ±y contributes factor 2.
    // Ixx = Σ dm*(y²+z²); for small dihedral z≈0 so Ixx ≈ Σ dm*y² (spanwise only).
    double Izz = 0.0, Ixx = 0.0;
    for (const auto& sm : strips) {
        double dx = sm.x_cen - mp.x_cg;
        Izz += 2.0 * sm.dm * (dx * dx + sm.y * sm.y);
        Ixx += 2.0 * sm.dm * (sm.y * sm.y + sm.z * sm.z);
    }
    // Point masses (motor + battery at CL; servos at ±y_srv).
    auto addIzz = [&](double m, double x, double y_pm) {
        double dx = x - mp.x_cg;
        Izz += m * (dx * dx + y_pm * y_pm);
        Ixx += m * y_pm * y_pm;
    };
    addIzz(m_motor, x_motor, 0.0);                        // centerline
    addIzz(m_batt,  x_batt,  0.0);                        // centerline
    addIzz(m_avi,   x_avi,   0.0);                        // ponytail: avi at y=0; lift if avi goes off-axis
    addIzz(m_srv, x_srv,  t_srv * w.semi_span);           // right servo
    addIzz(m_srv, x_srv, -t_srv * w.semi_span);           // left servo
    mp.Izz = Izz;
    mp.Ixx = Ixx;

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

    // Motor sits in a local blunt TE boss at the root (CAD plane-split; see
    // config note), NOT inside the sharp-TE airfoil thickness — so there is no
    // airfoil-thickness gate for it. hw_clearance is the avionics check only.

    // Avionics block at t=0.35, 30% chord: use lofted station for accurate section shape.
    double avi_hh = cfg.getd("avionics_half_h", 0.012);
    const Station* avi_st = w.stations.empty() ? nullptr : &w.stations[0];
    {
        double best_dt = 1e9;
        for (const auto& s : w.stations) {
            double dt = std::fabs(s.y / (w.semi_span > 0 ? w.semi_span : 1.0) - 0.35);
            if (dt < best_dt) { best_dt = dt; avi_st = &s; }
        }
    }
    double half_t_avi = avi_st
        ? 0.5 * (geom::cst_upper(avi_st->af, 0.30) - geom::cst_lower(avi_st->af, 0.30)) * avi_st->chord
        : 0.0;
    double avi_clear = half_t_avi - avi_hh;

    mp.hw_clearance = avi_clear;

    // ---- Pusher propeller keep-out ----------------------------------------
    // Disk is FIXED at x = root_chord + hub_gap, centered on thrust axis (y=z=0).
    // Every station inside the disk radius must have its TE forward of the disk face.
    // Previous code anchored to the wing's own aftmost TE → clearance was always ≥0.
    double prop_r    = 0.5 * cfg.getd("prop_diameter", 0.203);
    double hub_gap   = cfg.getd("prop_hub_gap", 0.010);
    double blade_clr = cfg.getd("prop_blade_clear", 0.005);
    double face = w.root_chord + hub_gap - blade_clr;  // fixed disk face (x)
    double prop_cl = 1.0;
    const Station* prop_bind = nullptr;
    for (const auto& s : w.stations) {
        if (std::sqrt(s.y*s.y + s.z*s.z) < prop_r) {
            double cl = face - (s.x_le + s.chord);   // >0 = TE forward of disk
            if (cl < prop_cl) { prop_cl = cl; prop_bind = &s; }
        }
    }
    mp.prop_clearance = prop_cl;
    if (std::getenv("AERO_CV_DIAG") && prop_bind)
        std::cerr << "[prop] clearance=" << prop_cl << " m  binding y="
                  << prop_bind->y << " z=" << prop_bind->z << "\n";

    // Battery box geometry for reporting/CAD (centerline placement, z=0).
    mp.batt_cx = x_batt; mp.batt_cy = 0.0; mp.batt_cz = 0.0;
    mp.batt_lx = batt_len; mp.batt_ly = batt_wid; mp.batt_lz = batt_hgt;

    return mp;
}

}  // namespace massprops
}  // namespace aero
