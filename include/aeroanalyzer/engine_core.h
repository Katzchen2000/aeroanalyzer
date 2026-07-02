// engine_core.h — shared data structures for the whole pipeline.
//
// Everything the GA optimizes flows through these types. Physics modules
// (geom / massprops / aero / stability) read and write them; the GA only ever
// touches Candidate. Keep this header dependency-light.
#pragma once
#include <array>
#include <vector>
#include <string>
#include "aeroanalyzer/linalg.h"

namespace aero {

// ---- physical constants & fixed operating point -------------------------
constexpr double PI       = 3.14159265358979323846;
constexpr double DEG2RAD  = PI / 180.0;
constexpr double RAD2DEG  = 180.0 / PI;
constexpr double RHO      = 1.225;     // air density, kg/m^3 (sea level ISA)
constexpr double MU       = 1.81e-5;   // dynamic viscosity, Pa.s
constexpr double GRAV     = 9.80665;   // m/s^2
constexpr double V_CRUISE = 15.0;      // m/s — target cruise (plan §3)

// Objective indices into Candidate::objectives (all MINIMIZED).
enum Obj { OBJ_DRAG = 0, OBJ_MASS = 1, OBJ_SM = 2, N_OBJ = 3 };

enum class ControlMode { Elevon, Split };

// ---- geometry -----------------------------------------------------------
// CST (Kulfan) airfoil: y(x) = C(x) * S(x) + x * te, with class function
// C = x^N1 (1-x)^N2 and shape S a Bernstein sum. Independent upper/lower
// weight vectors so reflex (essential for a tailless wing) is representable.
struct Airfoil {
    std::vector<double> wu;        // upper Bernstein weights
    std::vector<double> wl;        // lower Bernstein weights
    double te_thick = 0.002;       // trailing-edge thickness, fraction of chord
    double N1 = 0.5, N2 = 1.0;     // round LE, sharp TE
};

// One spanwise station after lofting.
struct Station {
    double y      = 0.0;   // spanwise coordinate, m (0 = root)
    double chord  = 0.0;   // m
    double x_le   = 0.0;   // leading-edge x, m (sweep offset)
    double z      = 0.0;   // dihedral offset, m
    double twist    = 0.0;   // local incidence (washout), rad
    double width    = 0.0;   // projected spanwise strip width, m (S_ref/MAC/control integrals)
    double ds       = 0.0;   // true arc-length strip width, m (mass/material — NOT cosine-projected)
    double dihedral = 0.0;   // local section normal tilt, rad (0=flat, π/2=vertical)
    double eta      = 0.0;   // normalised span parameter t=y_flat/semi_span
    Airfoil af;            // lofted section shape at this station (root->tip blend)
};

// Decoded design: a physical half-wing (mirror assumed).
// Every spanwise distribution is a smooth Bezier curve over eta in [0,1] —
// organic, no faceting. geom::chord_at/xle_at/twist_at/dihedral_at (geom.h)
// are the single source of truth; every consumer reads through them (never
// re-derive chord/x_le/twist/dihedral from these control-point vectors directly).
struct WingGeometry {
    double semi_span  = 0.6;   // m (half span, tip at +semi_span) — gene
    std::vector<double> chord_cp;   // NCP_CHORD Bezier pts, m (CP0 = root chord)
    std::vector<double> sweep_cp;   // NCP_SWEEP Bezier pts of x_le/semi_span (CP0 = 0)
    std::vector<double> twist_cp;   // NCP_TWIST Bezier pts, rad
    std::vector<double> dih_cp;     // NCP_DIH Bezier pts of dihedral angle, rad (CP0 = 0)
    // derived summaries (set by geom::decode/loft) — for analytical downstream
    // models that want one representative scalar rather than the full curve.
    double root_chord  = 0.25;  // m, = chord_at(0)
    double tip_chord   = 0.13;  // m, = chord_at(1)
    double le_sweep    = 0.0;   // rad, net representative sweep = atan2(xle_at(1), semi_span)
    double washout     = 0.0;   // rad, = twist_at(1) - twist_at(0)
    double z_tip       = 0.0;   // dimensionless, = z(tip)/semi_span
    std::vector<Airfoil> sections;  // K control sections; loft() blends piecewise
    ControlMode mode = ControlMode::Elevon;  // ponytail: fixed; G_MODE gene dropped
    double battery_x     = 0.05;  // battery-box CG x location, m (CG trim handle)
    double cs_chord_frac = 0.25;  // control-surface chord fraction [0.15, 0.35]
    double ail_span_frac = 0.60;  // aileron inboard edge, fraction of semi-span [0.40, 0.80]
    std::vector<Station> stations;  // filled by geom::loft()
};

// ---- mass properties ----------------------------------------------------
struct MassProps {
    double mass     = 0.0;  // kg, full wing (both halves + payload)
    double x_cg     = 0.0;  // m
    double mac      = 0.0;  // mean aerodynamic chord, m
    double x_mac_le = 0.0;  // x of the MAC leading edge, m
    double S_ref    = 0.0;  // reference planform area (full), m^2
    double b_full   = 0.0;  // full span, m
    double AR       = 0.0;  // aspect ratio
    double volume   = 0.0;  // structural volume (shell+infill), m^3
    double spar_clearance  = 1.0;  // min spar-to-OML clearance, m (neg = breach)
    double hw_clearance    = 1.0;  // min hardware-to-OML clearance, m (neg = breach)
    double prop_clearance  = 1.0;  // min wing-TE to prop disk, m (neg = intrusion)
    double Izz            = 0.0;  // yaw moment of inertia about CG, kg*m^2
    double Ixx            = 0.0;  // roll moment of inertia about CG, kg*m^2
    // Battery box geometry (report/CAD only; point-mass model unchanged).
    // ponytail: y/z are report-only; add a fit-check gate when those matter.
    double batt_cx = 0.0, batt_cy = 0.0, batt_cz = 0.0;  // box center, m
    double batt_lx = 0.0, batt_ly = 0.0, batt_lz = 0.0;  // L x W x H, m
};

// ---- aerodynamic operating-point result ---------------------------------
struct AeroState {
    double alpha = 0.0;       // trimmed angle of attack, rad
    double delta_e = 0.0;     // elevator/elevon deflection, rad
    double CL = 0.0, CD = 0.0, CDi = 0.0, CDp = 0.0, CM = 0.0;
    double e = 1.0;           // span efficiency
    double x_np = 0.0;        // neutral point, m
    double x_np_high = 0.0;   // neutral point re-evaluated at high alpha, m
    double static_margin = 0.0;  // (x_np - x_cg)/MAC
    double hinge_moment = 0.0;   // worst-case servo torque, kg-cm (pitch+roll)
    double cl_da = 0.0;          // roll control derivative, per rad
    double cl_p  = 0.0;          // roll damping derivative, per (pb/2V)
    double roll_helix = 0.0;     // steady roll helix pb/2V at max diff. deflection
    double cn_da = 0.0;          // aileron yaw derivative; >0 = adverse yaw
    double cn_beta = 0.0;        // weathercock stability, per rad (>0 = stable)
    double cn_r    = 0.0;        // yaw rate damping, per (rb/2V) (<0 = damped)
    double dutch_roll_zeta  = 0.0; // Dutch-roll damping ratio (<0 = divergent)
    double dutch_roll_omega = 0.0; // Dutch-roll natural frequency, rad/s
    double phugoid_zeta     = 0.0; // phugoid damping ratio (Lanchester)
    double polar_confidence = 1.0; // minimum section-polar confidence seen
    std::vector<double> cl_local;  // per-station local lift coefficient
    bool tip_stall = false;
    bool trimmed = false;     // did the trim solve converge?
    // Banked-turn predictions (n·CL_cruise; zero extra panel solve)
    double cn_beta_turn        = 0.0; // weathercock stability at turn CL
    double dutch_roll_zeta_turn = 0.0; // dutch-roll damping at turn CL
    bool   tip_stall_turn      = false; // tip-stall proxy at load factor
};

// ---- GA candidate -------------------------------------------------------
struct Candidate {
    std::vector<double> genes;            // normalized to [lo,hi] bounds
    std::array<double, N_OBJ> objectives{{0, 0, 0}};
    double cv = 0.0;          // total constraint violation (0 = feasible)
    // NSGA-II bookkeeping
    int rank = 0;
    double crowding = 0.0;
    bool evaluated = false;
    // diagnostics carried for reporting (not used by the GA operators)
    std::string note;
};

// Incumbent snapshot for the console dashboard.
struct EngineDashboard {
    int generation = 0;
    int pop_size = 0;
    int feasible_count = 0;
    int front0_size = 0;
    Candidate best_drag, best_mass, best_sm;
};

}  // namespace aero
