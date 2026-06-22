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
    double twist  = 0.0;   // local incidence (washout), rad
    double width  = 0.0;   // spanwise strip width for integration, m
};

// Decoded design: a physical half-wing (mirror assumed).
struct WingGeometry {
    double semi_span  = 0.6;   // m (half span, tip at +semi_span)
    double root_chord = 0.25;  // m
    double tip_chord  = 0.13;  // m
    double le_sweep   = 0.0;   // rad (leading-edge sweep)
    double washout    = 0.0;   // rad (tip twist, negative = washout)
    Airfoil section;           // constant section shape (scaled per station)
    ControlMode mode = ControlMode::Elevon;
    double battery_x  = 0.05;  // battery-box CG x location, m (CG trim handle)
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
    double spar_clearance = 1.0;  // min spar-to-OML clearance, m (neg = breach)
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
    double hinge_moment = 0.0;   // required servo torque, kg-cm
    std::vector<double> cl_local;  // per-station local lift coefficient
    bool tip_stall = false;
    bool trimmed = false;     // did the trim solve converge?
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
