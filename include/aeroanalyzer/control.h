// control.h — Shared flight-control derivatives (M6).
//
// Centralises the pitch/roll/hinge model that was previously duplicated
// verbatim in aero_potential.cpp and aero_panel.cpp.  Both backends call
// control::compute(); the GA, stability, and evaluate layers are unchanged.
#pragma once
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace control {

struct Derivs {
    // Pitch (same fields as the old local PitchControl, renamed for clarity)
    double CLde    = 0.0;   // CL per rad of symmetric deflection
    double Cmde    = 0.0;   // Cm per rad of symmetric deflection
    double Sf      = 0.0;   // pitch-surface area (full wing), m^2
    double cf_chord = 0.0;  // mean pitch-surface chord, m
    double x_cp    = 0.0;   // pitch-surface centre-of-pressure, m

    // Roll (M6: mode-aware, strip-theory)
    double Cl_da      = 0.0;   // roll control deriv, per rad (pos)
    double Cl_p       = 0.0;   // roll damping, per (pb/2V) (neg)
    double roll_helix = 0.0;   // steady pb/2V at max 4:1 differential deflection

    // Hinge moment (M6: Glauert thin-airfoil, worst-case pitch+roll combined)
    double hinge_moment = 0.0; // kg-cm
};

// Compute all control derivatives for one operating point.
// a      = 3-D lift-curve slope (rad^-1), from potential::lift_curve_slope_3d.
// alpha  = angle of attack (rad) — for hinge-moment scaling.
// delta_e = trimmed pitch deflection (rad) — for worst-case hinge moment.
// cfg    = baseline config (reads aileron_deflect_max_deg, aileron_diff_ratio,
//          motor_diameter is NOT read here — that lives in massprops).
Derivs compute(const WingGeometry& w, const MassProps& mp,
               double a, double alpha, double delta_e, const Config& cfg);

}  // namespace control
}  // namespace aero
