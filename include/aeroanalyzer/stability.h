// stability.h — Phase 5: trim, static-margin objective, control gates.
#pragma once
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace stability {

// Newton–Raphson trim: solve for (alpha, delta_e) such that CL = W/(qS) and
// Cm_cg = 0 (plan §7). Jacobian by finite difference so it stays valid when the
// reference aero is swapped for the nonlinear panel solver. Returns the trimmed
// AeroState; state.trimmed is false if it failed to converge.
AeroState trim(const WingGeometry& w, const MassProps& mp,
               const viscous::Surrogate& surr, const Config& cfg);

// SM objective (MINIMIZED): 0 inside [lo,hi], else distance to the nearest edge.
double sm_objective(double sm, double lo, double hi);

}  // namespace stability
}  // namespace aero
