// aero_potential.h — Phase 3: 3D lifting aerodynamics.
//
// IMPORTANT (Milestone status): solve() currently uses an analytic finite-wing
// REFERENCE MODEL (Prandtl lift slope + strip viscous coupling). It is the
// validation oracle and lets the whole GA pipeline run today. The Morino
// Dirichlet panel solver (plan §5) is Milestone 3 and MUST replace the body of
// solve() behind this exact signature — nothing downstream changes. Until it
// reproduces the analytic gates in tests/test_aero.cpp (lift slope, elliptic e),
// do not trust panel output.
#pragma once
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace potential {

// Prandtl 3D lift-curve slope from 2D slope a0, aspect ratio AR, span eff. e.
// With a0=2π and e=1 this reduces to 2π·AR/(AR+2) — the test gate.
double lift_curve_slope_3d(double a0, double AR, double e);

// Span efficiency estimate from planform taper (reference model).
double oswald_e(double AR, double taper);

// Wing aerodynamic center (≈ neutral point for the inviscid model):
// area-weighted quarter-chord, which captures sweep.
double wing_ac_x(const WingGeometry& w);

// Operating-point solve at prescribed alpha and pitch-control deflection.
AeroState solve(const WingGeometry& w, const MassProps& mp,
                const viscous::Surrogate& surr, const Config& cfg,
                double alpha, double delta_e);

}  // namespace potential
}  // namespace aero
