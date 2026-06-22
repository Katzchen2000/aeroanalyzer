// massprops.h — Phase 2: cross-sectional mass integration (NOT B-Rep).
//
// Per the review: the plan's "inline B-Rep" is really cross-sectional lofting +
// analytic clearance queries. We integrate 2D structural areas across the 20
// stations (trapezoidal) and do spar/OML clearance as a closed-form distance —
// no boolean geometry kernel, no voxels.
#pragma once
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace massprops {

// Normalized (chord=1) enclosed area of the section: ∫(upper-lower)dx.
double section_area_hat(const Airfoil& f);
// Normalized perimeter of the section (for thin-wall shell area).
double section_perimeter_hat(const Airfoil& f);

MassProps compute(const WingGeometry& w, const Config& cfg);

}  // namespace massprops
}  // namespace aero
