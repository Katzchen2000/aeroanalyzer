// avl_export.h — emit an AVL .avl deck (+ section .dat) for the regression oracle.
//
// AVL is the agreed cross-check (vortex lattice) for lift slope, induced drag,
// and Cm0/reflex. Run: avl <stem>.avl, then OPER -> 'x'. Compare CLa, e, Cm.
#pragma once
#include <string>
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace avl {

// Writes <stem>.avl and <stem>_s*.dat. Returns false on file error.
bool write_case(const std::string& stem, const WingGeometry& w,
                const MassProps& mp, const Config& cfg);

// Writes <stem>_3d.csv directly from lofted stations (no AVL intermediary).
// Produces one contour per station × 2 sides; usable in Fusion360/SolidWorks.
// Opens a blunt TE boss near the root (motor_boss_diameter/_span_frac in cfg),
// tapering to the sharp TE outboard.
bool write_3d_csv(const std::string& stem, const WingGeometry& w, const Config& cfg);

// Writes <stem>.stl: ASCII STL mesh from lofted stations (watertight half-wings
// + root/tip caps) plus a prop-disk triangle fan as a visual marker.
// prop_diameter / prop_hub_gap read from cfg.
bool write_stl(const std::string& stem, const WingGeometry& w, const Config& cfg);

}  // namespace avl
}  // namespace aero
