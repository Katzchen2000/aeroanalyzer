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

// Writes <stem>.avl and <stem>.dat. Returns false on file error.
bool write_case(const std::string& stem, const WingGeometry& w,
                const MassProps& mp, const Config& cfg);

}  // namespace avl
}  // namespace aero
