// airfoil_io.h — load .dat airfoils (Selig/Lednicer), generate NACA 4-digit,
// and fit them to CST (reusing geom::fit_cst). Feeds both the GA seeding
// (ancestry injection, plan §3) and the offline XFOIL surrogate sampler.
#pragma once
#include <string>
#include <vector>
#include <utility>
#include "aeroanalyzer/engine_core.h"

namespace aero {
namespace airfoil_io {

// Surface ordinates as (x, z) with x ascending in [0,1].
struct Coords {
    std::vector<std::pair<double, double>> upper;
    std::vector<std::pair<double, double>> lower;
    std::string name;
};

// Parse a .dat file. Auto-detects Selig (single TE->LE->TE loop) vs Lednicer
// (count header + upper block + lower block). Sets ok=false on failure.
Coords load_dat(const std::string& path, bool& ok);

// NACA 4-digit (e.g. "2412"), cosine-sampled with n points per surface.
Coords naca4(const std::string& code, int n = 80);

// Trailing-edge thickness (fraction of chord) implied by the coordinates.
double estimate_te(const Coords& c);

// Fit CST weights (order -> order+1 weights per surface) to the coordinates.
Airfoil to_airfoil(const Coords& c, int order, double te);

}  // namespace airfoil_io
}  // namespace aero
