// aero_viscous.h — Phase 4: shape-parameterized viscous surrogate bridge.
//
// The GA morphs airfoils, so a single fixed polar won't do. Offline, XFOIL
// (via tools/build_surrogate) calibrates a compact polar — cd0, k, cl_max,
// cl_min, cm0 (cd = cd0 + k·cl²) — for many sampled CST shapes × Reynolds
// numbers, written to data/surrogates/polar_coeffs.csv. At runtime we
// interpolate those coefficients across (shape, Re) and evaluate. If the table
// is absent we fall back to an analytic polar so the pipeline always runs.
//
// Validity = the convex region around the seeds (the GA's CST bounds are derived
// from the seeds), so morphing stays inside the trained hull; queries outside
// are flagged (clamped), never silently extrapolated.
#pragma once
#include <string>
#include <vector>
#include "aeroanalyzer/config.h"

namespace aero {
namespace viscous {

struct Polar {
    double cd = 0.02;
    double cm = 0.0;
    double cl_max = 1.15;
    double cl_min = -0.6;
    bool clamped = false;   // query outside the valid cl/Re/shape hull
};

class Surrogate {
public:
    // Loads data/surrogates/polar_coeffs.csv if present; else analytic fallback.
    // Returns true if real tables were loaded.
    bool load(const std::string& dir, const Config& cfg);

    // Query at section shape (CST weights), lift coefficient cl, chord Re.
    Polar query(const std::vector<double>& shape, double cl, double Re) const;

    bool using_tables() const { return tables_loaded_; }
    std::size_t sample_count() const { return samples_.size(); }

private:
    struct Sample {
        std::vector<double> shape;
        double Re = 0.0;
        double cd0 = 0.0, k = 0.0, cl_max = 0.0, cl_min = 0.0, cm0 = 0.0;
    };

    Polar analytic(double cl, double Re) const;   // table-free fallback

    bool tables_loaded_ = false;
    double ncrit_ = 4.0;
    std::vector<Sample> samples_;
    std::vector<double> hull_lo_, hull_hi_;   // per shape dimension
    double re_lo_ = 0.0, re_hi_ = 0.0;
};

}  // namespace viscous
}  // namespace aero
