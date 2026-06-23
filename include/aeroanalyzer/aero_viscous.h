// aero_viscous.h — shape-parameterized viscous surrogate bridge.
//
// The GA morphs airfoils, so a single fixed polar won't do. Every query returns
// a compact polar — cd0, k, cl_max, cl_min, cm0 (cd = cd0 + k·cl²) — for the
// section's CST shape and Reynolds number. Three interchangeable backends supply
// it (config viscous_backend):
//   - NeuralFoil (default): an in-process neural surrogate of XFOIL evaluated on
//     the fly (see aero_neuralfoil.h). Continuous, smooth, GPL-free.
//   - Table: a pre-computed polar CSV (data/surrogates/polar_coeffs.csv)
//     interpolated across (shape, Re) by inverse-distance weighting, with a
//     convex-hull guard (queries outside the trained region are flagged, never
//     silently extrapolated). Read-only here; the offline generator was retired.
//   - Analytic: a closed-form stand-in so the pipeline always runs.
#pragma once
#include <string>
#include <vector>
#include "aeroanalyzer/config.h"
#include "aeroanalyzer/aero_neuralfoil.h"

namespace aero {
namespace viscous {

struct Polar {
    double cd = 0.02;
    double cm = 0.0;
    double cl_max = 1.15;
    double cl_min = -0.6;
    bool clamped = false;   // query outside the valid cl/Re/shape hull
};

// Viscous engine behind the surrogate. NeuralFoil is the default (in-process
// neural surrogate of XFOIL); Table is the legacy offline XFOIL polar CSV with
// IDW interpolation; Analytic is the table-free closed-form fallback.
enum class Backend { NeuralFoil, Table, Analytic };

class Surrogate {
public:
    // Selects the backend from cfg (viscous_backend, default "neuralfoil"),
    // loading NeuralFoil weights or polar_coeffs.csv as needed. Falls back
    // gracefully (NeuralFoil -> Table -> Analytic). Returns true if a data-driven
    // backend (NeuralFoil weights or a real table) was loaded.
    bool load(const std::string& dir, const Config& cfg);

    // Query at section shape (CST weights), lift coefficient cl, chord Re. te is
    // the section trailing-edge thickness fraction (used only by the NeuralFoil
    // backend; the Table/Analytic paths ignore it).
    Polar query(const std::vector<double>& shape, double cl, double Re,
                double te = 0.005) const;

    bool using_tables() const { return tables_loaded_; }
    bool using_neuralfoil() const { return backend_ == Backend::NeuralFoil; }
    std::size_t sample_count() const { return samples_.size(); }

private:
    struct Sample {
        std::vector<double> shape;
        double Re = 0.0;
        double cd0 = 0.0, k = 0.0, cl_max = 0.0, cl_min = 0.0, cm0 = 0.0;
    };

    Polar analytic(double cl, double Re) const;   // table-free fallback

    // NeuralFoil path: synthesize the compact polar (cd0, k, cl_max, cl_min, cm0)
    // for a (shape, te, Re) by an alpha sweep + parabola fit. Sets low_conf when
    // the swept operating band falls out of the network's trust region.
    void nf_coeffs(const std::vector<double>& shape, double te, double Re,
                   double& cd0, double& k, double& cl_max, double& cl_min,
                   double& cm0, bool& low_conf) const;

    Backend backend_ = Backend::NeuralFoil;
    nf::NeuralFoil nf_;
    double xtr_upper_ = 1.0, xtr_lower_ = 1.0;  // forced-transition x/c (free=1)

    bool tables_loaded_ = false;
    double ncrit_ = 4.0;
    std::vector<Sample> samples_;
    std::vector<double> hull_lo_, hull_hi_;   // per shape dimension
    double re_lo_ = 0.0, re_hi_ = 0.0;
};

}  // namespace viscous
}  // namespace aero
