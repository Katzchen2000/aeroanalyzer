// geom.h — Phase 1: CST geometry, genome decode, lofting, thin-airfoil coeffs.
#pragma once
#include <vector>
#include <string>
#include <utility>
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace geom {

// Design-variable layout. Order matters — keep in sync with decode().
enum Gene {
    G_ROOT = 0,   // root chord, m
    G_TAPER,      // tip/root chord ratio
    G_SEMISPAN,   // half span, m
    G_SWEEP,      // leading-edge sweep, deg
    G_WASHOUT,    // tip twist, deg (negative = washout)
    G_BATTERY,    // battery-box CG x, m
    G_WU0, G_WU1, G_WU2, G_WU3,   // upper CST weights
    G_WL0, G_WL1, G_WL2, G_WL3,   // lower CST weights
    G_TE,         // trailing-edge thickness, fraction of chord
    G_MODE,       // <0.5 Elevon, else Split
    G_CS_CHORD,   // control-surface chord fraction [0.15, 0.35]
    G_AIL_SPAN,   // aileron inboard edge, fraction of semi-span [0.40, 0.80]
    N_GENES
};

struct GenomeSpec {
    std::vector<double> lo, hi;
    std::vector<std::string> names;
    std::size_t size() const { return lo.size(); }
};
GenomeSpec default_genome();

// ---- CST primitives -----------------------------------------------------
double bernstein(int i, int n, double x);
// Surface ordinate at chordwise fraction x in [0,1].
double cst_upper(const Airfoil& f, double x);
double cst_lower(const Airfoil& f, double x);
// Camber line z_c(x) = (upper+lower)/2 sampled at n cosine-clustered points.
std::vector<std::pair<double, double>> camber_line(const Airfoil& f, int n);

// Least-squares fit of CST weights to sampled (x,z) ordinates (normal eqns).
// Used by the round-trip validation gate and "ancestry injection" (plan §3).
Airfoil fit_cst(const std::vector<std::pair<double, double>>& upper,
                const std::vector<std::pair<double, double>>& lower,
                int order, double te_thick);

// Thin-airfoil-theory section coefficients derived from the camber line.
// Reflex shows up as a positive cm_ac, which is what lets a tailless wing trim.
struct ThinAirfoil {
    double alpha_L0 = 0.0;  // zero-lift angle, rad
    double cm_ac    = 0.0;  // moment about aerodynamic center (~c/4)
    double cl_alpha = 2.0 * PI;  // per rad
};
ThinAirfoil thin_airfoil(const Airfoil& f);

// ---- genome -> physical wing -------------------------------------------
WingGeometry decode(const std::vector<double>& genes, const GenomeSpec& spec);
// Cosine-spaced spanwise discretization (plan §3).
void loft(WingGeometry& w, int n_stations);
// Planform integrals from the lofted stations.
void planform(const WingGeometry& w, double& S_ref, double& mac,
              double& x_mac_le, double& b_full, double& AR);

}  // namespace geom
}  // namespace aero
