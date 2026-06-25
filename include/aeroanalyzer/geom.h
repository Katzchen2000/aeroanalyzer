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
// Planform genes (fixed indices 0-11):
constexpr int G_ROOT     = 0;   // root chord, m
constexpr int G_TAPER    = 1;   // tip/root chord ratio
constexpr int G_SEMISPAN = 2;   // half span, m
constexpr int G_SWEEP    = 3;   // leading-edge sweep, deg
constexpr int G_WASHOUT  = 4;   // tip twist, deg (negative = washout)
constexpr int G_BATTERY  = 5;   // battery-box CG x, m
constexpr int G_TE       = 6;   // trailing-edge thickness, fraction of chord
constexpr int G_MODE     = 7;   // <0.5 Elevon, else Split
constexpr int G_CS_CHORD = 8;   // control-surface chord fraction [0.15, 0.35]
constexpr int G_AIL_SPAN = 9;   // aileron inboard edge, fraction of semi-span [0.40, 0.80]
constexpr int G_LE_BOW   = 10;  // LE parabolic bow amplitude, m
constexpr int G_TE_BOW   = 11;  // TE chord augmentation bow, m

constexpr int N_PLANFORM    = 12;
constexpr int N_CST_PER_SEC = 8;   // 4 wu + 4 wl per section
constexpr int N_SECTIONS    = 5;
constexpr int N_GENES = N_PLANFORM + N_SECTIONS * N_CST_PER_SEC;  // 52

// Canonical η breakpoints for the 5 control sections.
constexpr double SECTION_ETA[5] = {0.0, 0.5, 0.75, 0.875, 1.0};

// Index of CST weight gene: sec 0..4, wl 0=upper/1=lower, i 0..3
inline int G_SEC(int sec, int wl, int i) {
    return N_PLANFORM + sec * N_CST_PER_SEC + wl * 4 + i;
}

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
