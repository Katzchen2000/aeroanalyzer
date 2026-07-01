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
// Every spanwise distribution (chord, LE-sweep, twist, dihedral) is a smooth
// Bezier curve over eta in [0,1] — organic, no faceting, no discrete winglet.
// The vertical shape is a continuous dihedral-ANGLE curve, arc-integrated in
// loft(); a smoothly rising angle near the tip is an organic raised tip
// (gull / bird-wing), with no crease and no separate winglet sub-model.
// Planform genes (fixed indices 0-26); G_SEC() addresses the CST block at 27+.
// Dropped: G_TE (sharp everywhere), G_MODE (fixed Elevon), power-law chord/sweep
// exponents, cubic gull coeffs, winglet cant/eta/blend/taper, bezier_te/fold toggles.
constexpr int NCP_CHORD = 6;   // chord(eta) Bezier control points (degree 5)
constexpr int NCP_SWEEP = 6;   // x_le/semi_span(eta) Bezier pts; CP0 pinned to 0 (not a gene)
constexpr int NCP_TWIST = 6;   // twist(eta) Bezier control points (degree 5)
constexpr int NCP_DIH   = 7;   // dihedral-angle(eta) Bezier pts; CP0 pinned to 0 (not a gene)

constexpr int G_SEMISPAN = 0;                             // half span, m
constexpr int G_CHORD_CP = G_SEMISPAN + 1;                // NCP_CHORD genes (CP0=root)
constexpr int G_SWEEP_CP = G_CHORD_CP + NCP_CHORD;        // NCP_SWEEP-1 genes (CP1..CPn)
constexpr int G_TWIST_CP = G_SWEEP_CP + (NCP_SWEEP - 1);  // NCP_TWIST genes
constexpr int G_DIH_CP   = G_TWIST_CP + NCP_TWIST;        // NCP_DIH-1 genes (CP1..CPn)
constexpr int G_BATTERY  = G_DIH_CP + (NCP_DIH - 1);      // battery-box CG x, m
constexpr int G_CS_CHORD = G_BATTERY + 1;                 // control-surface chord fraction
constexpr int G_AIL_SPAN = G_CS_CHORD + 1;                // aileron inboard edge, frac semi-span

constexpr int N_PLANFORM    = G_AIL_SPAN + 1;             // 27
constexpr int N_CST_PER_SEC = 8;   // 4 wu + 4 wl per section
constexpr int N_SECTIONS    = 5;
constexpr int N_GENES = N_PLANFORM + N_SECTIONS * N_CST_PER_SEC;  // 67

// Canonical η breakpoints for the 5 control sections (CST loft only).
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

// ---- smooth spanwise evaluators (single source of truth) ----------------
// Every consumer of chord/x_le/twist/dihedral reads these — never re-derive
// the wing shape from raw genes/control points elsewhere.
double chord_at(const WingGeometry& w, double eta);      // m
double xle_at(const WingGeometry& w, double eta);         // m
double twist_at(const WingGeometry& w, double eta);       // rad
double dihedral_at(const WingGeometry& w, double eta);    // rad (local curve angle)

// Test/scratch convenience: build control-point curves that exactly reproduce
// a straight-line planform (root->tip chord, 0->tip sweep, root->tip twist,
// flat dihedral) via Bezier's linear-precision property (equally-spaced,
// arithmetic-progression control points reduce to the straight line exactly).
// Production genomes go through decode(); this is for hand-built WingGeometry
// in tests. Assumes w.semi_span is already set. Also fills the derived summary
// fields (root_chord/tip_chord/le_sweep/washout) for consistency.
void set_linear_planform(WingGeometry& w, double root_chord, double tip_chord,
                         double le_sweep_rad, double washout_rad);

// ---- genome -> physical wing -------------------------------------------
GenomeSpec default_genome(const Config& cfg = {});
WingGeometry decode(const std::vector<double>& genes, const GenomeSpec& spec, const Config& cfg = {});
// Cosine-spaced spanwise discretization (plan §3).
void loft(WingGeometry& w, int n_stations);
// Planform integrals from the lofted stations.
void planform(const WingGeometry& w, double& S_ref, double& mac,
              double& x_mac_le, double& b_full, double& AR);

}  // namespace geom
}  // namespace aero
