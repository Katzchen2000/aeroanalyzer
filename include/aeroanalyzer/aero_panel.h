// aero_panel.h — Milestone 3: Morino (Dirichlet) constant-strength
// source/doublet panel method. Developed alongside the interim VLM and selected
// via config (aero_model = panel). Exposes the influence kernels for unit
// testing plus the solver entry point that potential::solve() dispatches to.
#pragma once
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include "aeroanalyzer/linalg.h"
#include <array>

namespace aero {
namespace panel {

// A flat (low-order) quadrilateral panel, corners ordered CCW about +normal.
struct Quad {
    std::array<Vec3, 4> c;
};

// Per-unit-strength potential influence of a constant source / doublet panel,
// evaluated at field point P. Conventions (verified by the kernel probe):
//   doublet():  Phi_d / mu     = -Omega / (4*pi)   (Omega = signed solid angle)
//               -> interior-face self influence = -1/2.
//   source():   Phi_s / sigma  = -(1/4pi) * INT_panel (1/r) dS
// These match Katz & Plotkin's low-order constant-strength panel relations.
double doublet_potential(const Quad& q, const Vec3& P);
double source_potential(const Quad& q, const Vec3& P);

// Geometry helpers (exposed for tests).
Vec3 quad_centroid(const Quad& q);
Vec3 quad_normal(const Quad& q);     // unit normal from the diagonals
double quad_area(const Quad& q);

// Mesh diagnostics (exposed for tests).
struct MeshStats {
    int n_panels = 0;
    int n_strips = 0;
    int nc = 0;             // chordwise panels per surface
    double wetted_area = 0; // sum of panel areas (~2x planform for a thin wing)
    double min_outward = 1; // min(panel.n . outward_ref), should be > 0
    double self_doublet = 0;// C_ii at a sample panel, should be ~ -0.5
};
MeshStats mesh_debug(const WingGeometry& w, const Config& cfg);

// Assembly/solve diagnostics (exposed for tests). Solves the Morino system at
// the given alpha (V=1, delta_e=0) and reports the Trefftz lift coefficient and
// the closed-body source flux (sum sigma*A, ~0 confirms watertightness).
struct PanelSolveStats {
    int n_panels = 0;
    double cl = 0;          // Trefftz CL from the wake circulation
    double sigma_flux = 0;  // sum sigma_j * A_j (should be ~0)
    double gamma_max = 0;   // peak spanwise circulation
    // Conditioning diagnostics (populated only when PANEL_DEBUG_COND is set in
    // the environment; the SVD is too costly for the GA hot path). cond_C is the
    // 2-norm condition number of the doublet+Kutta matrix; min_cp_gap is the
    // smallest distance between any two panel collocation points (the suspected
    // near-coincidence of the upper/lower trailing-edge panels at high nc).
    double cond_C = 0;
    double min_cp_gap = 0;
};
PanelSolveStats panel_solve_debug(const WingGeometry& w, const Config& cfg,
                                  double alpha);

// Dump the full-span spanwise loading (y, Gamma) for diagnostics.
struct LoadingDump { std::vector<double> y, gamma; };
LoadingDump panel_loading_debug(const WingGeometry& w, const Config& cfg,
                                double alpha);

// Neutral-point diagnostics: the legacy quarter-chord proxy vs the chordwise
// load-centre integration (see neutral_point_load) at a given alpha.
struct XnpDebug { double proxy = 0, load = 0, cl = 0, mac = 0; };
XnpDebug panel_xnp_debug(const WingGeometry& w, const Config& cfg, double alpha);

// Full operating-point solve (same outputs as the VLM path).
AeroState solve(const WingGeometry& w, const MassProps& mp,
                const viscous::Surrogate& surr, const Config& cfg,
                double alpha, double delta_e);

}  // namespace panel
}  // namespace aero
