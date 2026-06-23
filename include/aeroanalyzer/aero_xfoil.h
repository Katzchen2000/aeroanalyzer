// aero_xfoil.h - Milestone 4: native in-process viscous airfoil solver.
//
// A clean, modern C++17 reimplementation of the core XFOIL algorithm
// (M. Drela): a Hess-Smith source+vortex inviscid panel method strongly
// coupled to a two-equation integral boundary layer with envelope-e^n
// transition. The boundary-layer closure relations and the e^n amplification
// correlation are ported faithfully from Drela's Fortran (xblsys.f); the
// linear algebra is Eigen, and the legacy Fortran Gaussian elimination is gone.
//
// Design goals (per project requirements):
//   * Zero disk I/O. Airfoil ordinates come in via a CST Airfoil (or a raw
//     coordinate vector); polars are returned in a struct. Nothing is written.
//   * No heap churn in the hot path: a Solver owns all buffers and is reused
//     across an alpha/Re sweep, so a GA can call it from many OpenMP threads
//     (one Solver per thread) without allocator contention.
//   * Warm-starting: a converged boundary-layer state is retained and seeds the
//     next alpha/Re, dropping the Newton/interaction count for sorted sweeps.
//   * Coarse-mesh stable (N ~ 60-160 panels) and fast-failing: hard iteration
//     caps and separation watchdogs make an unphysical GA airfoil return a
//     failure flag immediately rather than spinning, so the caller can hit its
//     synthetic-surrogate fallback instantly.
//
// Incompressible only (Mach ~ 0.04 at the design point): the compressibility
// terms in the closures collapse (M^2 = 0). A hook is left for Karman-Tsien.
#pragma once
#include <vector>
#include "aeroanalyzer/engine_core.h"   // Airfoil

namespace aero {
namespace xfoil {

// One operating-point result. All coefficients are 2-D section values.
struct Result {
    double alpha   = 0.0;   // rad
    double cl      = 0.0;
    double cd      = 0.0;   // total profile drag (far-wake Squire-Young)
    double cdf     = 0.0;   // friction part
    double cdp     = 0.0;   // pressure part (cd - cdf)
    double cm      = 0.0;   // about quarter chord
    double xtr_top = 1.0;   // transition x/c, suction side (1 = none)
    double xtr_bot = 1.0;   // transition x/c, pressure side
    bool   converged = false;
    bool   separated = false;   // turbulent separation detected (near/at clmax)
};

struct Options {
    int    n_panel       = 120;   // panels per surface side (total ~2N)
    double Ncrit         = 9.0;   // e^n transition amplification (config: ncrit)
    int    max_newton    = 30;    // hard cap on global interaction iterations
    int    max_bl_newton = 12;    // hard cap on per-station BL Newton steps
    double tol           = 5e-5;  // interaction convergence (||d ue||/Vinf)
    double relax         = 1.0;   // base under-relaxation for the ue update
    bool   warm_start    = true;  // seed BL state from the previous solve
    // Viscous-inviscid coupling. Weak (default): the BL is marched on the
    // inviscid edge velocity and the section loads are the inviscid loads; this
    // is robust and matches XFOIL drag to a few percent across the unstalled
    // band. Strong: an explicit transpiration interaction also decambers the
    // section (closer cl), but the explicit fixed point is unstable in
    // laminar-bubble regimes and needs the implicit Newton -- experimental.
    bool   strong_coupling = false;
    // Fast-fail watchdogs:
    double hk_sep_limit  = 8.0;   // abort if kinematic shape Hk exceeds this
    double cd_sane_max   = 1.0;   // abort if drag leaves the physical range
};

// Solver owns all working memory and is meant to be reused across a sweep.
// NOT thread-safe: give each OpenMP thread its own Solver instance.
class Solver {
public:
    explicit Solver(const Options& opt = Options());
    ~Solver();

    // Set / change the geometry. Panels the CST airfoil in memory and builds +
    // factors the inviscid influence matrix once (reused across alpha). Returns
    // false if the geometry is degenerate. Clears any warm-start state.
    bool set_airfoil(const Airfoil& f);

    // Set geometry directly from ordered surface coordinates (TE -> upper -> LE
    // -> lower -> TE), e.g. to validate against a .dat. Same contract as above.
    bool set_coords(const std::vector<double>& x, const std::vector<double>& y);

    // Solve one operating point. Re is the chord Reynolds number. Uses the
    // retained BL state as the initial guess when opt.warm_start and the
    // previous solve converged. On non-convergence returns Result{converged=false}.
    Result solve(double alpha_rad, double Re);

    // Forget warm-start state (call when Re jumps or before an independent sweep).
    void reset_state();

    struct Impl;   // public so the .cpp's free helper functions can reach it

private:
    Impl* p_;      // pimpl keeps Eigen out of this header
};

// Convenience: a full alpha sweep with automatic warm-start continuation,
// marching outward from the lowest |alpha| so each step seeds the next. Re and
// Ncrit fixed. Non-converged points are returned with converged=false (the
// caller decides whether to drop them or fall back). Allocates one Solver.
std::vector<Result> sweep(const Airfoil& f, double a0_deg, double a1_deg,
                          double da_deg, double Re, const Options& opt = Options());

}  // namespace xfoil
}  // namespace aero
