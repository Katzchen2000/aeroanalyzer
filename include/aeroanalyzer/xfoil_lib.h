// xfoil_lib.h - C++ binding for the in-process XFOIL 6.99 Fortran core.
//
// This is the "accurate" surrogate engine: it links libxfoil.a (real XFOIL,
// compiled from tools/bin/Xfoil699src via build_mingw.ps1) and drives it
// through the iso_c_binding wrapper in tools/xfoil_lib/xfwrap.f. It gives
// xfoil.exe-EXACT viscous polars in-process with ZERO disk I/O, and it
// converges the laminar-separation-bubble cases that the native aero::xfoil
// weak-coupling solver can only approximate.
//
// IMPORTANT constraints baked into this wrapper:
//   * XFOIL keeps all state in Fortran COMMON blocks -> it is NOT thread-safe
//     and NOT reentrant. There is exactly one global solver per process. Never
//     call from more than one thread. To parallelize, fan out across PROCESSES
//     (see the self-spawning worker pool in tools/build_surrogate.cpp).
//   * xfl_init() must run exactly once per process. The Session ctor does it.
//   * Changing Re / Ncrit / Mach requires a per-condition re-init (a stale warm
//     boundary layer from a different condition diverges to NaN). set_condition()
//     tracks the active condition and re-inits only when it actually changes.
//   * XFOIL writes verbose iteration history to Fortran unit 6 (= stdout). Left
//     unsilenced this dominates runtime. Every solve here runs inside a
//     StdoutSilencer that redirects fd 1 to the null device for the duration.
//
// Only the offline surrogate generator links this. The runtime optimizer
// (aeroanalyzer.exe) never depends on XFOIL.
#pragma once
#include <vector>
#include <cmath>

#include "aeroanalyzer/engine_core.h"   // Airfoil
#include "aeroanalyzer/geom.h"          // cst_upper / cst_lower

#if defined(_WIN32)
#include <io.h>
#include <fcntl.h>
#define AERO_NULLDEV "NUL"
#define AERO_DUP _dup
#define AERO_DUP2 _dup2
#define AERO_OPEN _open
#define AERO_CLOSE _close
#define AERO_WRONLY (_O_WRONLY)
#else
#include <unistd.h>
#include <fcntl.h>
#define AERO_NULLDEV "/dev/null"
#define AERO_DUP dup
#define AERO_DUP2 dup2
#define AERO_OPEN open
#define AERO_CLOSE close
#define AERO_WRONLY (O_WRONLY)
#endif

namespace aero {
namespace xfoil_lib {

// ---- the Fortran entry points (tools/xfoil_lib/xfwrap.f) -------------------
extern "C" {
void xfl_init();
void xfl_set_airfoil(double* x, double* y, int n);
void xfl_set_cond(double re, double ncrit, double mach, double xtrt, double xtrb);
void xfl_solve(double adeg, double* cl, double* cd, double* cm, int* conv);
}

// One viscous operating point (2-D section coefficients).
struct Point {
    double alpha = 0.0;   // deg
    double cl = 0.0, cd = 0.0, cm = 0.0;
    bool converged = false;
};

// RAII guard: silence Fortran unit 6 (= fd 1, stdout) for its lifetime, then
// restore it. XFOIL's per-iteration dump is enormous; without this a single
// polar costs seconds of console I/O. fd 2 (stderr) is left alone for our own
// diagnostics. Best-effort: if the dup/redirect fails we simply do not silence.
class StdoutSilencer {
public:
    StdoutSilencer() {
        std::fflush(stdout);
        saved_ = AERO_DUP(1);
        if (saved_ < 0) return;
        nul_ = AERO_OPEN(AERO_NULLDEV, AERO_WRONLY);
        if (nul_ < 0) { AERO_CLOSE(saved_); saved_ = -1; return; }
        AERO_DUP2(nul_, 1);
    }
    ~StdoutSilencer() {
        if (saved_ < 0) return;
        std::fflush(stdout);
        AERO_DUP2(saved_, 1);
        AERO_CLOSE(saved_);
        if (nul_ >= 0) AERO_CLOSE(nul_);
    }
    StdoutSilencer(const StdoutSilencer&) = delete;
    StdoutSilencer& operator=(const StdoutSilencer&) = delete;
private:
    int saved_ = -1;
    int nul_ = -1;
};

// A single global XFOIL solver (one per process). Construct one, then drive it
// airfoil-by-airfoil. The ctor performs the one-time xfl_init().
class Session {
public:
    Session() {
        StdoutSilencer q;
        xfl_init();
        cur_re_ = cur_ncrit_ = cur_mach_ = -1.0;
        have_airfoil_ = false;
    }

    // Panel a CST airfoil into XFOIL's buffer. Coordinates run TE(upper) ->
    // LE -> TE(lower), the standard XFOIL/Selig order. A geometry change
    // invalidates everything; the next solve must follow a set_condition().
    // npts is points per surface (~80 matches xfoil.exe's PANE default of 160).
    void set_airfoil(const Airfoil& f, int npts = 80) {
        x_.clear();
        y_.clear();
        x_.reserve(2 * npts + 1);
        y_.reserve(2 * npts + 1);
        const double PI = 3.14159265358979323846;
        // upper surface, TE (x=1) -> LE (x=0)
        for (int i = npts; i >= 1; --i) {
            double x = 0.5 * (1.0 - std::cos(PI * i / npts));
            x_.push_back(x);
            y_.push_back(geom::cst_upper(f, x));
        }
        // leading edge
        x_.push_back(0.0);
        y_.push_back(0.0);
        // lower surface, LE (x=0) -> TE (x=1)
        for (int i = 1; i <= npts; ++i) {
            double x = 0.5 * (1.0 - std::cos(PI * i / npts));
            x_.push_back(x);
            y_.push_back(geom::cst_lower(f, x));
        }
        StdoutSilencer q;
        xfl_set_airfoil(x_.data(), y_.data(), static_cast<int>(x_.size()));
        have_airfoil_ = true;
        // a fresh geometry drops the active condition (xfl_set_airfoil clears
        // the cached BL); force the next set_condition to actually re-init.
        cur_re_ = cur_ncrit_ = cur_mach_ = -1.0;
    }

    // Establish the viscous operating condition. Re-inits the boundary layer
    // ONLY when the condition changed (or after a geometry change), so an
    // alpha sweep at fixed (Re, Ncrit) keeps its warm start.
    void set_condition(double Re, double Ncrit, double Mach = 0.0,
                       double xtr_top = 1.0, double xtr_bot = 1.0) {
        if (Re == cur_re_ && Ncrit == cur_ncrit_ && Mach == cur_mach_) return;
        StdoutSilencer q;
        xfl_set_cond(Re, Ncrit, Mach, xtr_top, xtr_bot);
        cur_re_ = Re;
        cur_ncrit_ = Ncrit;
        cur_mach_ = Mach;
    }

    // Solve one alpha (degrees) at the active geometry+condition. conv=false
    // (and finite-checked) on any non-convergence; never throws.
    Point solve(double alpha_deg) {
        Point p;
        p.alpha = alpha_deg;
        double cl = 0.0, cd = 0.0, cm = 0.0;
        int conv = 0;
        {
            StdoutSilencer q;
            xfl_solve(alpha_deg, &cl, &cd, &cm, &conv);
        }
        p.cl = cl;
        p.cd = cd;
        p.cm = cm;
        p.converged = (conv != 0) && std::isfinite(cl) && std::isfinite(cd) &&
                      std::isfinite(cm) && cd > 0.0 && cd < 1.0;
        return p;
    }

    // Convenience: a full alpha sweep at one (Re, Ncrit). Re-inits the
    // condition, then marches alpha ascending so each converged point warm-
    // starts the next (XFOIL retains its BL between solves). Non-converged
    // points are returned with converged=false for the caller to drop.
    std::vector<Point> polar(const Airfoil& f, double Re, double Ncrit,
                             double a0_deg, double a1_deg, double da_deg,
                             double Mach = 0.0) {
        set_airfoil(f);
        set_condition(Re, Ncrit, Mach);
        std::vector<Point> out;
        for (double a = a0_deg; a <= a1_deg + 1e-9; a += da_deg)
            out.push_back(solve(a));
        return out;
    }

private:
    std::vector<double> x_, y_;
    double cur_re_, cur_ncrit_, cur_mach_;
    bool have_airfoil_;
};

}  // namespace xfoil_lib
}  // namespace aero
