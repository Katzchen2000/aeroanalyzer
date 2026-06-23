// aero_xfoil.cpp - native in-process viscous airfoil solver (see aero_xfoil.h).
//
// Stage map (mirrors the milestone-4 plan):
//   [1] paneling from CST            -- panel_geometry()
//   [2] Hess-Smith inviscid panel    -- build_inviscid(), inviscid_ue()
//   [3] integral-BL closures (Drela) -- the lam_*/turb_*/amp_* helpers
//   [4] BL march + e^n transition    -- march_side()
//   [5] viscous-inviscid coupling    -- solve()
//   [6] far-wake Squire-Young drag   -- solve()
//
// The boundary-layer closures and the e^n amplification correlation are ported
// from Drela's xblsys.f (DAMPL, DIL/CFL/HSL, DIT/CFT/HST, DILW, HKIN) with the
// incompressible (M=0) simplifications folded in. The legacy Fortran Gaussian
// elimination is replaced by Eigen's PartialPivLU, factored once per geometry
// and reused across alpha and across the viscous-inviscid interaction loop.
#include "aeroanalyzer/aero_xfoil.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/engine_core.h"

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdlib>

namespace aero {
namespace xfoil {

namespace {
constexpr double TWO_PI = 6.283185307179586;
constexpr double PI_    = 3.141592653589793;

// Drela BL constants (xbl.f INISET / BLPAR.INC).
constexpr double SCCON = 5.6;
constexpr double GACON = 6.70;
constexpr double GBCON = 0.75;
constexpr double GCCON = 18.0;
constexpr double DLCON = 0.9;
constexpr double DUXCON = 1.0;
constexpr double CTCON = 0.5 / (GACON * GACON * GBCON);   // 0.01485...
}  // namespace

// ===========================================================================
//  Solver state (pimpl)
// ===========================================================================
struct Solver::Impl {
    Options opt;

    // ---- geometry (paneled surface) ----
    int N = 0;                       // number of panels
    std::vector<double> X, Y;        // node coords, size N+1 (TE..upper..LE..lower..TE)
    std::vector<double> xc, yc;      // panel midpoints, size N
    std::vector<double> sint, cost;  // panel orientation, size N
    std::vector<double> len;         // panel lengths, size N
    std::vector<double> sarc;        // surface arc length at midpoints, size N

    // ---- inviscid system (Hess-Smith: N sources + 1 vortex) ----
    Eigen::PartialPivLU<Eigen::MatrixXd> lu;   // factored (N+1)x(N+1)
    Eigen::MatrixXd At;    // tangential source influence At(i,j), size N x N
    Eigen::VectorXd Avt;   // tangential vortex influence summed over panels, size N
    Eigen::MatrixXd Mresp; // edge-velocity response to unit wall blowing, N x N
                           //   dUe_i = sum_j Mresp(i,j) * wn_j  (transpiration)
    // two unit inviscid solutions (alpha = 0 and 90 deg)
    Eigen::VectorXd uea, ueb;   // inviscid ue at midpoints, size N
    bool geom_ok = false;

    // ---- warm-start BL state (per surface node, retained across solves) ----
    bool have_state = false;
    double state_Re = -1.0;
    std::vector<double> st_theta, st_dstar, st_ue;  // size N (indexed by panel)

    explicit Impl(const Options& o) : opt(o) {}
};

// ===========================================================================
//  [1] Paneling from CST
// ===========================================================================
// Fill panel midpoints / orientation / arc length from the node arrays m.X,m.Y.
// Enforces a clockwise node order (airfoil interior on the right) so that the
// outward normal n = (-sint, cost) and the source self-term (+1/2) are
// consistent. Returns false if degenerate.
static bool finalize_geometry(Solver::Impl& m) {
    int nn = static_cast<int>(m.X.size());
    int N = nn - 1;
    if (N < 40) return false;
    // signed area (shoelace): > 0 means counter-clockwise -> reverse to CW.
    double area2 = 0.0;
    for (int i = 0; i < N; ++i)
        area2 += m.X[i] * m.Y[i + 1] - m.X[i + 1] * m.Y[i];
    if (area2 > 0.0) {
        std::reverse(m.X.begin(), m.X.end());
        std::reverse(m.Y.begin(), m.Y.end());
    }
    m.N = N;
    m.xc.assign(N, 0.0); m.yc.assign(N, 0.0);
    m.sint.assign(N, 0.0); m.cost.assign(N, 0.0); m.len.assign(N, 0.0);
    m.sarc.assign(N, 0.0);
    double acc = 0.0;
    for (int i = 0; i < N; ++i) {
        double dx = m.X[i + 1] - m.X[i];
        double dy = m.Y[i + 1] - m.Y[i];
        double L = std::sqrt(dx * dx + dy * dy);
        if (L < 1e-12) return false;
        m.len[i] = L;
        m.cost[i] = dx / L;
        m.sint[i] = dy / L;
        m.xc[i] = 0.5 * (m.X[i] + m.X[i + 1]);
        m.yc[i] = 0.5 * (m.Y[i] + m.Y[i + 1]);
        m.sarc[i] = acc + 0.5 * L;
        acc += L;
    }
    return true;
}

// Build nodes from a CST airfoil: TE -> upper -> LE -> lower -> TE, cosine
// clustered. finalize_geometry() then normalizes orientation.
static bool panel_geometry(Solver::Impl& m, const Airfoil& f, int nside) {
    nside = std::max(20, nside);
    m.X.clear(); m.Y.clear();
    m.X.reserve(2 * nside + 1);
    m.Y.reserve(2 * nside + 1);
    for (int j = nside; j >= 1; --j) {
        double x = 0.5 * (1.0 - std::cos(PI_ * j / nside));
        m.X.push_back(x);
        m.Y.push_back(geom::cst_upper(f, x));
    }
    m.X.push_back(0.0);
    m.Y.push_back(0.0);
    for (int j = 1; j <= nside; ++j) {
        double x = 0.5 * (1.0 - std::cos(PI_ * j / nside));
        m.X.push_back(x);
        m.Y.push_back(geom::cst_lower(f, x));
    }
    return finalize_geometry(m);
}

// ===========================================================================
//  [2] Hess-Smith inviscid panel method
// ===========================================================================
// Influence of a unit constant-strength source AND unit constant-strength
// vortex on panel j, evaluated at point (px,py), returned as global-frame
// velocity components. Uses the source/vortex duality (vortex = source rotated
// 90deg). i==j self term handled by the caller passing is_self.
static inline void panel_influence(const Solver::Impl& m, int j, double px,
                                   double py, bool is_self,
                                   double& su, double& sv,   // source -> (u,v) global
                                   double& vu, double& vv) {  // vortex -> (u,v) global
    double c = m.cost[j], s = m.sint[j], L = m.len[j];
    double dx = px - m.X[j], dy = py - m.Y[j];
    // into panel-local frame
    double xp = dx * c + dy * s;
    double yp = -dx * s + dy * c;
    double ln, beta;
    if (is_self) {
        ln = 0.0;          // r1 == r2 at the midpoint
        beta = PI_;        // control point just outside the surface
    } else {
        double r1 = std::sqrt(xp * xp + yp * yp);
        double r2 = std::sqrt((xp - L) * (xp - L) + yp * yp);
        ln = std::log(r1 / r2);
        double th1 = std::atan2(yp, xp);
        double th2 = std::atan2(yp, xp - L);
        beta = th2 - th1;
    }
    double inv = 1.0 / TWO_PI;
    // local-frame velocities
    double sul = ln * inv,  svl = beta * inv;      // source (u,w)_local
    double vul = -beta * inv, vvl = ln * inv;      // vortex (u,w)_local = (-w_s,u_s)
    // rotate local -> global:  ug = ul*c - wl*s ; vg = ul*s + wl*c
    su = sul * c - svl * s;  sv = sul * s + svl * c;
    vu = vul * c - vvl * s;  vv = vul * s + vvl * c;
}

// Build the (N+1)x(N+1) Hess-Smith system, factor it, and precompute the two
// unit inviscid surface-velocity solutions (alpha = 0 and 90deg).
static bool build_inviscid(Solver::Impl& m) {
    const int N = m.N;
    Eigen::MatrixXd A(N + 1, N + 1);
    A.setZero();
    m.At.resize(N, N);
    m.Avt.resize(N);
    m.Avt.setZero();

    for (int i = 0; i < N; ++i) {
        double nx = -m.sint[i], ny = m.cost[i];   // panel-i outward normal
        double tx = m.cost[i],  ty = m.sint[i];   // tangent
        double vsum = 0.0;        // vortex normal influence summed over j
        for (int j = 0; j < N; ++j) {
            double su, sv, vu, vv;
            panel_influence(m, j, m.xc[i], m.yc[i], i == j, su, sv, vu, vv);
            double s_n = su * nx + sv * ny;   // source normal infl
            double s_t = su * tx + sv * ty;   // source tangential infl
            double v_n = vu * nx + vv * ny;   // vortex normal infl
            double v_t = vu * tx + vv * ty;   // vortex tangential infl
            A(i, j) = s_n;                    // tangency: source contributions
            vsum += v_n;
            m.At(i, j) = s_t;
            m.Avt(i) += v_t;
        }
        A(i, N) = vsum;                       // coefficient of the single gamma
    }
    // Kutta condition: tangential velocity equal & opposite at the two TE panels.
    {
        int p1 = 0, p2 = N - 1;               // first (upper TE) and last (lower TE) panels
        double t1x = m.cost[p1], t1y = m.sint[p1];
        double t2x = m.cost[p2], t2y = m.sint[p2];
        double vsum = 0.0;
        for (int j = 0; j < N; ++j) {
            double su, sv, vu, vv;
            panel_influence(m, j, m.xc[p1], m.yc[p1], p1 == j, su, sv, vu, vv);
            double s_t1 = su * t1x + sv * t1y, v_t1 = vu * t1x + vv * t1y;
            panel_influence(m, j, m.xc[p2], m.yc[p2], p2 == j, su, sv, vu, vv);
            double s_t2 = su * t2x + sv * t2y, v_t2 = vu * t2x + vv * t2y;
            A(N, j) = s_t1 + s_t2;
            vsum += v_t1 + v_t2;
        }
        A(N, N) = vsum;
    }

    m.lu.compute(A);

    // two unit RHS solutions
    auto solve_unit = [&](double vx, double vy, Eigen::VectorXd& ue) {
        Eigen::VectorXd rhs(N + 1);
        for (int i = 0; i < N; ++i) {
            double nx = -m.sint[i], ny = m.cost[i];
            rhs(i) = -(vx * nx + vy * ny);
        }
        int p1 = 0, p2 = N - 1;
        rhs(N) = -((vx * m.cost[p1] + vy * m.sint[p1]) +
                   (vx * m.cost[p2] + vy * m.sint[p2]));
        Eigen::VectorXd sol = m.lu.solve(rhs);
        double gam = sol(N);
        ue.resize(N);
        for (int i = 0; i < N; ++i) {
            double base = vx * m.cost[i] + vy * m.sint[i];   // freestream tangential
            double s = (m.At.row(i) * sol.head(N))(0);
            ue(i) = base + s + gam * m.Avt(i);
        }
        if (std::getenv("XF_DEBUG"))
            std::fprintf(stderr, "[xf] unit(%.0f,%.0f): gamma=%.5f ue[0]=%.4f ue[N-1]=%.4f sigma0=%.4f\n",
                         vx, vy, gam, ue(0), ue(N - 1), sol(0));
    };
    solve_unit(1.0, 0.0, m.uea);
    solve_unit(0.0, 1.0, m.ueb);

    // Edge-velocity response to wall transpiration (the viscous coupling matrix).
    // A unit outward blowing velocity at panel j enters the tangency RHS as e_j;
    // back-substituting through the factored LU gives the source/vortex response,
    // whose tangential trace is column j of Mresp. Built once, reused every alpha
    // and every interaction sweep. (Kutta RHS unaffected by blowing.)
    m.Mresp.resize(N, N);
    Eigen::VectorXd rhs(N + 1), sol(N + 1);
    for (int j = 0; j < N; ++j) {
        rhs.setZero();
        rhs(j) = 1.0;
        sol = m.lu.solve(rhs);
        double gam = sol(N);
        for (int i = 0; i < N; ++i)
            m.Mresp(i, j) = (m.At.row(i) * sol.head(N))(0) + gam * m.Avt(i);
    }
    return true;
}

// Inviscid surface tangential velocity at angle alpha (linear superposition).
static inline void inviscid_ue(const Solver::Impl& m, double alpha,
                               std::vector<double>& ue) {
    double ca = std::cos(alpha), sa = std::sin(alpha);
    ue.resize(m.N);
    for (int i = 0; i < m.N; ++i) ue[i] = ca * m.uea(i) + sa * m.ueb(i);
}

// Integrate surface pressure to get cl, cm about the quarter chord. ue here is
// the (signed) surface tangential velocity at panel midpoints; Cp = 1 - ue^2.
static void forces_from_cp(const Solver::Impl& m, const std::vector<double>& ue,
                           double alpha, double& cl, double& cm) {
    double ca = std::cos(alpha), sa = std::sin(alpha);
    double fx = 0.0, fy = 0.0, mq = 0.0;
    for (int i = 0; i < m.N; ++i) {
        double cp = 1.0 - ue[i] * ue[i];
        double nx = -m.sint[i], ny = m.cost[i];          // outward normal
        // pressure force on panel = -Cp * n * len  (per unit chord, q=1)
        double dfx = -cp * nx * m.len[i];
        double dfy = -cp * ny * m.len[i];
        fx += dfx; fy += dfy;
        // moment about (0.25, 0) of the pressure force
        double rx = m.xc[i] - 0.25, ry = m.yc[i];
        mq += rx * dfy - ry * dfx;
    }
    // rotate body-axis force into wind axes: lift is perpendicular to freestream
    cl = fy * ca - fx * sa;
    cm = mq;
}

// ===========================================================================
//  Public Solver wrappers
// ===========================================================================
Solver::Solver(const Options& opt) : p_(new Impl(opt)) {}
Solver::~Solver() { delete p_; }
void Solver::reset_state() { p_->have_state = false; }

bool Solver::set_airfoil(const Airfoil& f) {
    p_->geom_ok = false;
    p_->have_state = false;
    if (!panel_geometry(*p_, f, p_->opt.n_panel)) return false;
    if (!build_inviscid(*p_)) return false;
    p_->geom_ok = true;
    return true;
}

bool Solver::set_coords(const std::vector<double>& x,
                        const std::vector<double>& y) {
    p_->geom_ok = false;
    p_->have_state = false;
    if (x.size() != y.size() || x.size() < 41) return false;
    Solver::Impl& m = *p_;
    m.X = x; m.Y = y;
    if (!finalize_geometry(m)) return false;
    if (!build_inviscid(m)) return false;
    m.geom_ok = true;
    return true;
}

// Forward declaration of the viscous driver (stage 4-6); defined below.
static Result viscous_solve(Solver::Impl& m, double alpha, double Re);

Result Solver::solve(double alpha_rad, double Re) {
    if (!p_->geom_ok) return Result{};
    return viscous_solve(*p_, alpha_rad, Re);
}

std::vector<Result> sweep(const Airfoil& f, double a0_deg, double a1_deg,
                          double da_deg, double Re, const Options& opt) {
    std::vector<Result> out;
    Solver solver(opt);
    if (!solver.set_airfoil(f)) return out;
    const double d2r = PI_ / 180.0;
    if (da_deg <= 0) da_deg = 0.5;
    // march outward from the smallest |alpha| so warm-start always flows from a
    // well-behaved point into the harder ones.
    std::vector<double> alphas;
    for (double a = a0_deg; a <= a1_deg + 1e-9; a += da_deg) alphas.push_back(a);
    // order by |alpha| ascending
    std::sort(alphas.begin(), alphas.end(),
              [](double p, double q) { return std::fabs(p) < std::fabs(q); });
    std::vector<Result> tmp(alphas.size());
    solver.reset_state();
    for (std::size_t i = 0; i < alphas.size(); ++i)
        tmp[i] = solver.solve(alphas[i] * d2r, Re);
    // return sorted by alpha ascending for the caller's convenience
    std::vector<std::size_t> idx(alphas.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b) { return alphas[a] < alphas[b]; });
    out.reserve(alphas.size());
    for (std::size_t k : idx) out.push_back(tmp[k]);
    return out;
}

// ===========================================================================
//  [3] Integral-BL closure relations (Drela xblsys.f), value-only ports.
//  Incompressible (M=0): Hk == H, density-shape H** == 0, the compressibility
//  fudges in HST/CFT collapse. We carry only the values; the per-station Newton
//  forms its Jacobian numerically, so the analytic _HK/_RT sensitivities in the
//  Fortran are deliberately not ported.
// ===========================================================================
namespace {

// Laminar H* (KE-shape) correlation -- HSL.
double cl_hsl(double hk) {
    if (hk < 4.35) {
        double t = hk - 4.35;
        return 0.0111 * t * t / (hk + 1.0) - 0.0278 * t * t * t / (hk + 1.0)
               + 1.528 - 0.0002 * (t * hk) * (t * hk);
    }
    return 0.015 * (hk - 4.35) * (hk - 4.35) / hk + 1.528;
}

// Turbulent H* correlation -- HST (M=0 so the FM compressibility factor is 1).
double cl_hst(double hk, double rt) {
    const double HSMIN = 1.5, DHSINF = 0.015;
    double ho = (rt > 400.0) ? 3.0 + 400.0 / rt : 4.0;
    double rtz = (rt > 200.0) ? rt : 200.0;
    if (hk < ho) {
        double hr = (ho - hk) / (ho - 1.0);
        return (2.0 - HSMIN - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5)
               + HSMIN + 4.0 / rtz;
    }
    double grt = std::log(rtz);
    double hdif = hk - ho;
    double rtmp = hk - ho + 4.0 / grt;
    double htmp = 0.007 * grt / (rtmp * rtmp) + DHSINF / hk;
    return hdif * hdif * htmp + HSMIN + 4.0 / rtz;
}

// Laminar skin friction -- CFL  ( returns Cf ).
double cl_cfl(double hk, double rt) {
    if (hk < 5.5) {
        double t = (5.5 - hk) * (5.5 - hk) * (5.5 - hk) / (hk + 1.0);
        return (0.0727 * t - 0.07) / rt;
    }
    double t = 1.0 - 1.0 / (hk - 4.5);
    return (0.015 * t * t - 0.07) / rt;
}

// Turbulent skin friction -- CFT  ( returns Cf, M=0 so FC=1 ).
double cl_cft(double hk, double rt) {
    double grt = std::log(rt);
    if (grt < 3.0) grt = 3.0;
    double gex = -1.74 - 0.31 * hk;
    double arg = -1.33 * hk;
    if (arg < -20.0) arg = -20.0;
    double thk = std::tanh(4.0 - hk / 0.875);
    double cfo = 0.3 * std::exp(arg) * std::pow(grt / 2.3026, gex);
    return cfo + 1.1e-4 * (thk - 1.0);
}

// Laminar dissipation 2CD/H* -- DIL.
double cl_dil(double hk, double rt) {
    if (hk < 4.0)
        return (0.00205 * std::pow(4.0 - hk, 5.5) + 0.207) / rt;
    double hkb = hk - 4.0;
    double den = 1.0 + 0.02 * hkb * hkb;
    return (-0.0016 * hkb * hkb / den + 0.207) / rt;
}

// Laminar wake dissipation 2CD/H* -- DILW.
double cl_dilw(double hk, double rt) {
    double hs = cl_hsl(hk);
    double rcd = 1.10 * (1.0 - 1.0 / hk) * (1.0 - 1.0 / hk) / hk;
    return 2.0 * rcd / (hs * rt);
}

// Envelope e^n amplification rate dN/dxi -- DAMPL.
double cl_dampl(double hk, double th, double rt) {
    const double DGR = 0.08;
    double hmi = 1.0 / (hk - 1.0);
    double aa = 2.492 * std::pow(hmi, 0.43);
    double bb = std::tanh(14.0 * hmi - 9.24);
    double grcrit = aa + 0.7 * (bb + 1.0);
    double gr = std::log10(rt);
    if (gr < grcrit - DGR) return 0.0;
    double rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    double rfac = (rnorm >= 1.0) ? 1.0 : 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
    double arg = 3.87 * hmi - 2.52;
    double ex = std::exp(-arg * arg);
    double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    double af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
    return af * dadr / th * rfac;
}

// Average amplification over an interval -- AXSET (rms blend + N->Ncrit nudge).
double cl_axset(double hk1, double t1, double rt1, double a1,
                double hk2, double t2, double rt2, double a2, double acrit) {
    double ax1 = cl_dampl(hk1, t1, rt1);
    double ax2 = cl_dampl(hk2, t2, rt2);
    double axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    double axa = (axsq <= 0.0) ? 0.0 : std::sqrt(axsq);
    double arg = 20.0 * (acrit - 0.5 * (a1 + a2));
    if (arg > 20.0) arg = 20.0;
    double exn = (arg <= 0.0) ? 1.0 : std::exp(-arg);
    double dax = exn * 0.002 / (t1 + t2);
    return axa + dax;
}

// Secondary "2" variables from the primary state (Ue, theta, dstar, s) for the
// given station type (1 lam, 2 turb, 3 wake). At M=0, Hk == H. Re is the chord
// Reynolds number so Rtheta = Ue*theta*Re.  `s` is the amplification N for a
// laminar station (unused here) or the sqrt(Ctau) lag variable when turbulent.
struct Aux {
    double h, hk, rt, hs, us, cq, cf, di, de;
};
Aux bl_aux(double u, double t, double d, double s, int ityp, double Re) {
    Aux a;
    a.h = d / t;
    double hk = a.h;
    hk = (ityp == 3) ? std::max(hk, 1.00005) : std::max(hk, 1.05);
    a.hk = hk;
    a.rt = u * t * Re;
    if (a.rt < 1.0) a.rt = 1.0;
    a.hs = (ityp == 1) ? cl_hsl(hk) : cl_hst(hk, a.rt);

    a.us = 0.5 * a.hs * (1.0 - (hk - 1.0) / (GBCON * a.h));
    if (ityp <= 2 && a.us > 0.95) a.us = 0.98;
    if (ityp == 3 && a.us > 0.99995) a.us = 0.99995;

    double gcc = (ityp == 2) ? GCCON : 0.0;
    double hkc = hk - 1.0 - gcc / a.rt;
    if (hkc < 0.01) hkc = 0.01;
    double hkb = hk - 1.0, usb = 1.0 - a.us;
    a.cq = std::sqrt(CTCON * a.hs * hkb * hkc * hkc / (usb * a.h * hk * hk));

    if (ityp == 3) {
        a.cf = 0.0;
    } else if (ityp == 1) {
        a.cf = cl_cfl(hk, a.rt);
    } else {
        a.cf = cl_cft(hk, a.rt);
        double cfl = cl_cfl(hk, a.rt);
        if (cfl > a.cf) a.cf = cfl;          // laminar floor at tiny Rtheta
    }

    if (ityp == 1) {
        a.di = cl_dil(hk, a.rt);
    } else {
        double di;
        if (ityp == 2) {
            double cft = cl_cft(hk, a.rt);
            di = (0.5 * cft * a.us) * 2.0 / a.hs;     // wall contribution
            double grt = std::log(a.rt);
            double hmin = 1.0 + 2.1 / grt;
            double fl = (hk - 1.0) / (hmin - 1.0);
            di *= 0.5 + 0.5 * std::tanh(fl);          // low-Hk wall fudge (DFAC)
        } else {
            di = 0.0;                                  // wake: no wall shear
        }
        // turbulent outer-layer dissipation (S = sqrt(Ctau))
        di += s * s * (0.995 - a.us) * 2.0 / a.hs;
        di += 0.15 * (0.995 - a.us) * (0.995 - a.us) / a.rt * 2.0 / a.hs;
        if (ityp == 2) {
            double dl = cl_dil(hk, a.rt);
            if (dl > di) di = dl;                       // laminar floor
        } else {
            double dl = cl_dilw(hk, a.rt);
            if (dl > di) di = dl;
            di *= 2.0;                                  // both wake halves
        }
        a.di = di;
    }

    a.de = (3.15 + 1.72 / (hk - 1.0)) * t + d;
    if (a.de > 12.0 * t) a.de = 12.0 * t;
    return a;
}

// One marched BL station: primary state + type + cached arc position.
struct Stn {
    double x;    // xi  (arc length from stagnation)
    double u;    // Ue
    double t;    // theta
    double d;    // dstar
    double s;    // amplification N (laminar) or sqrt(Ctau) (turbulent)
    bool   turb; // false = laminar, true = turbulent
};

// Midpoint skin friction for the momentum-equation Cf term (BLMID).
double cf_mid(double hka, double rta, bool turb) {
    if (!turb) return cl_cfl(hka, rta);
    double cf = cl_cft(hka, rta), cfl = cl_cfl(hka, rta);
    return (cfl > cf) ? cfl : cf;
}

// The three interval residuals (REZC amplification/lag, REZT momentum, REZH
// shape), faithfully from BLDIF at M=0 (Hca=0, Hwa=0). ityp: 1 lam, 2 turb,
// 3 wake. Station 1 is upstream (known), station 2 the current unknown.
void bl_residual(const Stn& s1, const Stn& s2, int ityp, double Re,
                 double acrit, double R[3]) {
    Aux a1 = bl_aux(s1.u, s1.t, s1.d, s1.s, ityp, Re);
    Aux a2 = bl_aux(s2.u, s2.t, s2.d, s2.s, ityp, Re);

    double ulog = std::log(s2.u / s1.u);
    double xlog = std::log(s2.x / s1.x);
    double tlog = std::log(s2.t / s1.t);
    double hlog = std::log(a2.hs / a1.hs);

    // local upwinding factor UPW
    double hdcon = (ityp == 3 ? 1.0 : 5.0) / (a2.hk * a2.hk);
    double arg = std::fabs((a2.hk - 1.0) / (a1.hk - 1.0));
    double hl = std::log(arg);
    double hlsq = std::min(hl * hl, 15.0);
    double upw = 1.0 - 0.5 * std::exp(-hlsq * hdcon);

    // ---- momentum (REZT) ----
    double ha = 0.5 * (a1.h + a2.h);
    double ta = 0.5 * (s1.t + s2.t), xa = 0.5 * (s1.x + s2.x);
    double cfm = cf_mid(0.5 * (a1.hk + a2.hk), 0.5 * (a1.rt + a2.rt), s2.turb);
    double cfx_m = 0.5 * cfm * xa / ta
                 + 0.25 * (a1.cf * s1.x / s1.t + a2.cf * s2.x / s2.t);
    R[1] = tlog + (ha + 2.0) * ulog - xlog * 0.5 * cfx_m;

    // ---- shape / KE (REZH) ----
    double xot1 = s1.x / s1.t, xot2 = s2.x / s2.t;
    double dix = (1.0 - upw) * a1.di * xot1 + upw * a2.di * xot2;
    double cfx_h = (1.0 - upw) * a1.cf * xot1 + upw * a2.cf * xot2;
    R[2] = hlog + (1.0 - ha) * ulog + xlog * (0.5 * cfx_h - dix);

    // ---- first equation (REZC) ----
    if (ityp == 1) {
        double ax = cl_axset(a1.hk, s1.t, a1.rt, s1.s,
                             a2.hk, s2.t, a2.rt, s2.s, acrit);
        R[0] = s2.s - s1.s - ax * (s2.x - s1.x);
    } else {
        double ald = (ityp == 3) ? DLCON : 1.0;
        double sa = (1.0 - upw) * s1.s + upw * s2.s;
        double cqa = (1.0 - upw) * a1.cq + upw * a2.cq;
        double cfa = (1.0 - upw) * a1.cf + upw * a2.cf;
        double hka = (1.0 - upw) * a1.hk + upw * a2.hk;
        double usa = 0.5 * (a1.us + a2.us);
        double dea = 0.5 * (a1.de + a2.de);
        double da = 0.5 * (s1.d + s2.d);
        double rta = 0.5 * (a1.rt + a2.rt);
        double gcc = (ityp == 2) ? GCCON : 0.0;
        double hkc = hka - 1.0 - gcc / rta;
        if (hkc < 0.01) hkc = 0.01;
        double hr = hkc / (GACON * ald * hka);
        double uq = (0.5 * cfa - hr * hr) / (GBCON * da);
        double scc = SCCON * 1.333 / (1.0 + usa);
        double slog = std::log(s2.s / s1.s);
        double dxi = s2.x - s1.x;
        R[0] = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog
             + dea * 2.0 * (uq * dxi - ulog) * DUXCON;
    }
}

// Solve one downstream station with a 4x4 Newton (numerical Jacobian), mirroring
// MRCHUE.  Unknowns: (theta, dstar, s, Ue). The 4th equation is either the
// direct closure (Ue = the prescribed inviscid value) or, once the kinematic
// shape would exceed HMAX, the inverse closure (Hk = HTARG) which lets Ue float
// locally so a laminar/turbulent separation bubble is traversed instead of
// killing the march. Returns false only on a hard numerical failure.
bool march_station(const Stn& s1, Stn& s2, int ityp, double Re, double acrit,
                   int max_it, double hk_lim) {
    const double ue_presc = s2.u;                  // prescribed inviscid edge vel
    const double HMAX = (ityp == 1) ? 3.8 : 2.5;   // direct/inverse switch shape
    bool inverse = false;
    double htarg = 0.0;

    for (int it = 0; it < max_it; ++it) {
        double R0[4];
        bl_residual(s1, s2, ityp, Re, acrit, R0);
        double hk2 = s2.d / s2.t;
        R0[3] = inverse ? (hk2 - htarg) : (s2.u - ue_presc);

        Eigen::Matrix4d J;
        for (int c = 0; c < 4; ++c) {
            Stn sp = s2;
            double base = (c == 0 ? s2.t : c == 1 ? s2.d : c == 2 ? s2.s : s2.u);
            double h = 1e-6 * std::max(std::fabs(base), 1e-7);
            if (c == 0) sp.t += h; else if (c == 1) sp.d += h;
            else if (c == 2) sp.s += h; else sp.u += h;
            double Rp[4];
            bl_residual(s1, sp, ityp, Re, acrit, Rp);
            double hkp = sp.d / sp.t;
            Rp[3] = inverse ? (hkp - htarg) : (sp.u - ue_presc);
            for (int rrow = 0; rrow < 4; ++rrow) J(rrow, c) = (Rp[rrow] - R0[rrow]) / h;
        }
        Eigen::Vector4d rhs(-R0[0], -R0[1], -R0[2], -R0[3]);
        Eigen::Vector4d dx = J.partialPivLu().solve(rhs);
        if (!dx.allFinite()) return false;

        // under-relax on the largest fractional change (MRCHUE)
        double sref = s2.turb ? std::max(s2.s, 1e-3) : 10.0;
        double dmax = std::max({std::fabs(dx(0) / s2.t), std::fabs(dx(1) / s2.d),
                                std::fabs(dx(2) / sref), std::fabs(dx(3) / s2.u)});
        double rlx = (dmax > 0.3) ? 0.3 / dmax : 1.0;

        // decide direct vs inverse from the would-be shape (only in direct mode)
        if (!inverse) {
            double htest = (s2.d + rlx * dx(1)) / (s2.t + rlx * dx(0));
            if (htest > HMAX) {
                inverse = true;
                double hk1 = s1.d / s1.t;
                htarg = (ityp == 1) ? hk1 + 0.03 * (s2.x - s1.x) / s1.t
                                    : hk1 - 0.15 * (s2.x - s1.x) / s1.t;
                if (htarg < HMAX) htarg = HMAX;
                continue;                          // re-linearise in inverse mode
            }
        }

        s2.t += rlx * dx(0);
        s2.d += rlx * dx(1);
        s2.s += rlx * dx(2);
        s2.u += rlx * dx(3);
        if (s2.turb) s2.s = std::min(std::max(s2.s, 1e-7), 0.30);
        if (s2.t <= 0.0 || s2.d <= 0.0 || s2.u <= 0.0 || !std::isfinite(s2.t))
            return false;
        if (s2.d < 1.02 * s2.t) s2.d = 1.02 * s2.t;   // DSLIM-style Hk floor

        if (dmax * rlx <= 1e-5) break;
    }
    double hk = s2.d / s2.t;
    // Laminar separation bubbles legitimately reach high Hk before transition;
    // only the turbulent watchdog gates genuine (trailing-edge) separation.
    double sep_lim = (ityp == 1) ? 16.0 : hk_lim;
    if (!std::isfinite(hk) || hk > sep_lim) return false;
    return true;
}

// March one surface side from the stagnation point to the trailing edge.
// `xs/us` are arc length (from stagnation) and edge speed at each station,
// ordered outward. Fills theta/dstar/H at the TE and the transition x-fraction.
// Returns false if the side separates / fails to converge before the TE.
struct SideOut {
    double theta_te = 0, dstar_te = 0, h_te = 0, ue_te = 0;
    double xtr = 1.0;
    double cdf = 0.0;        // friction-drag contribution (sum Cf*Ue^2*dx_proj)
    bool   ok = false;
    bool   separated = false;
    std::vector<double> xi, dstar;   // per station (stagnation -> TE order)
};

// Compute the side outputs (friction drag, TE momentum/shape/Ue for the
// Squire-Young far-wake drag) from the stations marched up to `last`. If `last`
// is short of the TE the flow separated; theta/dstar are then extrapolated to
// the TE (theta ~ xi^0.5) so the wake formula still yields a finite drag.
SideOut finish_side(const std::vector<Stn>& st, const std::vector<double>& xs,
                    int last, double xte, bool separated, double xtr, double Re) {
    SideOut out;
    out.xtr = xtr;
    out.separated = separated;
    // per-station displacement thickness for the transpiration source; the
    // separated tail (beyond `last`) is extrapolated as dstar ~ xi^0.5.
    int ns = static_cast<int>(xs.size());
    out.xi.assign(ns, 0.0);
    out.dstar.assign(ns, 0.0);
    for (int j = 0; j < ns; ++j) {
        out.xi[j] = xs[j];
        if (j <= last) out.dstar[j] = st[j].d;
        else out.dstar[j] = st[last].d * std::sqrt(std::max(xs[j], 1e-9)
                                                   / std::max(xs[last], 1e-9));
    }
    // friction drag: integral of Cf*Ue^2 along the surface (q_inf = 1).
    double cdf = 0.0;
    for (int j = 1; j <= last; ++j) {
        int ityp = st[j].turb ? 2 : 1;
        Aux a0 = bl_aux(st[j - 1].u, st[j - 1].t, st[j - 1].d, st[j - 1].s,
                        st[j - 1].turb ? 2 : 1, Re);
        Aux a1 = bl_aux(st[j].u, st[j].t, st[j].d, st[j].s, ityp, Re);
        double tau0 = a0.cf * st[j - 1].u * st[j - 1].u;
        double tau1 = a1.cf * st[j].u * st[j].u;
        cdf += 0.5 * (tau0 + tau1) * (st[j].x - st[j - 1].x);
    }
    out.cdf = cdf;

    double th = st[last].t, ds = st[last].d, ue = st[last].u, xl = st[last].x;
    if (separated && xl < xte && xl > 1e-6) {
        double rr = std::sqrt(xte / xl);
        th *= rr;                     // laminar/turbulent theta ~ xi^0.5 growth
        ds *= rr;
    }
    out.theta_te = th;
    out.dstar_te = ds;
    out.h_te = ds / th;
    out.ue_te = ue;
    out.ok = (last >= 2);
    return out;
}

SideOut march_side(const std::vector<double>& xs, const std::vector<double>& us,
                   const std::vector<double>& xcs,  // chordwise x at each station
                   double xte, double Re, const Options& opt) {
    int ns = static_cast<int>(xs.size());
    if (ns < 3) return SideOut{};
    double acrit = opt.Ncrit;
    double hk_lim = opt.hk_sep_limit;

    std::vector<Stn> st(ns);
    // Thwaites similarity initialisation at the first station (BULE = 1).
    double x1 = std::max(xs[0], 1e-6), u1 = std::max(us[0], 1e-4);
    double tsq = 0.075 * x1 / (u1 * Re);
    st[0].x = x1; st[0].u = u1;
    st[0].t = std::sqrt(std::max(tsq, 1e-12));
    st[0].d = 2.2 * st[0].t;
    st[0].s = 0.0;            // amplification N
    st[0].turb = false;

    bool turb = false;
    double xtr = 1.0;
    for (int j = 1; j < ns; ++j) {
        st[j].x = std::max(xs[j], st[j - 1].x + 1e-9);
        st[j].u = std::max(us[j], 1e-4);
        st[j].t = st[j - 1].t;       // continuity guess from upstream
        st[j].d = st[j - 1].d;
        st[j].s = st[j - 1].s;
        st[j].turb = turb;

        bool ok;
        if (!turb) {
            ok = march_station(st[j - 1], st[j], 1, Re, acrit,
                               opt.max_bl_newton, hk_lim);
            if (ok && st[j].s >= acrit) {     // transition -> re-solve turbulent
                turb = true;
                xtr = xcs[j];
                Aux ap = bl_aux(st[j - 1].u, st[j - 1].t, st[j - 1].d, 0.0, 2, Re);
                st[j - 1].s = ap.cq;          // seed upstream Ctau (equilibrium)
                st[j - 1].turb = true;
                st[j].turb = true;
                st[j].s = 0.03;
                ok = march_station(st[j - 1], st[j], 2, Re, acrit,
                                   opt.max_bl_newton, hk_lim);
            }
        } else {
            ok = march_station(st[j - 1], st[j], 2, Re, acrit,
                               opt.max_bl_newton, hk_lim);
        }
        if (std::getenv("XF_BL")) {
            Aux a = bl_aux(st[j].u, st[j].t, st[j].d, st[j].s, st[j].turb ? 2 : 1, Re);
            std::fprintf(stderr, "  j=%3d x=%.4f xc=%.4f ue=%.4f th=%.5f H=%.3f Rt=%.1f n/ct=%.4f %s %s\n",
                         j, st[j].x, xcs[j], st[j].u, st[j].t, a.h, a.rt, st[j].s,
                         st[j].turb ? "T" : "L", ok ? "" : "FAIL");
        }
        if (!ok) return finish_side(st, xs, j - 1, xte, true, xtr, Re);
    }
    return finish_side(st, xs, ns - 1, xte, false, xtr, Re);
}

}  // namespace

// ===========================================================================
//  [4-6] viscous driver: edge velocity -> BL march -> Squire-Young drag.
//  Weak coupling by default (BL marched on the inviscid Ue; inviscid loads),
//  with the strong transpiration interaction available behind opt.strong_coupling.
// ===========================================================================

// March both BL sides given a signed edge-velocity field. Fills the per-panel
// displacement thickness, returns the two SideOuts, and reports success.
static bool bl_pass(Solver::Impl& m, const std::vector<double>& ue, double Re,
                    std::vector<double>& dstar, SideOut& sA, SideOut& sB,
                    std::vector<int>& idxA, std::vector<int>& idxB) {
    const int N = m.N;
    int kstag = -1; double best_xc = 1e9;
    for (int i = 1; i < N; ++i)
        if (ue[i - 1] * ue[i] < 0.0) {
            double xcm = 0.5 * (m.xc[i - 1] + m.xc[i]);
            if (xcm < best_xc) { best_xc = xcm; kstag = i; }
        }
    if (kstag < 2 || kstag > N - 2) return false;
    double a1 = std::fabs(ue[kstag - 1]), a2 = std::fabs(ue[kstag]);
    double frac = (a1 + a2 > 0) ? a1 / (a1 + a2) : 0.5;
    double s_stag = m.sarc[kstag - 1] + frac * (m.sarc[kstag] - m.sarc[kstag - 1]);

    idxA.clear(); idxB.clear();
    std::vector<double> xA, uA, cA, xB, uB, cB;
    for (int j = kstag - 1; j >= 0; --j) {
        idxA.push_back(j); xA.push_back(s_stag - m.sarc[j]);
        uA.push_back(std::fabs(ue[j])); cA.push_back(m.xc[j]);
    }
    for (int j = kstag; j < N; ++j) {
        idxB.push_back(j); xB.push_back(m.sarc[j] - s_stag);
        uB.push_back(std::fabs(ue[j])); cB.push_back(m.xc[j]);
    }
    sA = march_side(xA, uA, cA, xA.back(), Re, m.opt);
    sB = march_side(xB, uB, cB, xB.back(), Re, m.opt);
    if (!sA.ok || !sB.ok) return false;
    for (std::size_t k = 0; k < idxA.size(); ++k) dstar[idxA[k]] = sA.dstar[k];
    for (std::size_t k = 0; k < idxB.size(); ++k) dstar[idxB[k]] = sB.dstar[k];
    return true;
}

static Result viscous_solve(Solver::Impl& m, double alpha, double Re) {
    Result r;
    r.alpha = alpha;
    const int N = m.N;
    std::vector<double> ue_inv(N), ue(N), dstar(N, 0.0);
    inviscid_ue(m, alpha, ue_inv);
    ue = ue_inv;                                  // signed tangential edge velocity

    SideOut sA, sB;
    std::vector<int> idxA, idxB;

    if (!m.opt.strong_coupling) {
        // ---- weak coupling: one BL pass on the inviscid edge velocity ----
        if (!bl_pass(m, ue_inv, Re, dstar, sA, sB, idxA, idxB)) {
            r.converged = false; return r;
        }
        forces_from_cp(m, ue_inv, alpha, r.cl, r.cm);   // inviscid section loads
    } else {
        // ---- strong coupling: under-relaxed transpiration interaction ----
        // EXPERIMENTAL: the explicit fixed point is unstable where a laminar
        // separation bubble drives a steep d(dstar)/ds; robust use needs the
        // implicit Newton. Kept for non-bubble (higher-Re) regimes and study.
        Eigen::VectorXd wn(N);
        double relax = (m.opt.relax > 0.0 && m.opt.relax <= 1.0) ? m.opt.relax : 0.3;
        bool ok = false;
        for (int iter = 0; iter < m.opt.max_newton; ++iter) {
            if (!bl_pass(m, ue, Re, dstar, sA, sB, idxA, idxB)) {
                r.converged = false; return r;
            }
            wn.setZero();
            auto blow = [&](const std::vector<int>& idx, const SideOut& s) {
                int ns = static_cast<int>(idx.size());
                for (int k = 0; k < ns; ++k) {
                    int p = idx[k];
                    int km = std::max(k - 1, 0), kp = std::min(k + 1, ns - 1);
                    double mm = std::fabs(ue[idx[km]]) * dstar[idx[km]];
                    double mp = std::fabs(ue[idx[kp]]) * dstar[idx[kp]];
                    double dxi = s.xi[kp] - s.xi[km];
                    wn(p) = (dxi > 1e-9) ? (mp - mm) / dxi : 0.0;
                }
            };
            blow(idxA, sA); blow(idxB, sB);
            Eigen::VectorXd duv = m.Mresp * wn;
            double dmax = 0.0;
            for (int i = 0; i < N; ++i) {
                double unew = ue_inv[i] + duv(i);
                double step = relax * (unew - ue[i]);
                double cap = 0.1 * std::max(std::fabs(ue_inv[i]), 0.1);
                step = std::min(std::max(step, -cap), cap);
                dmax = std::max(dmax, std::fabs(unew - ue[i]));
                ue[i] += step;
            }
            ok = true;
            if (dmax < m.opt.tol) break;
        }
        if (!ok) { r.converged = false; return r; }
        forces_from_cp(m, ue, alpha, r.cl, r.cm);       // decambered loads
    }

    // ---- far-wake Squire-Young profile drag, summed over both surfaces ----
    auto sqy = [](const SideOut& s) {
        if (s.theta_te <= 0.0 || s.ue_te <= 0.0) return 0.0;
        return 2.0 * s.theta_te * std::pow(s.ue_te, 0.5 * (5.0 + s.h_te));
    };
    double cd = sqy(sA) + sqy(sB);
    double cdf = sA.cdf + sB.cdf;

    r.cd  = cd;
    r.cdf = std::min(cdf, cd);
    r.cdp = std::max(cd - r.cdf, 0.0);
    r.xtr_top = sA.xtr;
    r.xtr_bot = sB.xtr;
    r.separated = sA.separated || sB.separated;
    r.converged = std::isfinite(cd) && cd > 0.0 && cd < m.opt.cd_sane_max
                  && !r.separated;
    return r;
}

}  // namespace xfoil
}  // namespace aero
