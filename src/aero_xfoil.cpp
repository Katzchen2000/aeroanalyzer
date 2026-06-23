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

// One marched BL station: primary state + type + cached arc position. Defined
// here (ahead of Solver::Impl) so the coupled-Newton buffers can hold an array
// of them.  x = xi (arc length from stagnation), u = Ue, t = theta, d = dstar,
// s = amplification N (laminar) or sqrt(Ctau) (turbulent).
struct Stn {
    double x = 0, u = 0, t = 0, d = 0, s = 0;
    bool   turb = false;
};
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
    // full (source,vortex) strength vectors for those two unit solutions, and
    // the per-panel blowing-response strengths (column j = A^{-1} e_j). These let
    // the wake coupling evaluate the inviscid field and the blowing response at
    // off-surface (wake) points.  Sizes: sol_a/sol_b = N+1, SrcResp = (N+1) x N.
    Eigen::VectorXd sol_a, sol_b;
    Eigen::MatrixXd SrcResp;
    bool geom_ok = false;

    // ---- warm-start BL state (per surface node, retained across solves) ----
    bool have_state = false;
    double state_Re = -1.0;
    std::vector<double> st_theta, st_dstar, st_ue;  // size N (indexed by panel)

    // ---- coupled viscous-inviscid Newton working buffers ----
    //   All sized once in build_inviscid() to the panel count N (an upper bound
    //   on the surface-station count), so coupled_solve() never allocates in its
    //   assembly/solve loops -- safe to run one Solver per OpenMP thread.
    std::vector<int>    cp_pan;     // full station -> panel index
    std::vector<double> cp_vti;     // full station -> BL-frame sign (+-1)
    std::vector<double> cp_uinv;    // full station -> inviscid edge SPEED (>=0)
    std::vector<int>    cp_side;    // full station -> 0/1 (which side)
    std::vector<int>    cp_lpos;    // full station -> position along its side
    std::vector<int>    cp_prev;    // full station -> upstream full-station index
    std::vector<int>    cp_unk;     // full station -> Newton unknown index (-1 = frozen)
    std::vector<Stn>    cp_st;      // full station -> current state (u,t,d,s,turb,x)
    Eigen::MatrixXd     cp_Mr;      // vti-scaled, reordered Mresp (transpiration)
    Eigen::MatrixXd     cp_G;       // d/d(xi) operator: w = G * m (per side)
    Eigen::MatrixXd     cp_C;       // BL-frame edge-velocity coupling = cp_Mr*cp_G
    Eigen::MatrixXd     cp_J;       // global Jacobian, 3*Nunk x 3*Nunk
    Eigen::VectorXd     cp_R, cp_dx;
    Eigen::PartialPivLU<Eigen::MatrixXd> cp_lu;

    // ---- minimal wake (rebuilt per alpha; geometry depends on flow direction) ----
    //   A short source line streaming off the TE that carries the BL mass defect
    //   downstream so the surface transpiration is no longer amputated at the TE.
    int wk_n = 0;                          // number of wake panels in use
    std::vector<double> wk_x, wk_y;        // wake node coords, size wk_n+1
    std::vector<double> wk_xc, wk_yc;      // wake panel midpoints, size wk_n
    std::vector<double> wk_c, wk_s, wk_l;  // wake panel tangent (c,s) + length
    std::vector<double> wk_uinv;           // inviscid tangential vel at midpoints

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

// Global-frame velocity at (px,py) from a unit constant-strength SOURCE panel
// described by start node (x0,y0), unit tangent (c,s) and length L. Same kernel
// as panel_influence but for an arbitrary (e.g. wake) panel not in the surface
// arrays. Never self-evaluated (wake points are off the surface).
static inline void source_panel_vel(double x0, double y0, double c, double s,
                                    double L, double px, double py,
                                    double& su, double& sv) {
    double dx = px - x0, dy = py - y0;
    double xp = dx * c + dy * s;
    double yp = -dx * s + dy * c;
    double r1 = std::sqrt(xp * xp + yp * yp);
    double r2 = std::sqrt((xp - L) * (xp - L) + yp * yp);
    double ln = std::log(std::max(r1, 1e-12) / std::max(r2, 1e-12));
    double beta = std::atan2(yp, xp - L) - std::atan2(yp, xp);
    double inv = 1.0 / TWO_PI;
    double sul = ln * inv, svl = beta * inv;
    su = sul * c - svl * s;
    sv = sul * s + svl * c;
}

// Global-frame velocity at an off-surface point (px,py) induced by a strength
// vector sol (N surface sources + 1 trailing vortex, size N+1) plus a freestream
// (vx_inf,vy_inf). Used to evaluate the inviscid field and the blowing response
// at wake collocation points.
static inline void field_velocity(const Solver::Impl& m, const Eigen::VectorXd& sol,
                                  double vx_inf, double vy_inf,
                                  double px, double py, double& u, double& v) {
    u = vx_inf; v = vy_inf;
    double gam = sol(m.N);
    for (int j = 0; j < m.N; ++j) {
        double su, sv, vu, vv;
        panel_influence(m, j, px, py, false, su, sv, vu, vv);
        u += sol(j) * su + gam * vu;
        v += sol(j) * sv + gam * vv;
    }
}

// Build the minimal wake geometry off the trailing edge for the current flow
// direction, and evaluate the inviscid tangential edge velocity along it. The
// wake leaves the TE along the freestream direction with geometrically growing
// panels out to ~1 chord. Idempotent; sizes are bounded by the panel count.
static void build_wake(Solver::Impl& m, double alpha) {
    const int N = m.N;
    int NWP = N / 8 + 2;
    const double ca = std::cos(alpha), sa = std::sin(alpha);
    const double dirx = ca, diry = sa;            // wake streams along freestream
    double tex = 0.5 * (m.X[0] + m.X[N]);         // trailing-edge point
    double tey = 0.5 * (m.Y[0] + m.Y[N]);

    // geometric spacing: first panel ~ last surface panel, total length ~1 chord
    double dl0 = std::max({m.len[0], m.len[N - 1], 1e-4});
    const double Lw = 1.0;
    double r = 1.0;                               // bisection for growth ratio
    {
        double lo = 1.0 + 1e-6, hi = 3.0;
        auto tot = [&](double rr) {
            return (std::fabs(rr - 1.0) < 1e-9) ? dl0 * NWP
                                                : dl0 * (std::pow(rr, NWP) - 1.0) / (rr - 1.0);
        };
        if (tot(lo) > Lw) { r = lo; }
        else if (tot(hi) < Lw) { r = hi; }
        else { for (int it = 0; it < 80; ++it) { double mid = 0.5 * (lo + hi);
                   if (tot(mid) < Lw) lo = mid; else hi = mid; } r = 0.5 * (lo + hi); }
    }

    m.wk_n = NWP;
    m.wk_x.assign(NWP + 1, 0.0); m.wk_y.assign(NWP + 1, 0.0);
    m.wk_xc.assign(NWP, 0.0);    m.wk_yc.assign(NWP, 0.0);
    m.wk_c.assign(NWP, dirx);    m.wk_s.assign(NWP, diry);
    m.wk_l.assign(NWP, 0.0);     m.wk_uinv.assign(NWP, 0.0);
    m.wk_x[0] = tex; m.wk_y[0] = tey;
    double dl = dl0;
    for (int k = 0; k < NWP; ++k) {
        m.wk_x[k + 1] = m.wk_x[k] + dl * dirx;
        m.wk_y[k + 1] = m.wk_y[k] + dl * diry;
        m.wk_xc[k] = 0.5 * (m.wk_x[k] + m.wk_x[k + 1]);
        m.wk_yc[k] = 0.5 * (m.wk_y[k] + m.wk_y[k + 1]);
        m.wk_l[k] = dl;
        dl *= r;
    }
    // inviscid tangential edge velocity along the wake (alpha-blended strengths)
    Eigen::VectorXd sol = ca * m.sol_a + sa * m.sol_b;
    for (int k = 0; k < NWP; ++k) {
        double u, v;
        field_velocity(m, sol, ca, sa, m.wk_xc[k], m.wk_yc[k], u, v);
        m.wk_uinv[k] = u * m.wk_c[k] + v * m.wk_s[k];
    }
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
    auto solve_unit = [&](double vx, double vy, Eigen::VectorXd& ue,
                          Eigen::VectorXd& solout) {
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
        solout = sol;                                        // retain full strengths
        if (std::getenv("XF_DEBUG"))
            std::fprintf(stderr, "[xf] unit(%.0f,%.0f): gamma=%.5f ue[0]=%.4f ue[N-1]=%.4f sigma0=%.4f\n",
                         vx, vy, gam, ue(0), ue(N - 1), sol(0));
    };
    solve_unit(1.0, 0.0, m.uea, m.sol_a);
    solve_unit(0.0, 1.0, m.ueb, m.sol_b);

    // Edge-velocity response to wall transpiration (the viscous coupling matrix).
    // A unit outward blowing velocity at panel j enters the tangency RHS as e_j;
    // back-substituting through the factored LU gives the source/vortex response,
    // whose tangential trace is column j of Mresp. Built once, reused every alpha
    // and every interaction sweep. (Kutta RHS unaffected by blowing.)
    m.Mresp.resize(N, N);
    m.SrcResp.resize(N + 1, N);
    Eigen::VectorXd rhs(N + 1), sol(N + 1);
    for (int j = 0; j < N; ++j) {
        rhs.setZero();
        rhs(j) = 1.0;
        sol = m.lu.solve(rhs);
        m.SrcResp.col(j) = sol;                  // strengths for blowing at panel j
        double gam = sol(N);
        for (int i = 0; i < N; ++i)
            m.Mresp(i, j) = (m.At.row(i) * sol.head(N))(0) + gam * m.Avt(i);
    }

    // Preallocate the coupled-Newton buffers to the panel count (>= station
    // count). Done once per geometry so coupled_solve() is allocation-free.
    m.cp_pan.assign(N, 0);  m.cp_vti.assign(N, 0.0);  m.cp_uinv.assign(N, 0.0);
    m.cp_side.assign(N, 0); m.cp_lpos.assign(N, 0);   m.cp_prev.assign(N, -1);
    m.cp_unk.assign(N, -1); m.cp_st.assign(N, Stn{});
    m.cp_Mr.resize(N, N);
    m.cp_G.resize(N, N);
    m.cp_C.resize(N, N);
    m.cp_J.resize(3 * N, 3 * N);
    m.cp_R.resize(3 * N);
    m.cp_dx.resize(3 * N);
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
    // Full marched state, retained so the coupled Newton can linearise about it.
    std::vector<Stn> stn;            // stations 0..last (0 = stagnation IC)
    int  last  = 0;                  // index of the last marched station
    int  itran = 1 << 30;            // station index of transition (>= size = none)
};

// Compute the side outputs (friction drag, TE momentum/shape/Ue for the
// Squire-Young far-wake drag) from the stations marched up to `last`. If `last`
// is short of the TE the flow separated; theta/dstar are then extrapolated to
// the TE (theta ~ xi^0.5) so the wake formula still yields a finite drag.
SideOut finish_side(const std::vector<Stn>& st, const std::vector<double>& xs,
                    int last, double xte, bool separated, double xtr, double Re,
                    int itran = (1 << 30)) {
    SideOut out;
    out.xtr = xtr;
    out.separated = separated;
    // Retain the marched state (stations 0..last) for the coupled Newton.
    out.stn.assign(st.begin(), st.begin() + (last + 1));
    out.last = last;
    out.itran = itran;
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
    int itran = ns;                  // station index where flow turns turbulent
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
                itran = j;
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
        if (!ok) return finish_side(st, xs, j - 1, xte, true, xtr, Re, itran);
    }
    return finish_side(st, xs, ns - 1, xte, false, xtr, Re, itran);
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

// ===========================================================================
//  Strong viscous-inviscid coupling: the simultaneous implicit Newton.
//
//  Faithful to Drela's XFOIL VISCAL/SETBL/BLSOLV/UPDATE, surface-only (no
//  marched wake; far-wake drag stays Squire-Young). The defining difference
//  from the old explicit fixed point is that the boundary-layer conservation
//  equations and the displacement->edge-velocity closure are solved together
//  in ONE Newton system, so the steep d(dstar)/ds of a laminar bubble can no
//  longer run the iteration away.
//
//  Unknowns per surface station k (>=1, the k=0 stagnation point is a frozen
//  Thwaites start): (s_k, theta_k, m_k) with m = Ue*dstar the mass defect and
//      Ue_a = Uinv_a + sum_c C(a,c) * m_c,      dstar_a = m_a / Ue_a.
//  C = (vti * Mresp) * (d/dxi): the transpiration model. A growing displacement
//  thickness is an outward wall blowing w = d(Ue*dstar)/dxi; Mresp maps that
//  blowing to the inviscid edge velocity (built once in build_inviscid), and
//  d/dxi is a per-side finite-difference operator. This is the same physics as
//  the explicit loop, but here C enters the Jacobian analytically (the dense
//  VM block of BLSOLV) instead of being chased by under-relaxed substitution.
//
//  The residual is the 3 BL interval equations (BLDIF) per station; its local
//  sensitivities are taken numerically (reusing bl_residual), and the global
//  m-coupling is chained through C analytically. Newton steps are damped by the
//  XFOIL trust-region limiter (worst normalised change clamped to [-0.5,+1.5]).
// ===========================================================================
namespace {
constexpr double COUPLING_SIGN = 1.0;  // Mresp built with +1 = outward blowing

// R0 plus the eight local residual sensitivities for the interval (s1 -> s2):
// d R / d{s1,t1,d1,u1, s2,t2,d2,u2}. Forward differences; bl_residual is cheap.
struct LinBlk {
    double R0[3];
    double Ps1[3], Pt1[3], Pd1[3], Pu1[3];   // wrt station 1 (upstream)
    double Ps2[3], Pt2[3], Pd2[3], Pu2[3];   // wrt station 2 (current)
};
void linearize(const Stn& s1, const Stn& s2, int ityp, double Re, double acrit,
               LinBlk& L) {
    bl_residual(s1, s2, ityp, Re, acrit, L.R0);
    auto col = [&](const Stn& a, const Stn& b, int which, double Stn::*fld,
                   double P[3]) {
        Stn ap = a, bp = b;
        double base = (which == 1 ? a.*fld : b.*fld);
        double h = 1e-6 * std::max(std::fabs(base), 1e-7);
        if (which == 1) ap.*fld += h; else bp.*fld += h;
        double Rp[3];
        bl_residual(ap, bp, ityp, Re, acrit, Rp);
        for (int i = 0; i < 3; ++i) P[i] = (Rp[i] - L.R0[i]) / h;
    };
    col(s1, s2, 1, &Stn::s, L.Ps1); col(s1, s2, 1, &Stn::t, L.Pt1);
    col(s1, s2, 1, &Stn::d, L.Pd1); col(s1, s2, 1, &Stn::u, L.Pu1);
    col(s1, s2, 2, &Stn::s, L.Ps2); col(s1, s2, 2, &Stn::t, L.Pt2);
    col(s1, s2, 2, &Stn::d, L.Pd2); col(s1, s2, 2, &Stn::u, L.Pu2);
}
}  // namespace

static Result coupled_solve(Solver::Impl& m, double alpha, double Re) {
    Result r;
    r.alpha = alpha;
    const int N = m.N;
    const Options& opt = m.opt;
    const double acrit = opt.Ncrit;

    // ---- 1) initial weak march: stagnation split, transition, state guess ----
    std::vector<double> ue_inv(N), dstar(N, 0.0);
    inviscid_ue(m, alpha, ue_inv);                 // signed panel tangential vel
    SideOut s0[2];
    std::vector<int> idx[2];
    if (!bl_pass(m, ue_inv, Re, dstar, s0[0], s0[1], idx[0], idx[1])) {
        r.converged = false; return r;
    }

    // ---- 2) build the full station list (both sides, k=0 = stagnation) ----
    int F = 0;                                     // number of full stations
    int itranF[2];
    for (int is = 0; is < 2; ++is) {
        const SideOut& s = s0[is];
        int nst = s.last + 1;
        itranF[is] = -1;
        for (int k = 0; k < nst; ++k) {
            int a = F++;
            int p = idx[is][k];
            m.cp_pan[a]  = p;
            m.cp_vti[a]  = (ue_inv[p] >= 0.0) ? 1.0 : -1.0;
            m.cp_uinv[a] = std::fabs(ue_inv[p]);
            m.cp_side[a] = is;
            m.cp_lpos[a] = k;
            m.cp_prev[a] = (k == 0) ? -1 : (a - 1);
            m.cp_st[a]   = s.stn[k];
            if (k >= s.itran && itranF[is] < 0) itranF[is] = a;
        }
    }
    // assign Newton unknown indices to every k>=1 station (k=0 frozen)
    int Nu = 0;
    for (int a = 0; a < F; ++a)
        m.cp_unk[a] = (m.cp_lpos[a] >= 1) ? Nu++ : -1;
    if (Nu < 4) { r.converged = false; return r; }
    const int NX = 3 * Nu;

    // arc length per full station (BL frame, increasing downstream)
    std::vector<double> xi(F);
    for (int a = 0; a < F; ++a) xi[a] = m.cp_st[a].x;

    // minimal wake geometry + inviscid edge velocity along it
    build_wake(m, alpha);
    if (std::getenv("XF_NEWT")) {
        std::fprintf(stderr, "[wake] NWP=%d  TE=(%.4f,%.4f) dir=(%.3f,%.3f) len0=%.5f lenW=%.5f total=%.4f\n",
                     m.wk_n, m.wk_x[0], m.wk_y[0], m.wk_c[0], m.wk_s[0],
                     m.wk_l[0], m.wk_l[m.wk_n - 1],
                     m.wk_x[m.wk_n] - m.wk_x[0]);
        std::fprintf(stderr, "[wake] uinv:");
        for (int k = 0; k < m.wk_n; k += std::max(1, m.wk_n / 8))
            std::fprintf(stderr, " [%d]=%.4f", k, m.wk_uinv[k]);
        std::fprintf(stderr, "  (last=%.4f)\n", m.wk_uinv[m.wk_n - 1]);
    }

    // ---- 3) coupling operator  C = (vti . Mresp) * (d/dxi) ----
    // Mr(a,b) = sign * vti_a * Mresp(pan_a, pan_b); G is the per-side central
    // difference  w_b = d m / d xi  so that  Ue_a = Uinv_a + (Mr*G m)_a.
    auto Mr = m.cp_Mr.topLeftCorner(F, F);
    auto G  = m.cp_G.topLeftCorner(F, F);
    for (int a = 0; a < F; ++a)
        for (int b = 0; b < F; ++b)
            Mr(a, b) = COUPLING_SIGN * m.cp_vti[a] * m.Mresp(m.cp_pan[a], m.cp_pan[b]);
    // Floor the arc spacing used for the transpiration gradient. The cosine
    // mesh puts ~5e-4-chord panels at the TE, so a raw d/dxi over a marginal
    // (near-separation) dstar profile produces a nonphysical source spike. The
    // displacement source physically cannot vary faster than ~1% chord; flooring
    // dd here bounds the source the same way a wake's larger spacing would, and
    // keeps C linear (the floor is part of the operator, exact in the Jacobian).
    double dd_floor = 0.02;
    if (const char* e = std::getenv("XF_DSFLOOR")) dd_floor = std::atof(e);
    G.setZero();
    for (int a = 0; a < F; ++a) {
        int prev = m.cp_prev[a];
        // next station on the same side (contiguous indexing within a side)
        int nxt = (a + 1 < F && m.cp_side[a + 1] == m.cp_side[a]) ? a + 1 : -1;
        if (prev >= 0 && nxt >= 0) {                // central difference
            double dd = std::max(xi[nxt] - xi[prev], dd_floor);
            G(a, nxt) += 1.0 / dd;  G(a, prev) -= 1.0 / dd;
        } else if (nxt >= 0) {                       // forward (k=0 endpoint)
            double dd = std::max(xi[nxt] - xi[a], dd_floor);
            G(a, nxt) += 1.0 / dd;  G(a, a) -= 1.0 / dd;
        } else if (prev >= 0) {                      // backward (TE endpoint)
            double dd = std::max(xi[a] - xi[prev], dd_floor);
            G(a, a) += 1.0 / dd;    G(a, prev) -= 1.0 / dd;
        }
    }
    m.cp_C.topLeftCorner(F, F).noalias() = Mr * G;

    // ---- 4) initial unknown state: m = Ue*dstar from the weak march ----
    std::vector<double> mfull(F), uefull(F), dsfull(F);
    for (int a = 0; a < F; ++a) {
        uefull[a] = std::max(m.cp_st[a].u, 1e-6);
        dsfull[a] = m.cp_st[a].d;
        mfull[a]  = uefull[a] * dsfull[a];
    }

    if (std::getenv("XF_NEWT")) {  // residual on the RAW marched state (floated Ue)
        double rsum = 0.0, rmx = 0.0; int armx = -1;
        for (int a = 0; a < F; ++a) {
            if (m.cp_lpos[a] < 1) continue;
            int pa = m.cp_prev[a];
            int ityp = m.cp_st[a].turb ? 2 : 1;
            double R3[3];
            bl_residual(m.cp_st[pa], m.cp_st[a], ityp, Re, acrit, R3);
            for (int i = 0; i < 3; ++i) { rsum += R3[i]*R3[i];
                if (std::fabs(R3[i]) > rmx) { rmx = std::fabs(R3[i]); armx = a; } }
        }
        std::fprintf(stderr, "[rawR] ||R(marched state, floated Ue)||=%.3e  max=%.3e at a=%d (lpos %d turb %d H=%.2f)\n",
                     std::sqrt(rsum), rmx, armx, armx>=0?m.cp_lpos[armx]:-1,
                     armx>=0?(int)m.cp_st[armx].turb:-1,
                     armx>=0?dsfull[armx]/std::max(m.cp_st[armx].t,1e-9):0.0);
    }
    if (std::getenv("XF_NEWT")) {  // raw weak-march profile near both TEs
        int aHmax = 0; double Hmax = 0;
        for (int a = 0; a < F; ++a) {
            double H = dsfull[a] / std::max(m.cp_st[a].t, 1e-9);
            if (H > Hmax) { Hmax = H; aHmax = a; }
        }
        std::fprintf(stderr, "[wkmrch] RAW max H=%.3f at a=%d (side %d lpos %d x=%.4f turb=%d ue=%.4f th=%.6f ds=%.6f)\n",
                     Hmax, aHmax, m.cp_side[aHmax], m.cp_lpos[aHmax], m.cp_st[aHmax].x,
                     (int)m.cp_st[aHmax].turb, uefull[aHmax], m.cp_st[aHmax].t, dsfull[aHmax]);
        for (int is = 0; is < 2; ++is) {
            int last = -1;
            for (int b = 0; b < F; ++b) if (m.cp_side[b] == is) last = b;
            std::fprintf(stderr, "[wkmrch] side %d last 5 stations (raw weak march):\n", is);
            for (int a = std::max(0, last - 4); a <= last; ++a)
                std::fprintf(stderr, "[wkmrch]   a=%d lpos=%d x=%.5f ue=%.4f th=%.6f ds=%.6f H=%.3f m=%.6f turb=%d\n",
                             a, m.cp_lpos[a], m.cp_st[a].x, uefull[a], m.cp_st[a].t,
                             dsfull[a], dsfull[a] / std::max(m.cp_st[a].t, 1e-9),
                             mfull[a], (int)m.cp_st[a].turb);
        }
    }

    // Continuation factor on the displacement coupling: lam ramps 0 -> 1 so the
    // mass defect relaxes gradually instead of being slammed with the (steep,
    // un-relaxed) weak-march dm/ds spike in one Newton step. Ue = Uinv + lam*C*m.
    double lam = 1.0;

    // recompute Ue (=Uinv + lam*C m) and dstar for every station from the m field
    auto refresh_state = [&]() {
        for (int a = 0; a < F; ++a) {
            double ue = m.cp_uinv[a];
            double acc = 0.0;
            for (int b = 0; b < F; ++b) acc += m.cp_C(a, b) * mfull[b];
            ue += lam * acc;
            if (ue < 1e-5) ue = 1e-5;              // keep edge speed positive
            uefull[a] = ue;
            m.cp_st[a].u = ue;
            if (m.cp_lpos[a] == 0) {
                // frozen stagnation: refresh Thwaites IC from the new Ue
                double x1 = std::max(xi[a], 1e-6);
                double tsq = 0.075 * x1 / (ue * Re);
                m.cp_st[a].t = std::sqrt(std::max(tsq, 1e-12));
                m.cp_st[a].d = 2.2 * m.cp_st[a].t;
                m.cp_st[a].s = 0.0;
                dsfull[a] = m.cp_st[a].d;
                mfull[a]  = ue * dsfull[a];
            } else {
                double ds = mfull[a] / ue;
                if (ds < 1.02 * m.cp_st[a].t) ds = 1.02 * m.cp_st[a].t;
                dsfull[a] = ds;
                m.cp_st[a].d = ds;
            }
        }
    };
    refresh_state();

    if (std::getenv("XF_NEWT")) {
        // explicit blowing field w = G m and its arc spacing near LE/TE
        std::vector<double> w(F, 0.0);
        for (int a = 0; a < F; ++a)
            for (int b = 0; b < F; ++b) w[a] += m.cp_G(a, b) * mfull[b];
        double wmn = 1e9, wmx = -1e9; int awmx = 0;
        for (int a = 0; a < F; ++a) { wmn = std::min(wmn, w[a]);
            if (w[a] > wmx) { wmx = w[a]; awmx = a; } }
        int te = F - 1;
        double dd_le = xi[2] - xi[0], dd_te = xi[te] - xi[te - 2];
        std::fprintf(stderr, "[init] w range [%.3f,%.3f]  dxi_le=%.5f dxi_te=%.5f  uinv_te=%.4f\n",
                     wmn, wmx, dd_le, dd_te, m.cp_uinv[te]);
        // locate the w spike and dump m / dxi in its neighborhood
        std::fprintf(stderr, "[wspk] argmax|w| a=%d side=%d lpos=%d pan=%d turb=%d\n",
                     awmx, m.cp_side[awmx], m.cp_lpos[awmx], m.cp_pan[awmx],
                     (int)m.cp_st[awmx].turb);
        for (int a = std::max(0, awmx - 2); a <= std::min(F - 1, awmx + 2); ++a)
            std::fprintf(stderr, "[wspk]   a=%d lpos=%d xi=%.6f m=%.6f w=%.4f uinv=%.4f ue=%.4f th=%.6f ds=%.6f\n",
                         a, m.cp_lpos[a], xi[a], mfull[a], w[a], m.cp_uinv[a],
                         uefull[a], m.cp_st[a].t, dsfull[a]);
        double cmax = 0.0;
        for (int a = 0; a < F; ++a) for (int b = 0; b < F; ++b)
            cmax = std::max(cmax, std::fabs(m.cp_C(a, b)));
        for (int a : {1, F / 2, F - 1}) {
            double crow = 0.0, cabs = 0.0;
            for (int b = 0; b < F; ++b) { crow += m.cp_C(a, b) * mfull[b];
                                          cabs = std::max(cabs, std::fabs(m.cp_C(a, b))); }
            std::fprintf(stderr, "[init] a=%d pan=%d uinv=%.4f Cm=%.4f ue=%.4f mrow_absmax=%.3f m=%.5f\n",
                         a, m.cp_pan[a], m.cp_uinv[a], crow, uefull[a], cabs, mfull[a]);
        }
        std::fprintf(stderr, "[init] global max|C|=%.4f  F=%d Nu=%d\n", cmax, F, Nu);

        // ---- Mresp conditioning probe (separate pathology from 1/ds physics) ----
        // Apply blowing fields that are SMOOTH IN PHYSICAL (panel-arc) space,
        // directly through the raw Mresp (panel-indexed), so the station
        // reordering does not introduce a spurious discontinuity.
        const double PI = 3.14159265358979;
        double mr_glob = 0.0;
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
            mr_glob = std::max(mr_glob, std::fabs(m.Mresp(i, j)));
        double Ltot = m.sarc[N - 1] + 0.5 * m.len[N - 1];
        double rz_mn = 1e9, rz_mx = -1e9, rn_mn = 1e9, rn_mx = -1e9, netz = 0.0;
        int iz_mx = 0, in_mx = 0;
        for (int i = 0; i < N; ++i) {
            double duz = 0.0, dun = 0.0;
            for (int j = 0; j < N; ++j) {
                double wz = std::sin(2.0 * PI * m.sarc[j] / Ltot); // smooth in arc
                duz += m.Mresp(i, j) * wz;
                dun += m.Mresp(i, j) * 1.0;                        // uniform
            }
            rz_mn = std::min(rz_mn, duz);
            if (duz > rz_mx) { rz_mx = duz; iz_mx = i; }
            if (std::fabs(dun) > std::fabs(rn_mn) || in_mx == 0) { in_mx = i; }
            rn_mx = std::max(rn_mx, dun);
            rn_mn = std::min(rn_mn, dun);
        }
        std::fprintf(stderr, "[mrsp] netzero spike at panel i=%d (TE pans are 0 and %d); uniform |spike| at i=%d\n",
                     iz_mx, N - 1, in_mx);
        for (int j = 0; j < N; ++j) netz += std::sin(2.0 * PI * m.sarc[j] / Ltot);
        std::fprintf(stderr, "[mrsp] max|Mresp|=%.3f\n", mr_glob);
        std::fprintf(stderr, "[mrsp] dUe(arc-smooth net-zero w=sin) range [%.4f,%.4f] (sum w=%.3f)\n",
                     rz_mn, rz_mx, netz);
        std::fprintf(stderr, "[mrsp] dUe(uniform w=1) range [%.4f,%.4f]\n", rn_mn, rn_mx);
    }

    // ---- 5) Newton iteration (VISCAL driver), wrapped in lam-continuation ----
    auto Jblk = m.cp_J.topLeftCorner(NX, NX);
    auto Rv   = m.cp_R.head(NX);
    auto Dv   = m.cp_dx.head(NX);
    bool converged = false;
    const double DHI = 1.5, DLO = -0.5;

    // continuation schedule for the displacement coupling strength
    std::vector<double> lam_sched = {0.15, 0.30, 0.50, 0.70, 0.85, 1.0};
    if (const char* e = std::getenv("XF_CONT")) {
        int ns = std::atoi(e);
        if (ns == 1) lam_sched = {1.0};                 // no continuation (debug)
        else if (ns > 1) { lam_sched.clear();
            for (int i = 1; i <= ns; ++i) lam_sched.push_back(double(i) / ns); }
    }

    for (size_t ls = 0; ls < lam_sched.size(); ++ls) {
        lam = lam_sched[ls];
        refresh_state();                                // re-seed Ue at this lam
        bool lam_conv = false;
    for (int iter = 0; iter < opt.max_newton; ++iter) {
        Jblk.setZero();
        Rv.setZero();

        for (int a = 0; a < F; ++a) {
            int g = m.cp_unk[a];
            if (g < 0) continue;                    // frozen stagnation station
            int pa = m.cp_prev[a];
            // Determine lam/turb from the transition index, NOT cp_st[a].turb:
            // the marcher's transition re-solve flags the last laminar (bubble-
            // peak) station turbulent while leaving its laminar H, which would
            // mislabel that interval's residual. itranF[side] is the first
            // turbulent full-station; <0 means the side stayed laminar.
            int itr = itranF[m.cp_side[a]];
            int ityp = (itr >= 0 && a >= itr) ? 2 : 1;

            LinBlk L;
            linearize(m.cp_st[pa], m.cp_st[a], ityp, Re, acrit, L);

            double uea = uefull[a], dsa = dsfull[a];
            double uep = uefull[pa], dsp = dsfull[pa];
            bool prev_unk = (m.cp_unk[pa] >= 0);
            int row = 3 * g;
            for (int i = 0; i < 3; ++i) Rv(row + i) = -L.R0[i];

            // station a's own s, theta columns
            for (int i = 0; i < 3; ++i) {
                Jblk(row + i, 3 * g + 0) += L.Ps2[i];
                Jblk(row + i, 3 * g + 1) += L.Pt2[i];
            }
            // upstream station's s, theta columns (if it is an unknown)
            if (prev_unk) {
                int gp = m.cp_unk[pa];
                for (int i = 0; i < 3; ++i) {
                    Jblk(row + i, 3 * gp + 0) += L.Ps1[i];
                    Jblk(row + i, 3 * gp + 1) += L.Pt1[i];
                }
            }
            // dense m-coupling: chain dUe/dm = C, ddstar/dm through both stations
            for (int c = 0; c < F; ++c) {
                int gc = m.cp_unk[c];
                if (gc < 0) continue;
                double Cac = lam * m.cp_C(a, c);
                double dd2 = ((c == a) ? 1.0 / uea : 0.0) - (dsa / uea) * Cac;
                double coef_d2 = dd2, coef_u2 = Cac;
                double coef_d1 = 0.0, coef_u1 = 0.0;
                if (prev_unk) {
                    double Cpc = lam * m.cp_C(pa, c);
                    coef_u1 = Cpc;
                    coef_d1 = ((c == pa) ? 1.0 / uep : 0.0) - (dsp / uep) * Cpc;
                }
                int mc = 3 * gc + 2;
                for (int i = 0; i < 3; ++i)
                    Jblk(row + i, mc) += L.Pd2[i] * coef_d2 + L.Pu2[i] * coef_u2
                                       + L.Pd1[i] * coef_d1 + L.Pu1[i] * coef_u1;
            }
        }

        if (!Rv.allFinite()) { r.converged = false; return r; }
        double rnorm = Rv.norm();
        if (std::getenv("XF_RDUMP") && iter == 0 && ls == 0) {
            // localize the largest interval residuals at the weak-march start
            struct RR { double v; int a, eq; };
            std::vector<RR> rr;
            for (int a = 0; a < F; ++a) { int g = m.cp_unk[a]; if (g < 0) continue;
                for (int eq = 0; eq < 3; ++eq) rr.push_back({std::fabs(Rv(3*g+eq)), a, eq}); }
            std::sort(rr.begin(), rr.end(), [](const RR&x,const RR&y){return x.v>y.v;});
            std::fprintf(stderr, "[rdump] top interval residuals (|R|=%.2e):\n", rnorm);
            for (int t = 0; t < 8 && t < (int)rr.size(); ++t) {
                int a = rr[t].a;
                std::fprintf(stderr, "[rdump]  R=%.3e eq=%d a=%d side=%d lpos=%d turb=%d x=%.5f ue=%.4f th=%.6f ds=%.6f H=%.3f\n",
                    rr[t].v, rr[t].eq, a, m.cp_side[a], m.cp_lpos[a], (int)m.cp_st[a].turb,
                    m.cp_st[a].x, uefull[a], m.cp_st[a].t, dsfull[a],
                    dsfull[a]/std::max(m.cp_st[a].t,1e-9));
            }
        }
        m.cp_lu.compute(Jblk);
        Dv.noalias() = m.cp_lu.solve(Rv);
        if (!Dv.allFinite()) { r.converged = false; return r; }
        if (std::getenv("XF_NEWT") && iter < 3)
            std::fprintf(stderr, "      |R|=%.3e  (lam=%.2f it=%d)\n", rnorm, lam, iter);

        // ---- UPDATE: trust-region under-relaxation (XFOIL UPDATE) ----
        double rlx = 1.0;
        double dnmax = 0.0; int amax = -1; char vmax = '?';
        for (int a = 0; a < F; ++a) {
            int g = m.cp_unk[a];
            if (g < 0) continue;
            double ds_var = Dv(3 * g + 0), dth = Dv(3 * g + 1), dm = Dv(3 * g + 2);
            // dUe from the mass-defect change at this station's row of C
            double due = 0.0;
            for (int c = 0; c < F; ++c) {
                int gc = m.cp_unk[c];
                if (gc >= 0) due += lam * m.cp_C(a, c) * Dv(3 * gc + 2);
            }
            double ddstr = (dm - dsfull[a] * due) / uefull[a];
            double dn1 = m.cp_st[a].turb ? ds_var / std::max(m.cp_st[a].s, 1e-4)
                                         : ds_var / 10.0;
            double dn2 = dth / std::max(m.cp_st[a].t, 1e-9);
            double dn3 = ddstr / std::max(dsfull[a], 1e-9);
            double dn4 = std::fabs(due) / 0.25;
            auto clamp = [&](double dn, char v) {
                if (std::fabs(dn) > dnmax) { dnmax = std::fabs(dn); amax = a; vmax = v; }
                if (rlx * dn > DHI) rlx = DHI / dn;
                if (rlx * dn < DLO) rlx = DLO / dn;
            };
            clamp(dn1, 'c'); clamp(dn2, 't'); clamp(dn3, 'd'); (void)dn4;
        }
        if (std::getenv("XF_NEWT") && amax >= 0)
            std::fprintf(stderr, "      dnmax=%.2e v=%c a=%d (lpos=%d ue=%.4f th=%.5f ds=%.5f turb=%d)\n",
                         dnmax, vmax, amax, m.cp_lpos[amax], uefull[amax],
                         m.cp_st[amax].t, dsfull[amax], (int)m.cp_st[amax].turb);

        // apply the under-relaxed step, accumulate rms change
        double rms = 0.0;
        int cnt = 0;
        for (int a = 0; a < F; ++a) {
            int g = m.cp_unk[a];
            if (g < 0) continue;
            double ds_var = rlx * Dv(3 * g + 0);
            double dth    = rlx * Dv(3 * g + 1);
            double dm     = rlx * Dv(3 * g + 2);
            m.cp_st[a].s += ds_var;
            m.cp_st[a].t += dth;
            mfull[a]     += dm;
            if (m.cp_st[a].t < 1e-7) m.cp_st[a].t = 1e-7;
            if (mfull[a]   < 1e-9)   mfull[a]   = 1e-9;
            if (m.cp_st[a].turb) m.cp_st[a].s = std::min(std::max(m.cp_st[a].s, 1e-7), 0.3);
            double dn2 = dth / std::max(m.cp_st[a].t, 1e-9);
            rms += dn2 * dn2;
            ++cnt;
        }
        rms = (cnt > 0) ? std::sqrt(rms / cnt) : 1e9;
        refresh_state();

        if (std::getenv("XF_NEWT"))
            std::fprintf(stderr, "[newt] lam=%.2f it=%2d rlx=%.3f rms=%.2e Nu=%d\n",
                         lam, iter, rlx, rms, Nu);

        if (rms < 1e-4) { lam_conv = true; break; }
    }
        if (!lam_conv) {            // failed to converge at this continuation step
            if (std::getenv("XF_NEWT"))
                std::fprintf(stderr, "[cont] stalled at lam=%.2f\n", lam);
            converged = false; break;
        }
        converged = (lam >= 1.0 - 1e-9);
    }

    // watchdog: a blown-up shape factor means a real (turbulent) separation
    if (converged)
        for (int a = 0; a < F; ++a) {
            if (!m.cp_st[a].turb) continue;
            double hk = m.cp_st[a].d / m.cp_st[a].t;
            if (!std::isfinite(hk) || hk > opt.hk_sep_limit) { r.separated = true; }
        }

    if (!converged) { r.converged = false; return r; }

    // ---- 6) forces from the coupled (decambered) edge velocity ----
    std::vector<double> ue_sgn(N);
    for (int i = 0; i < N; ++i) ue_sgn[i] = ue_inv[i];   // default (any gap)
    for (int a = 0; a < F; ++a)
        ue_sgn[m.cp_pan[a]] = m.cp_vti[a] * uefull[a];
    forces_from_cp(m, ue_sgn, alpha, r.cl, r.cm);

    // ---- far-wake Squire-Young drag from the corrected TE state ----
    auto te_state = [&](int is, double& th, double& ds, double& ue) {
        int a = -1;
        for (int b = 0; b < F; ++b) if (m.cp_side[b] == is) a = b;  // last on side
        th = m.cp_st[a].t; ue = uefull[a]; ds = dsfull[a];
    };
    double cd = 0.0;
    for (int is = 0; is < 2; ++is) {
        double th, ds, ue;
        te_state(is, th, ds, ue);
        double h = ds / th;
        if (th > 0.0 && ue > 0.0) cd += 2.0 * th * std::pow(ue, 0.5 * (5.0 + h));
    }
    r.cd  = cd;
    r.cdf = 0.0;                       // (friction split not recomputed here)
    r.cdp = cd;
    r.xtr_top = s0[0].xtr;
    r.xtr_bot = s0[1].xtr;
    r.converged = std::isfinite(cd) && cd > 0.0 && cd < opt.cd_sane_max
                  && !r.separated;
    return r;
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

    if (m.opt.strong_coupling) {
        // ---- strong coupling: the simultaneous implicit Newton (coupled_solve) ----
        return coupled_solve(m, alpha, Re);
    }

    // ---- weak coupling: one BL pass on the inviscid edge velocity ----
    if (!bl_pass(m, ue_inv, Re, dstar, sA, sB, idxA, idxB)) {
        r.converged = false; return r;
    }
    forces_from_cp(m, ue_inv, alpha, r.cl, r.cm);   // inviscid section loads

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
