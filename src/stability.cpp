#include "aeroanalyzer/stability.h"
#include "aeroanalyzer/aero_potential.h"
#include <cmath>
#include <algorithm>

namespace aero {
namespace stability {

double sm_objective(double sm, double lo, double hi) {
    if (sm < lo) return lo - sm;
    if (sm > hi) return sm - hi;
    return 0.0;
}

AeroState trim(const WingGeometry& w, const MassProps& mp,
               const viscous::Surrogate& surr, const Config& cfg) {
    const double V = cfg.getd("v_cruise", V_CRUISE);
    const double q = 0.5 * RHO * V * V;
    const double W = mp.mass * GRAV;
    const double CL_req = (q * mp.S_ref > 0) ? W / (q * mp.S_ref) : 0.0;

    double alpha = 2.0 * DEG2RAD;   // initial guess
    double delta = 0.0;

    AeroState st = potential::solve(w, mp, surr, cfg, alpha, delta);
    const int max_iter = 50;
    const double tol = 1e-8;
    const double h = 1e-6;          // FD step (rad)

    for (int it = 0; it < max_iter; ++it) {
        st = potential::solve(w, mp, surr, cfg, alpha, delta);
        double r0 = st.CL - CL_req;
        double r1 = st.CM;
        if (std::fabs(r0) < tol && std::fabs(r1) < tol) {
            st.trimmed = true;
            break;
        }
        // finite-difference Jacobian J = d[r]/d[alpha,delta]
        AeroState sa = potential::solve(w, mp, surr, cfg, alpha + h, delta);
        AeroState sd = potential::solve(w, mp, surr, cfg, alpha, delta + h);
        double j00 = (sa.CL - st.CL) / h;   // dCL/dalpha
        double j01 = (sd.CL - st.CL) / h;   // dCL/ddelta
        double j10 = (sa.CM - st.CM) / h;   // dCm/dalpha
        double j11 = (sd.CM - st.CM) / h;   // dCm/ddelta
        double det = j00 * j11 - j01 * j10;
        if (std::fabs(det) < 1e-12) break;  // singular -> leave untrimmed
        double da = (-r0 * j11 + r1 * j01) / det;
        double dd = (-j00 * r1 + j10 * r0) / det;
        // step clamp for robustness against the nonlinear panel model later
        double maxstep = 5.0 * DEG2RAD;
        da = std::max(-maxstep, std::min(maxstep, da));
        dd = std::max(-maxstep, std::min(maxstep, dd));
        alpha += da;
        delta += dd;
    }
    st = potential::solve(w, mp, surr, cfg, alpha, delta);
    st.trimmed = st.trimmed ||
        (std::fabs(st.CL - CL_req) < 1e-5 && std::fabs(st.CM) < 1e-5);
    return st;
}

}  // namespace stability
}  // namespace aero
