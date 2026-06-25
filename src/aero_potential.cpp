// aero_potential.cpp — Milestone 3: Horseshoe Vortex Lattice Method (VLM)
//
// Architecture overview:
//   • solve() produces AeroState for a given (alpha, delta_e).
//   • The AIC matrix depends only on geometry, so it is factored once per
//     unique WingGeometry address and cached in a thread_local struct.
//     trim() calls solve() ~100× with the same &w — only the first call pays
//     the factorisation cost.
//   • delta_e remains modelled via the existing linear pitch-control
//     derivative heuristic (plan §5); see below for rationale.
//   • CL, x_np, and strip cl_local come from the VLM; CDi = CL²/(π·e·AR);
//     CDp and viscous corrections are unchanged from Phase 1.
//
// VLM layout (per spanwise strip i, i = 0..N-1):
//   bound vortex  → x_bv_i  = x_le_i + 0.25·c_i
//   collocation   → x_cp_i  = x_le_i + 0.75·c_i
//   strip edges   → y_{i,L} = y_i − width_i/2,  y_{i,R} = y_i + width_i/2
//
// Boundary condition (Dirichlet on downwash):
//   Σ_j AIC[i,j] · Γ_j = -(α + twist_i − α_L0)   ← enforced at 3c/4
//
// Lift coefficient:
//   CL = 4 · Σ_i (Γ_i · dy_i) / S_ref              ← Kutta–Joukowski, both halves

#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/control.h"
#include "aeroanalyzer/geom.h"
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstring>

namespace aero {
namespace potential {

// ---- Public analytic helpers (tests call these directly) -----------------

double lift_curve_slope_3d(double a0, double AR, double e) {
    if (AR <= 0.0) return 0.0;
    return a0 / (1.0 + a0 / (PI * e * AR));
}

double oswald_e(double AR, double taper) {
    double e = 0.95 - 0.6 * (taper - 0.40) * (taper - 0.40);
    return std::max(0.60, std::min(0.98, e));
}

// Planform-area-weighted aerodynamic centre (c/4 locus).
double wing_ac_x(const WingGeometry& w) {
    double num = 0.0, den = 0.0;
    for (const auto& s : w.stations) {
        double area = s.chord * s.width;
        num += (s.x_le + 0.25 * s.chord) * area;
        den += area;
    }
    return (den > 0) ? num / den : 0.25 * w.root_chord;
}

namespace {

// ---- VLM influence kernel ------------------------------------------------
//
// z-velocity (positive = up) induced at (x_cp, y_cp, z=0)
// by a horseshoe vortex of unit strength Γ=1 whose:
//   • bound vortex runs along x = x_bv from y = y_L to y = y_R  (in +y direction)
//   • right trailing vortex exits at y_R → +∞x  (strength +1)
//   • left  trailing vortex exits at y_L → +∞x  (strength −1)
//
// Derivation: Biot-Savart for finite + semi-infinite filaments.
// All three contributions have the same sign (downwash) for a lifting
// horseshoe when the collocation is downstream of the bound vortex.
static constexpr double FOUR_PI = 4.0 * PI;
static constexpr double VLM_EPS = 1.0e-9;  // singularity guard

double vlm_w(double x_cp, double y_cp,
             double x_bv, double y_L, double y_R) {
    const double dx  = x_cp - x_bv;
    const double dyL = y_cp - y_L;   // > 0 when collocation is to the right of y_L
    const double dyR = y_cp - y_R;   // < 0 when collocation is between y_L and y_R

    const double r1 = std::sqrt(dx * dx + dyL * dyL) + VLM_EPS;
    const double r2 = std::sqrt(dx * dx + dyR * dyR) + VLM_EPS;

    // Bound vortex (finite segment)
    double w_bv = 0.0;
    if (std::fabs(dx) > VLM_EPS) {
        w_bv = -1.0 / (FOUR_PI * dx) * (dyL / r1 - dyR / r2);
    }

    // Right trailing vortex at y_R (strength +1, goes to +∞)
    double w_tr = 0.0;
    if (std::fabs(dyR) > VLM_EPS) {
        w_tr += 1.0 / (FOUR_PI * dyR) * (1.0 + dx / r2);
    }

    // Left trailing vortex at y_L (strength −1, goes to +∞)
    if (std::fabs(dyL) > VLM_EPS) {
        w_tr -= 1.0 / (FOUR_PI * dyL) * (1.0 + dx / r1);
    }

    return w_bv + w_tr;
}

// ---- VLM strip descriptor ------------------------------------------------
struct VLMStrip {
    double x_bv, x_cp;  // bound vortex and collocation x-coords
    double y_c;          // strip centre y
    double y_L, y_R;     // strip edges (y_L < y_R)
    double dy;           // strip width
    double chord;
    double twist;        // local geometric twist, rad
};

// ---- Geometry fingerprint ------------------------------------------------
// The AIC factorisation may only be reused when the geometry is *identical*.
// Keying the cache on &w is unsound: a new WingGeometry can reuse a previous
// one's address (e.g. the per-candidate EvalResult stack local in the GA loop),
// which silently served a stale factorisation. We key on the content instead:
// a 64-bit FNV-1a hash over every field that enters the AIC. Within a single
// trim() the geometry is constant across the ~100 solve() calls, so the
// factorisation is still built once; across candidates the hash differs and the
// AIC is correctly rebuilt.
std::uint64_t geom_signature(const WingGeometry& w) {
    std::uint64_t h = 1469598103934665603ull;  // FNV-1a offset basis
    auto mix = [&h](double d) {
        std::uint64_t b;
        std::memcpy(&b, &d, sizeof(b));
        h = (h ^ b) * 1099511628211ull;        // FNV-1a prime
    };
    mix(w.root_chord); mix(w.tip_chord); mix(w.semi_span);
    mix(w.le_sweep); mix(w.washout); mix(w.le_bow); mix(w.te_bow);
    for (const auto& sec : w.sections) {
        mix(sec.te_thick);
        for (double v : sec.wu) mix(v);
        for (double v : sec.wl) mix(v);
    }
    mix(static_cast<double>(w.stations.size()));
    for (const auto& s : w.stations) { mix(s.y); mix(s.chord); mix(s.x_le); mix(s.twist); }
    return h;
}

// ---- Thread-local AIC factorisation cache --------------------------------
struct VLMCache {
    std::uint64_t sig = 0;        // geometry fingerprint of the factored AIC
    bool valid = false;          // sig/fact/strips populated for a real geometry
    Eigen::PartialPivLU<Eigen::MatrixXd> fact;   // factored AIC, reused per RHS
    bool factored = false;
    std::vector<VLMStrip> strips;
};

thread_local VLMCache t_cache;

void build_vlm(const WingGeometry& w) {
    t_cache.factored = false;
    t_cache.valid = false;
    const int N = static_cast<int>(w.stations.size()) - 1;
    if (N < 1) return;

    t_cache.strips.resize(N);
    for (int i = 0; i < N; ++i) {
        const Station& si = w.stations[i];
        const Station& si1 = w.stations[i + 1];
        VLMStrip& st = t_cache.strips[i];

        // Strip centroid between adjacent stations
        double y_c = 0.5 * (si.y + si1.y);
        double c_c = 0.5 * (si.chord + si1.chord);
        double x_le_c = 0.5 * (si.x_le + si1.x_le);
        double tw_c = 0.5 * (si.twist + si1.twist);

        st.y_c  = y_c;
        st.y_L  = si.y;
        st.y_R  = si1.y;
        st.dy   = si1.y - si.y;
        st.chord = c_c;
        st.twist = tw_c;
        st.x_bv  = x_le_c + 0.25 * c_c;
        st.x_cp  = x_le_c + 0.75 * c_c;
    }

    // Assemble AIC: A[i,j] = downwash at strip i from unit horseshoe j
    // (right half) + mirror horseshoe j (left half)
    Eigen::MatrixXd A(N, N);
    for (int i = 0; i < N; ++i) {
        double xc = t_cache.strips[i].x_cp;
        double yc = t_cache.strips[i].y_c;
        for (int j = 0; j < N; ++j) {
            const VLMStrip& sj = t_cache.strips[j];
            // Right half-wing horseshoe
            double w_r = vlm_w(xc, yc, sj.x_bv, sj.y_L, sj.y_R);
            // Mirror (left half-wing) horseshoe — edges reflected about y=0.
            // The mirror bound vortex runs from −y_R to −y_L (same +y sense,
            // same circulation sense for symmetric lift), so swap sign of edges:
            double w_m = vlm_w(xc, yc, sj.x_bv, -sj.y_R, -sj.y_L);
            A(i, j) = w_r + w_m;
        }
    }

    t_cache.fact.compute(A);   // partial-pivot LU, reused across trim RHS
    t_cache.factored = true;
    t_cache.sig = geom_signature(w);
    t_cache.valid = true;
}

}  // anonymous namespace

// ---- Main aerodynamic solve ----------------------------------------------

AeroState solve(const WingGeometry& w, const MassProps& mp,
                const viscous::Surrogate& surr, const Config& cfg,
                double alpha, double delta_e) {
    // Milestone 3: the Morino panel solver is the default model of record
    // (validated vs AVL within a few percent on the knee/min_drag/min_mass decks).
    // The VLM below is the documented analytic fallback (aero_model = vlm).
    if (cfg.gets("aero_model", "panel") != "vlm")
        return panel::solve(w, mp, surr, cfg, alpha, delta_e);

    AeroState st;
    st.alpha   = alpha;
    st.delta_e = delta_e;

    // Section aerodynamics from thin-airfoil theory (camber line integration)
    geom::ThinAirfoil ta = geom::thin_airfoil(w.sections.empty() ? Airfoil{} : w.sections[0]);

    // Prandtl 3-D correction for control-derivative reference slope
    const double AR    = mp.AR;
    const double taper = (w.root_chord > 0) ? w.tip_chord / w.root_chord : 1.0;
    const double e_ref = oswald_e(AR, taper);
    st.e = e_ref;
    const double a_ref = lift_curve_slope_3d(ta.cl_alpha, AR, e_ref);

    control::Derivs pc = control::compute(w, mp, a_ref, alpha, delta_e, cfg);

    // ---- Build / reuse VLM factorisation --------------------------------
    // Reuse only when the cached factorisation belongs to this exact geometry
    // (content hash, NOT pointer identity — see geom_signature).
    if (!t_cache.valid || !t_cache.factored || t_cache.strips.empty() ||
        t_cache.sig != geom_signature(w)) {
        build_vlm(w);
    }

    const int N = static_cast<int>(t_cache.strips.size());

    // ---- Solve A·Γ = rhs (reuse factored AIC) ---------------------------
    // RHS: -(α_eff per strip) where α_eff = α + twist − α_L0
    std::vector<double> Gamma(N, 0.0);
    bool ok = t_cache.factored && N > 0;
    if (ok) {
        Eigen::VectorXd rhs(N);
        for (int i = 0; i < N; ++i)
            rhs(i) = -(alpha + t_cache.strips[i].twist - ta.alpha_L0);
        Eigen::VectorXd G = t_cache.fact.solve(rhs);
        for (int i = 0; i < N; ++i) Gamma[i] = G(i);
    }

    // ---- Lift & neutral-point from Kutta-Joukowski ----------------------
    // CL = (4 · Σ Γ_i·dy_i) / S_ref  [both half-wings, normalised V=1]
    double sum_G  = 0.0;   // Σ Γ·dy
    double sum_Gx = 0.0;   // Σ Γ·dy·x_cp  (for x_np)
    if (ok) {
        for (int i = 0; i < N; ++i) {
            double Gdy = Gamma[i] * t_cache.strips[i].dy;
            sum_G  += Gdy;
            // Neutral point lives at the 0.25c bound-vortex locus (x_bv), NOT
            // the 0.75c collocation point (x_cp). Weighting by x_cp put x_np a
            // half-chord too far aft, corrupting every static_margin.
            sum_Gx += Gdy * t_cache.strips[i].x_bv;
        }
    }

    double S_ref = mp.S_ref > 0 ? mp.S_ref : 1.0;
    double panel_CL = (ok && sum_G > 0.0) ? 4.0 * sum_G / S_ref : 0.0;

    // Neutral point (lift-weighted collocation x, falls back to planform AC)
    double x_np = (panel_CL > 1e-12) ? sum_Gx / sum_G : wing_ac_x(w);
    st.x_np     = x_np;
    st.static_margin = (mp.mac > 0) ? (x_np - mp.x_cg) / mp.mac : 0.0;

    // Total CL including control surface
    st.CL = panel_CL + pc.CLde * delta_e;

    // ---- Induced drag (Trefftz-plane result for span efficiency e_ref) --
    st.CDi = (AR > 0) ? st.CL * st.CL / (PI * e_ref * AR) : 0.0;

    // ---- Strip viscous coupling (sweep-normal) --------------------------
    const double V      = cfg.getd("v_cruise", V_CRUISE);
    const double cosL   = std::cos(w.le_sweep);
    const double sweep_deg = w.le_sweep * RAD2DEG;
    const double xflow_deg = cfg.getd("sweep_crossflow_deg", 25.0);
    const double xflow_fac = cfg.getd("crossflow_factor", 1.15);
    const double Shalf  = 0.5 * S_ref;

    std::vector<double> shape;  // rebuilt per station from the lofted section

    st.cl_local.assign(w.stations.size(), 0.0);
    double CDp_num  = 0.0;
    bool   tip_stall = false;
    double min_conf = 1.0;

    // Induced-angle estimate for strip cl (elliptic approximation for the
    // viscous query; the VLM already gave us the global CL accurately).
    double a_ind = (AR > 0) ? st.CL / (PI * e_ref * AR) : 0.0;

    for (std::size_t si = 0; si < w.stations.size(); ++si) {
        const Station& s = w.stations[si];
        double t = (w.semi_span > 0) ? s.y / w.semi_span : 0.0;

        // Local section cl: VLM per-strip where available, else lifting-line
        double cl_i;
        if (ok && si < static_cast<std::size_t>(N)) {
            // Strip cl from Kutta-Joukowski: cl = 2Γ/(V·c), V=1
            cl_i = (s.chord > 0) ? 2.0 * Gamma[si] / s.chord : 0.0;
        } else {
            cl_i = ta.cl_alpha * (alpha + s.twist - ta.alpha_L0 - a_ind);
        }
        st.cl_local[si] = cl_i;

        double Re       = RHO * V * cosL * s.chord / MU;
        double cl_norm  = cl_i / (cosL * cosL);
        shape.clear();
        shape.insert(shape.end(), s.af.wu.begin(), s.af.wu.end());
        shape.insert(shape.end(), s.af.wl.begin(), s.af.wl.end());
        viscous::Polar pol = surr.query(shape, cl_norm, Re, s.af.te_thick);
        min_conf = std::min(min_conf, pol.confidence);

        // Crossflow penalty on outer span for swept leading edges
        double ramp = 0.0;
        if (sweep_deg > xflow_deg - 3.0) {
            double u = (sweep_deg - (xflow_deg - 3.0)) / 6.0;
            ramp = std::max(0.0, std::min(1.0, u));
        }
        double outer  = std::max(0.0, (t - 0.5) / 0.5);
        double factor = 1.0 + (xflow_fac - 1.0) * ramp * outer;

        CDp_num += pol.cd * factor * s.chord * s.width;

        if (pol.clamped && cl_i > 0.0 && t > 0.6) tip_stall = true;
    }
    st.CDp    = (Shalf > 0) ? CDp_num / Shalf : 0.0;
    st.CD     = st.CDi + st.CDp;
    st.tip_stall = tip_stall;
    st.polar_confidence = min_conf;
    st.cn_da  = control::adverse_yaw_cn_da(w, mp, surr, st.cl_local, a_ref, cfg);

    // ---- Pitching moment about CG ----------------------------------------
    // CM = CM_0 (section camber) + CL·(x_cg − x_np)/mac + Cm_δe + Cm_thrust
    double Cm0  = ta.cm_ac;
    double Cm_cl = st.CL * ((mp.x_cg - x_np) / (mp.mac > 0 ? mp.mac : 1.0));
    double Cm_de = pc.Cmde * delta_e;

    double q      = 0.5 * RHO * V * V;
    double T      = st.CD * q * S_ref;
    double zt     = cfg.getd("thrust_z_offset", 0.020);
    double denom  = q * S_ref * (mp.mac > 0 ? mp.mac : 1.0);
    double Cm_thr = (denom > 0) ? (T * zt) / denom : 0.0;

    st.CM = Cm0 + Cm_cl + Cm_de - Cm_thr;

    // ---- Hinge moment (Glauert thin-airfoil, worst-case pitch+roll) -----
    st.hinge_moment = pc.hinge_moment;

    // ---- Roll control / damping / helix (M6) ----------------------------
    st.cl_da      = pc.Cl_da;
    st.cl_p       = pc.Cl_p;
    st.roll_helix = pc.roll_helix;

    // ---- High-alpha NP migration heuristic (VLM fallback; no strip cl_max) -
    double tip_unload     = 0.12 * std::sin(w.le_sweep);
    double washout_relief = std::max(0.0, -w.washout) * 3.0;
    double fwd = std::max(0.0, tip_unload - washout_relief);
    st.x_np_high = x_np - fwd * (mp.mac > 0 ? mp.mac : 0.0);

    return st;
}

}  // namespace potential
}  // namespace aero
