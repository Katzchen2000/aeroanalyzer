#include "aeroanalyzer/geom.h"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

namespace aero {
namespace geom {

namespace {
// Local dihedral angle above this threshold marks a station "in_winglet" for
// the relaxed chord floor (evaluate.cpp constraint 14). Organic tips have no
// discrete winglet device, so this is just a steep-curve marker, not a gene.
constexpr double WINGLET_DIH_THRESHOLD_DEG = 60.0;
}  // namespace

GenomeSpec default_genome(const Config& cfg) {
    (void)cfg;  // no toggles left; every gene is unconditional
    GenomeSpec g;
    g.lo.resize(N_GENES);
    g.hi.resize(N_GENES);
    g.names.resize(N_GENES);
    auto set = [&](int i, const char* nm, double lo, double hi) {
        g.names[i] = nm; g.lo[i] = lo; g.hi[i] = hi;
    };
    set(G_SEMISPAN, "semi_span_m", 0.45, 0.80);

    // Chord curve: CP0 = root chord (m); CP1..CP5 free organic control points (m).
    // Positive control points keep the curve inside the convex hull -> chord
    // stays positive everywhere without a monotone straitjacket.
    set(G_CHORD_CP + 0, "chord_cp0_root_m", 0.18, 0.35);
    for (int i = 1; i < NCP_CHORD; ++i) {
        char nm[24];
        std::snprintf(nm, sizeof(nm), "chord_cp%d_m", i);
        set(G_CHORD_CP + i, nm, 0.03, 0.35);
    }

    // Sweep curve: x_le(eta)/semi_span Bezier; CP0 pinned to 0 (root LE at x=0),
    // CP1..CP5 are genes. Small negative allows a gentle forward inner LE;
    // larger positive lets the curve sweep aft/crescent toward the tip.
    for (int i = 1; i < NCP_SWEEP; ++i) {
        char nm[24];
        std::snprintf(nm, sizeof(nm), "sweep_cp%d", i);
        set(G_SWEEP_CP + (i - 1), nm, -0.05, 0.80);
    }

    // Twist curve: fully free per control point, deg. CP0 = root incidence,
    // CP5 = tip; interior points give non-linear washout (the whole point).
    for (int i = 0; i < NCP_TWIST; ++i) {
        char nm[24];
        std::snprintf(nm, sizeof(nm), "twist_cp%d_deg", i);
        set(G_TWIST_CP + i, nm, -8.0, 4.0);
    }

    // Dihedral curve: local curve angle Bezier, deg. CP0 pinned to 0 (root
    // flat), CP1..CP6 are genes. Capped at 30 deg (was 85) -- near-vertical
    // folds fall outside the Prandtl-Munk nonplanar CDi correction's valid
    // regime (aero_panel.cpp) and the GA exploited that gap. decode() also
    // enforces monotonicity so the curve can't fold then unfold.
    for (int i = 1; i < NCP_DIH; ++i) {
        char nm[24];
        std::snprintf(nm, sizeof(nm), "dih_cp%d_deg", i);
        set(G_DIH_CP + (i - 1), nm, -10.0, 30.0);
    }

    set(G_BATTERY,  "battery_x_m",   0.00, 0.22);
    set(G_CS_CHORD, "cs_chord_frac", 0.15, 0.35);
    set(G_AIL_SPAN, "ail_span_frac", 0.40, 0.80);

    // CST bounds identical for every section; aft lower weight drives reflex.
    const double ulo[] = {0.10, 0.10, 0.08, 0.06};
    const double uhi[] = {0.34, 0.34, 0.30, 0.26};
    const double llo[] = {-0.22, -0.20, -0.15, -0.05};
    const double lhi[] = {0.05,   0.08,  0.15,  0.20};
    const char* pfx[]  = {"s0", "s1", "s2", "s3", "s4"};
    const char* usuf[] = {"_wu0","_wu1","_wu2","_wu3"};
    const char* lsuf[] = {"_wl0","_wl1","_wl2","_wl3"};
    for (int k = 0; k < N_SECTIONS; ++k) {
        for (int i = 0; i < 4; ++i) {
            char nm[16];
            std::snprintf(nm, sizeof(nm), "%s%s", pfx[k], usuf[i]);
            set(G_SEC(k, 0, i), nm, ulo[i], uhi[i]);
            std::snprintf(nm, sizeof(nm), "%s%s", pfx[k], lsuf[i]);
            set(G_SEC(k, 1, i), nm, llo[i], lhi[i]);
        }
    }
    return g;
}

double bernstein(int i, int n, double x) {
    // binomial(n,i) x^i (1-x)^(n-i)
    double c = 1.0;
    for (int k = 0; k < i; ++k) c = c * (n - k) / (k + 1);
    return c * std::pow(x, i) * std::pow(1.0 - x, n - i);
}

namespace {
// Generic Bezier evaluator over an arbitrary control-point vector.
double bezier_eval(const std::vector<double>& cp, double t) {
    int n = static_cast<int>(cp.size()) - 1;
    if (n < 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i <= n; ++i) s += cp[i] * bernstein(i, n, t);
    return s;
}
}  // namespace

double chord_at(const WingGeometry& w, double eta) {
    return bezier_eval(w.chord_cp, eta);
}
double xle_at(const WingGeometry& w, double eta) {
    return w.semi_span * bezier_eval(w.sweep_cp, eta);
}
double twist_at(const WingGeometry& w, double eta) {
    return bezier_eval(w.twist_cp, eta);
}
double dihedral_at(const WingGeometry& w, double eta) {
    return bezier_eval(w.dih_cp, eta);
}

static double class_fn(const Airfoil& f, double x) {
    return std::pow(x, f.N1) * std::pow(1.0 - x, f.N2);
}

static double shape(const std::vector<double>& w, double x) {
    int n = static_cast<int>(w.size()) - 1;
    if (n < 0) return 0.0;
    double s = 0.0;
    for (int i = 0; i <= n; ++i) s += w[i] * bernstein(i, n, x);
    return s;
}

double cst_upper(const Airfoil& f, double x) {
    return class_fn(f, x) * shape(f.wu, x) + x * 0.5 * f.te_thick;
}
double cst_lower(const Airfoil& f, double x) {
    return class_fn(f, x) * shape(f.wl, x) - x * 0.5 * f.te_thick;
}

std::vector<std::pair<double, double>> camber_line(const Airfoil& f, int n) {
    std::vector<std::pair<double, double>> c(n);
    for (int i = 0; i < n; ++i) {
        double theta = PI * i / (n - 1);
        double x = 0.5 * (1.0 - std::cos(theta));     // cosine clustering
        double zc = 0.5 * (cst_upper(f, x) + cst_lower(f, x));
        c[i] = {x, zc};
    }
    return c;
}

Airfoil fit_cst(const std::vector<std::pair<double, double>>& upper,
                const std::vector<std::pair<double, double>>& lower,
                int order, double te_thick) {
    Airfoil f;
    f.te_thick = te_thick;
    auto fit_side = [&](const std::vector<std::pair<double, double>>& pts,
                        double te_sign) -> std::vector<double> {
        const int n = order;                 // weights = n+1
        const Eigen::Index m = static_cast<Eigen::Index>(pts.size());
        Eigen::MatrixXd A(m, n + 1);
        Eigen::VectorXd b(m);
        for (Eigen::Index r = 0; r < m; ++r) {
            double x = pts[r].first;
            double cf = std::pow(x, f.N1) * std::pow(1.0 - x, f.N2);
            for (int i = 0; i <= n; ++i) A(r, i) = cf * bernstein(i, n, x);
            b(r) = pts[r].second - te_sign * x * 0.5 * te_thick;
        }
        // Robust least squares directly on the (rectangular) system — more
        // stable than normal equations for the Bernstein basis.
        Eigen::VectorXd w = A.colPivHouseholderQr().solve(b);
        return std::vector<double>(w.data(), w.data() + w.size());
    };
    f.wu = fit_side(upper, +1.0);
    f.wl = fit_side(lower, -1.0);
    return f;
}

ThinAirfoil thin_airfoil(const Airfoil& f) {
    // Numerically integrate camber slope over theta in [0,pi] with x=(1-cos)/2.
    const int N = 200;
    auto dzdx = [&](double x) {
        double h = 1e-4;
        double xa = std::max(1e-5, x - h);
        double xb = std::min(1.0 - 1e-5, x + h);
        double za = 0.5 * (cst_upper(f, xa) + cst_lower(f, xa));
        double zb = 0.5 * (cst_upper(f, xb) + cst_lower(f, xb));
        return (zb - za) / (xb - xa);
    };
    double I_alpha0 = 0.0, A1 = 0.0, A2 = 0.0;  // trapezoidal in theta
    double prev_g0 = 0, prev_g1 = 0, prev_g2 = 0, prev_th = 0;
    for (int i = 0; i < N; ++i) {
        double th = PI * i / (N - 1);
        double x = 0.5 * (1.0 - std::cos(th));
        double s = dzdx(x);
        double g0 = s * (1.0 - std::cos(th));   // -> alpha_L0 = (1/pi)∫ g0 dθ
        double g1 = s * std::cos(th);
        double g2 = s * std::cos(2.0 * th);
        if (i > 0) {
            double dt = th - prev_th;
            I_alpha0 += 0.5 * (g0 + prev_g0) * dt;
            A1       += 0.5 * (g1 + prev_g1) * dt;
            A2       += 0.5 * (g2 + prev_g2) * dt;
        }
        prev_g0 = g0; prev_g1 = g1; prev_g2 = g2; prev_th = th;
    }
    ThinAirfoil t;
    t.alpha_L0 = I_alpha0 / PI;
    A1 *= 2.0 / PI;
    A2 *= 2.0 / PI;
    t.cm_ac = 0.25 * PI * (A2 - A1);   // >0 for reflex
    t.cl_alpha = 2.0 * PI;
    return t;
}

namespace {
std::vector<double> linear_cp(int n, double v0, double v1) {
    std::vector<double> cp(n);
    for (int i = 0; i < n; ++i)
        cp[i] = v0 + (v1 - v0) * (n > 1 ? double(i) / (n - 1) : 0.0);
    return cp;
}
}  // namespace

void set_linear_planform(WingGeometry& w, double root_chord, double tip_chord,
                         double le_sweep_rad, double washout_rad) {
    w.chord_cp = linear_cp(NCP_CHORD, root_chord, tip_chord);
    double sweep_tip_frac = std::tan(le_sweep_rad);
    w.sweep_cp = linear_cp(NCP_SWEEP, 0.0, sweep_tip_frac);
    w.twist_cp = linear_cp(NCP_TWIST, 0.0, washout_rad);
    if (w.dih_cp.empty()) w.dih_cp.assign(NCP_DIH, 0.0);
    w.root_chord = root_chord;
    w.tip_chord  = tip_chord;
    w.le_sweep   = le_sweep_rad;
    w.washout    = washout_rad;
}

WingGeometry decode(const std::vector<double>& g, const GenomeSpec& spec, const Config& cfg) {
    (void)cfg;  // no toggles left; genes fully determine the wing
    auto clamp = [&](int i) {
        double v = g[i];
        if (v < spec.lo[i]) v = spec.lo[i];
        if (v > spec.hi[i]) v = spec.hi[i];
        return v;
    };
    WingGeometry w;
    w.semi_span = clamp(G_SEMISPAN);

    // ponytail: monotonic control points -> monotonic Bezier curve (convex-hull
    // property) -- cheapest way to forbid chord flare-back / sweep swing-forward
    // / dihedral fold-then-unfold without a separate penalty/constraint system.
    w.chord_cp.resize(NCP_CHORD);
    for (int i = 0; i < NCP_CHORD; ++i) {
        w.chord_cp[i] = clamp(G_CHORD_CP + i);
        if (i > 0) w.chord_cp[i] = std::min(w.chord_cp[i], w.chord_cp[i - 1]);  // non-increasing (taper)
    }

    w.sweep_cp.resize(NCP_SWEEP);
    w.sweep_cp[0] = 0.0;   // root LE pinned at x_le=0
    for (int i = 1; i < NCP_SWEEP; ++i) {
        w.sweep_cp[i] = clamp(G_SWEEP_CP + (i - 1));
        w.sweep_cp[i] = std::max(w.sweep_cp[i], w.sweep_cp[i - 1]);  // non-decreasing (no forward swing)
    }

    w.twist_cp.resize(NCP_TWIST);
    for (int i = 0; i < NCP_TWIST; ++i) w.twist_cp[i] = clamp(G_TWIST_CP + i) * DEG2RAD;

    w.dih_cp.resize(NCP_DIH);
    w.dih_cp[0] = 0.0;   // root dihedral pinned flat
    for (int i = 1; i < NCP_DIH; ++i) {
        w.dih_cp[i] = clamp(G_DIH_CP + (i - 1)) * DEG2RAD;
        w.dih_cp[i] = std::max(w.dih_cp[i], w.dih_cp[i - 1]);  // non-decreasing (no fold-then-unfold)
    }

    w.battery_x     = clamp(G_BATTERY);
    // ponytail: box-fit clamp moved to massprops::compute() where battery_len_m is available
    w.mode          = ControlMode::Elevon;  // G_MODE dropped; always elevon
    w.cs_chord_frac = clamp(G_CS_CHORD);
    w.ail_span_frac = clamp(G_AIL_SPAN);

    // Derived summaries for analytical downstream models (single source of
    // truth: computed from the same curves loft() will sample).
    w.root_chord = chord_at(w, 0.0);
    w.tip_chord  = chord_at(w, 1.0);
    w.le_sweep   = std::atan2(xle_at(w, 1.0), w.semi_span > 0.0 ? w.semi_span : 1.0);
    w.washout    = twist_at(w, 1.0) - twist_at(w, 0.0);
    // z_tip / nonplanar_h need the arc-integrated station z; loft() fills them.

    w.sections.resize(N_SECTIONS);
    for (int k = 0; k < N_SECTIONS; ++k) {
        w.sections[k].wu = {clamp(G_SEC(k,0,0)), clamp(G_SEC(k,0,1)),
                            clamp(G_SEC(k,0,2)), clamp(G_SEC(k,0,3))};
        w.sections[k].wl = {clamp(G_SEC(k,1,0)), clamp(G_SEC(k,1,1)),
                            clamp(G_SEC(k,1,2)), clamp(G_SEC(k,1,3))};
        w.sections[k].te_thick = 0.0;  // sharp TE everywhere; motor pocket via CAD split-plane
    }
    return w;
}

void loft(WingGeometry& w, int n) {
    if (n < 2) n = 2;
    w.stations.assign(n, Station{});
    if (w.sections.empty()) w.sections.resize(1);
    const int K = (int)w.sections.size();
    // η breakpoints: canonical 5-section table when K==5; linear spacing otherwise.
    // (CST loft breakpoints only — independent of the smooth planform curves.)
    std::vector<double> ETA(K);
    if (K == N_SECTIONS) {
        for (int k = 0; k < K; ++k) ETA[k] = SECTION_ETA[k];
    } else {
        for (int k = 0; k < K; ++k) ETA[k] = (K > 1) ? double(k) / (K - 1) : 0.0;
    }
    const double b = w.semi_span;
    // Cosine-spaced flat arc coords
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double th = PI * i / (n - 1);
        y[i] = b * 0.5 * (1.0 - std::cos(th));
    }
    // Pass 1: physical (yp, z_arr) and local dihedral angle via arc integration
    // of the smooth dihedral-angle Bezier curve. Continuous angle -> no crease;
    // a smoothly rising angle near the tip is an organic raised tip.
    std::vector<double> yp(n), z_arr(n), phi_arr(n);
    yp[0] = 0.0;
    z_arr[0] = 0.0;
    phi_arr[0] = dihedral_at(w, 0.0);
    for (int i = 1; i < n; ++i) {
        double d_arc = y[i] - y[i-1];
        double t_prev = (b > 0.0) ? y[i-1] / b : 0.0;
        double t = (b > 0.0) ? y[i] / b : 0.0;
        double t_mid = 0.5 * (t_prev + t);
        double phi_mid = dihedral_at(w, t_mid);
        yp[i] = yp[i-1] + d_arc * std::cos(phi_mid);
        z_arr[i] = z_arr[i-1] + d_arc * std::sin(phi_mid);
        phi_arr[i] = dihedral_at(w, t);
    }
    double max_abs_z = 0.0;
    for (int i = 0; i < n; ++i) max_abs_z = std::max(max_abs_z, std::fabs(z_arr[i]));
    w.nonplanar_h = max_abs_z;
    w.z_tip = (b > 0.0) ? z_arr[n - 1] / b : 0.0;

    for (int i = 0; i < n; ++i) {
        double t = (b > 0) ? y[i] / b : 0.0;
        Station s;
        s.y        = yp[i];
        s.z        = z_arr[i];
        s.dihedral = phi_arr[i];
        s.eta      = t;
        s.chord    = std::max(0.005, chord_at(w, t));   // safety floor only
        s.x_le     = xle_at(w, t);
        s.twist    = twist_at(w, t);
        s.in_winglet = (phi_arr[i] * RAD2DEG) > WINGLET_DIH_THRESHOLD_DEG;

        // piecewise linear blend between K control sections (CST airfoils only)
        int seg = K - 2;
        for (int k = 0; k < K - 1; ++k)
            if (ETA[k + 1] >= t) { seg = k; break; }
        if (seg < 0) seg = 0;
        if (K > 1 && seg > K - 2) seg = K - 2;
        double tl = (ETA[seg+1] > ETA[seg])
                    ? (t - ETA[seg]) / (ETA[seg+1] - ETA[seg]) : 0.0;
        if (tl < 0.0) tl = 0.0;
        if (tl > 1.0) tl = 1.0;
        const Airfoil& a0 = w.sections[seg];
        const Airfoil& a1 = (K > 1) ? w.sections[seg + 1] : a0;
        std::size_t nw = std::max(a0.wu.size(), a1.wu.size());
        s.af.N1 = a0.N1; s.af.N2 = a0.N2;
        s.af.wu.resize(nw); s.af.wl.resize(nw);
        for (std::size_t j = 0; j < nw; ++j) {
            s.af.wu[j] = (1.0-tl) * (j < a0.wu.size() ? a0.wu[j] : 0.0)
                        + tl      * (j < a1.wu.size() ? a1.wu[j] : 0.0);
            s.af.wl[j] = (1.0-tl) * (j < a0.wl.size() ? a0.wl[j] : 0.0)
                        + tl      * (j < a1.wl.size() ? a1.wl[j] : 0.0);
        }
        s.af.te_thick = (1.0-tl) * a0.te_thick + tl * a1.te_thick;
        double y_lo = (i == 0)     ? yp[0]     : 0.5 * (yp[i-1] + yp[i]);
        double y_hi = (i == n - 1) ? yp[n - 1] : 0.5 * (yp[i]   + yp[i+1]);
        s.width = y_hi - y_lo;
        // Arc-length twin of width, using the pre-projection arc coord y[] --
        // real skin/spar material follows the curve, not its cos(dihedral) shadow.
        double s_lo = (i == 0)     ? y[0]     : 0.5 * (y[i-1] + y[i]);
        double s_hi = (i == n - 1) ? y[n - 1] : 0.5 * (y[i]   + y[i+1]);
        s.ds = s_hi - s_lo;
        w.stations[i] = s;
    }
}

void planform(const WingGeometry& w, double& S_ref, double& mac,
              double& x_mac_le, double& b_full, double& AR) {
    // Half-wing integrals from stations, then doubled.
    double S_half = 0.0, mac_num = 0.0, xle_num = 0.0;
    for (const auto& s : w.stations) {
        S_half  += s.chord * s.width;
        mac_num += s.chord * s.chord * s.width;
        xle_num += s.x_le * s.chord * s.width;
    }
    S_ref = 2.0 * S_half;
    mac = (S_half > 0) ? mac_num / S_half : w.root_chord;
    x_mac_le = (S_half > 0) ? xle_num / S_half : 0.0;
    // ponytail: b_full must be the projected span (station.y, arc-integrated
    // via cos(dihedral) in loft()), not 2*semi_span -- that's the arc-length
    // parametrization coordinate, which overstates span whenever dihedral is
    // nonzero and silently inflates AR/CDi credit for folded designs.
    b_full = w.stations.empty() ? 2.0 * w.semi_span : 2.0 * w.stations.back().y;
    AR = (S_ref > 0) ? (b_full * b_full) / S_ref : 0.0;
}

}  // namespace geom
}  // namespace aero
