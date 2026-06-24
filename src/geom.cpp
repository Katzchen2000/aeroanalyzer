#include "aeroanalyzer/geom.h"
#include <Eigen/Dense>
#include <cmath>

namespace aero {
namespace geom {

GenomeSpec default_genome() {
    GenomeSpec g;
    g.lo.resize(N_GENES);
    g.hi.resize(N_GENES);
    g.names.resize(N_GENES);
    auto set = [&](int i, const char* nm, double lo, double hi) {
        g.names[i] = nm; g.lo[i] = lo; g.hi[i] = hi;
    };
    set(G_ROOT,     "root_chord_m", 0.18, 0.35);
    set(G_TAPER,    "taper_ratio",  0.30, 0.90);
    set(G_SEMISPAN, "semi_span_m",  0.45, 0.80);
    set(G_SWEEP,    "le_sweep_deg",  8.0, 30.0);
    set(G_WASHOUT,  "washout_deg",  -6.0,  0.0);
    set(G_BATTERY,  "battery_x_m",   0.00, 0.22);
    set(G_WU0, "wu0", 0.10, 0.34);
    set(G_WU1, "wu1", 0.10, 0.34);
    set(G_WU2, "wu2", 0.08, 0.30);
    set(G_WU3, "wu3", 0.06, 0.26);
    set(G_WL0, "wl0", -0.22, 0.05);
    set(G_WL1, "wl1", -0.20, 0.08);
    set(G_WL2, "wl2", -0.15, 0.15);
    set(G_WL3, "wl3", -0.05, 0.20);   // aft lower weight drives reflex
    set(G_TE,       "te_frac",       0.002, 0.010);
    set(G_MODE,     "mode",          0.0,   1.0);
    set(G_CS_CHORD, "cs_chord_frac", 0.15,  0.35);
    set(G_AIL_SPAN, "ail_span_frac", 0.40,  0.80);
    // tip section mirrors the root bounds; loft() interpolates root->tip.
    set(G_TIP_WU0, "tip_wu0", 0.10, 0.34);
    set(G_TIP_WU1, "tip_wu1", 0.10, 0.34);
    set(G_TIP_WU2, "tip_wu2", 0.08, 0.30);
    set(G_TIP_WU3, "tip_wu3", 0.06, 0.26);
    set(G_TIP_WL0, "tip_wl0", -0.22, 0.05);
    set(G_TIP_WL1, "tip_wl1", -0.20, 0.08);
    set(G_TIP_WL2, "tip_wl2", -0.15, 0.15);
    set(G_TIP_WL3, "tip_wl3", -0.05, 0.20);
    return g;
}

double bernstein(int i, int n, double x) {
    // binomial(n,i) x^i (1-x)^(n-i)
    double c = 1.0;
    for (int k = 0; k < i; ++k) c = c * (n - k) / (k + 1);
    return c * std::pow(x, i) * std::pow(1.0 - x, n - i);
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

WingGeometry decode(const std::vector<double>& g, const GenomeSpec& spec) {
    auto clamp = [&](int i) {
        double v = g[i];
        if (v < spec.lo[i]) v = spec.lo[i];
        if (v > spec.hi[i]) v = spec.hi[i];
        return v;
    };
    WingGeometry w;
    w.root_chord = clamp(G_ROOT);
    w.tip_chord  = w.root_chord * clamp(G_TAPER);
    w.semi_span  = clamp(G_SEMISPAN);
    w.le_sweep   = clamp(G_SWEEP) * DEG2RAD;
    w.washout    = clamp(G_WASHOUT) * DEG2RAD;
    w.battery_x  = clamp(G_BATTERY);
    w.section.wu = {clamp(G_WU0), clamp(G_WU1), clamp(G_WU2), clamp(G_WU3)};
    w.section.wl = {clamp(G_WL0), clamp(G_WL1), clamp(G_WL2), clamp(G_WL3)};
    w.section.te_thick = clamp(G_TE);
    w.section_tip.wu = {clamp(G_TIP_WU0), clamp(G_TIP_WU1), clamp(G_TIP_WU2), clamp(G_TIP_WU3)};
    w.section_tip.wl = {clamp(G_TIP_WL0), clamp(G_TIP_WL1), clamp(G_TIP_WL2), clamp(G_TIP_WL3)};
    w.section_tip.te_thick = w.section.te_thick;  // evaluate() may override w/ te_thick_tip_frac
    w.mode           = (clamp(G_MODE) < 0.5) ? ControlMode::Elevon : ControlMode::Split;
    w.cs_chord_frac  = clamp(G_CS_CHORD);
    w.ail_span_frac  = clamp(G_AIL_SPAN);
    return w;
}

void loft(WingGeometry& w, int n) {
    if (n < 2) n = 2;
    w.stations.assign(n, Station{});
    // Tip section blends to root when unset (most scratch/test callers only
    // set w.section); decode() populates section_tip for the optimizer.
    const Airfoil& root = w.section;
    const Airfoil& tip  = w.section_tip.wu.empty() ? w.section : w.section_tip;
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double th = PI * i / (n - 1);
        y[i] = w.semi_span * 0.5 * (1.0 - std::cos(th));  // cluster root & tip
    }
    for (int i = 0; i < n; ++i) {
        double t = (w.semi_span > 0) ? y[i] / w.semi_span : 0.0;  // 0 root .. 1 tip
        Station s;
        s.y = y[i];
        s.chord = w.root_chord + (w.tip_chord - w.root_chord) * t;
        s.x_le  = y[i] * std::tan(w.le_sweep);            // sheared sweep
        s.twist = w.washout * t;
        s.z = 0.0;
        // lofted section: linear blend of root/tip CST weights + TE thickness.
        s.af.N1 = root.N1; s.af.N2 = root.N2;
        s.af.wu.resize(root.wu.size());
        s.af.wl.resize(root.wl.size());
        for (std::size_t j = 0; j < root.wu.size(); ++j)
            s.af.wu[j] = (1.0 - t) * root.wu[j] + t * tip.wu[j];
        for (std::size_t j = 0; j < root.wl.size(); ++j)
            s.af.wl[j] = (1.0 - t) * root.wl[j] + t * tip.wl[j];
        s.af.te_thick = (1.0 - t) * root.te_thick + t * tip.te_thick;
        // strip width via midpoints (trapezoidal control volumes)
        double y_lo = (i == 0) ? y[0] : 0.5 * (y[i - 1] + y[i]);
        double y_hi = (i == n - 1) ? y[n - 1] : 0.5 * (y[i] + y[i + 1]);
        s.width = y_hi - y_lo;
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
    b_full = 2.0 * w.semi_span;
    AR = (S_ref > 0) ? (b_full * b_full) / S_ref : 0.0;
}

}  // namespace geom
}  // namespace aero
