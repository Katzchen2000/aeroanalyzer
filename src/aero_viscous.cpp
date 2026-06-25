#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/geom.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <limits>
#include <iostream>

namespace aero {
namespace viscous {

Polar Surrogate::analytic(double cl, double Re) const {
    Polar p;
    p.confidence = 1.0;
    Re = std::max(Re, 1.0e4);
    double rough = std::max(0.0, (9.0 - ncrit_)) * 0.0006;
    double cd0 = 0.0065 + 0.012 * std::pow(1.0e5 / Re, 0.5) + rough;
    double k = 0.013;
    p.cl_max = 1.20 - 0.10 * std::max(0.0, (9.0 - ncrit_)) / 5.0;
    p.cl_min = -0.65;

    double cl_eff = cl;
    if (cl > p.cl_max)      { cl_eff = p.cl_max; p.clamped = true; }
    else if (cl < p.cl_min) { cl_eff = p.cl_min; p.clamped = true; }
    p.cd = cd0 + k * cl_eff * cl_eff;
    if (p.clamped) {
        double over = std::fabs(cl - cl_eff);
        p.cd += 0.9 * over * over + 0.05 * over;
    }
    p.cm = 0.0;
    return p;
}

bool Surrogate::load(const std::string& dir, const Config& cfg) {
    ncrit_ = cfg.getd("ncrit", 4.0);
    xtr_upper_ = cfg.getd("nf_xtr_upper", 1.0);
    xtr_lower_ = cfg.getd("nf_xtr_lower", 1.0);
    samples_.clear();
    tables_loaded_ = false;

    const std::string be = cfg.gets("viscous_backend", "neuralfoil");
    if (be == "analytic") {
        backend_ = Backend::Analytic;
        std::cout << "[surrogate] analytic backend (no data tables)\n";
        return false;
    }
    if (be == "neuralfoil") {
        const std::string nfdir = cfg.gets("neuralfoil_dir", "data/Neurafoilbin");
        if (nf_.load(nfdir)) {
            backend_ = Backend::NeuralFoil;
            std::cout << "[surrogate] NeuralFoil backend loaded from " << nfdir
                      << "\n";
            return true;
        }
        std::cerr << "[surrogate] NeuralFoil weights unavailable in " << nfdir
                  << "; falling back to polar table\n";
    }

    // ---- Table backend: legacy offline XFOIL polar CSV + IDW interpolation ----
    backend_ = Backend::Table;
    std::ifstream f(dir + "/polar_coeffs.csv");
    if (!f) { backend_ = Backend::Analytic; return false; }

    std::string line;
    std::size_t shape_dim = 0;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        // header row: starts with a non-digit token
        if (!(std::isdigit((unsigned char)line[0]) || line[0] == '-' ||
              line[0] == '+' || line[0] == '.'))
            continue;
        std::vector<double> vals;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            try { vals.push_back(std::stod(tok)); } catch (...) { }
        }
        // layout: shape... , Re, cd0, k, cl_max, cl_min, cm0  (>= 7 columns)
        if (vals.size() < 7) continue;
        std::size_t sd = vals.size() - 6;
        if (shape_dim == 0) shape_dim = sd;
        if (sd != shape_dim) continue;
        Sample s;
        s.shape.assign(vals.begin(), vals.begin() + sd);
        s.Re     = vals[sd + 0];
        s.cd0    = vals[sd + 1];
        s.k      = vals[sd + 2];
        s.cl_max = vals[sd + 3];
        s.cl_min = vals[sd + 4];
        s.cm0    = vals[sd + 5];
        samples_.push_back(std::move(s));
    }
    if (samples_.empty()) { backend_ = Backend::Analytic; return false; }

    // per-dimension hull for the extrapolation guard
    hull_lo_.assign(shape_dim, std::numeric_limits<double>::infinity());
    hull_hi_.assign(shape_dim, -std::numeric_limits<double>::infinity());
    re_lo_ = std::numeric_limits<double>::infinity();
    re_hi_ = -std::numeric_limits<double>::infinity();
    for (const auto& s : samples_) {
        for (std::size_t d = 0; d < shape_dim; ++d) {
            hull_lo_[d] = std::min(hull_lo_[d], s.shape[d]);
            hull_hi_[d] = std::max(hull_hi_[d], s.shape[d]);
        }
        re_lo_ = std::min(re_lo_, s.Re);
        re_hi_ = std::max(re_hi_, s.Re);
    }
    tables_loaded_ = true;
    std::cout << "[surrogate] loaded " << samples_.size() << " samples ("
              << shape_dim << "-D shape) from " << dir << "/polar_coeffs.csv\n";
    return true;
}

// FNV-1a-style hash of the section shape (8 CST weights) + te, quantized so that
// numerically-identical candidates share a cache slot.
static std::size_t shape_sig(const std::vector<double>& shape, double te) {
    std::size_t h = 1469598103934665603ull;
    auto mix = [&](double v) {
        std::int64_t q = static_cast<std::int64_t>(std::llround(v * 1.0e6));
        h ^= static_cast<std::size_t>(q);
        h *= 1099511628211ull;
    };
    for (std::size_t i = 0; i < 8 && i < shape.size(); ++i) mix(shape[i]);
    mix(te);
    return h;
}

void Surrogate::nf_coeffs(const std::vector<double>& shape, double te, double Re,
                          double& cd0, double& k, double& cl_max, double& cl_min,
                          double& cm0, double& confidence) const {
    // ---- 4+4 project CST -> 8+8 Kulfan refit. Shared by every station of a
    //      candidate (same shape), so memoize; the GA runs under OpenMP, hence
    //      thread_local (no shared-state race on the const query path). ----
    struct Kulfan {
        std::array<double, 8> up{}, lo{};
        std::size_t sig = 0;
        bool valid = false;
    };
    thread_local Kulfan kf;
    const std::size_t sig = shape_sig(shape, te);
    if (!kf.valid || kf.sig != sig) {
        Airfoil src;
        src.wu.assign(shape.begin(), shape.begin() + 4);
        src.wl.assign(shape.begin() + 4, shape.begin() + 8);
        src.te_thick = te; src.N1 = 0.5; src.N2 = 1.0;
        const int M = 60;
        std::vector<std::pair<double, double>> up_pts(M), lo_pts(M);
        for (int i = 0; i < M; ++i) {
            double x = 0.5 * (1.0 - std::cos(PI * i / (M - 1)));
            up_pts[i] = {x, geom::cst_upper(src, x)};
            lo_pts[i] = {x, geom::cst_lower(src, x)};
        }
        Airfoil g = geom::fit_cst(up_pts, lo_pts, 7, te);
        for (int i = 0; i < 8; ++i) {
            kf.up[i] = (i < static_cast<int>(g.wu.size())) ? g.wu[i] : 0.0;
            kf.lo[i] = (i < static_cast<int>(g.wl.size())) ? g.wl[i] : 0.0;
        }
        kf.sig = sig; kf.valid = true;
    }

    // ---- polar cache keyed on (shape sig, Re bucket). Must exceed the spanwise
    //      station count: trim() re-runs solve() ~15x per candidate (Newton +
    //      finite-difference Jacobian), and every station hits a distinct Re
    //      bucket, so a too-small ring thrashes and re-sweeps the whole net on
    //      every solve. 64 comfortably holds one candidate's working set. ----
    struct Coeff {
        std::size_t key = 0;
        double cd0 = 0, k = 0, clmax = 0, clmin = 0, cm0 = 0;
        double confidence = 1.0;
        bool valid = false;
    };
    constexpr int NC = 64;
    thread_local std::array<Coeff, NC> cache{};
    thread_local int rr = 0;
    const std::size_t rebkt = static_cast<std::size_t>(std::llround(Re / 2000.0));
    const std::size_t key = sig ^ (rebkt * 2654435761ull);
    for (const auto& c : cache)
        if (c.valid && c.key == key) {
            cd0 = c.cd0; k = c.k; cl_max = c.clmax; cl_min = c.clmin;
            cm0 = c.cm0; confidence = c.confidence;
            return;
        }

    // ---- alpha sweep through NeuralFoil + parabola fit (cd = cd0 + k cl^2) ----
    const double a0 = -8.0, a1 = 16.0, da = 1.0;
    const double Rc = std::max(Re, 1.0e4);
    double S0 = 0, Sx = 0, Sxx = 0, Sy = 0, Sxy = 0;
    double clmax = -1e9, clmin = 1e9, cm_at_min = 0.0, min_abscl = 1e9;
    double conf_band = 1.0;
    for (double a = a0; a <= a1 + 1e-9; a += da) {
        nf::Aero r = nf_.eval(kf.up, kf.lo, 0.0, te, a, Rc, ncrit_,
                              xtr_upper_, xtr_lower_);
        if (!std::isfinite(r.cl) || !std::isfinite(r.cd)) continue;
        if (r.confidence > 0.5) {            // stall envelope where trusted
            clmax = std::max(clmax, r.cl);
            clmin = std::min(clmin, r.cl);
        }
        if (std::fabs(r.cl) < 1.0) {         // unstalled parabola band
            double x = r.cl * r.cl;
            S0 += 1; Sx += x; Sxx += x * x; Sy += r.cd; Sxy += x * r.cd;
            if (std::fabs(r.cl) < min_abscl) {
                min_abscl = std::fabs(r.cl); cm_at_min = r.cm;
            }
            conf_band = std::min(conf_band, r.confidence);
        }
    }
    const double det = S0 * Sxx - Sx * Sx;
    if (S0 >= 3 && std::fabs(det) > 1e-12) {
        cd0 = (Sy * Sxx - Sx * Sxy) / det;
        k   = (S0 * Sxy - Sx * Sy) / det;
    } else {
        cd0 = 0.012; k = 0.02;               // degenerate sweep -> safe defaults
    }
    if (k < 0.0) k = 0.0;
    if (cd0 < 1e-4) cd0 = 1e-4;
    cl_max = (clmax > -1e8) ? clmax : 1.2;
    cl_min = (clmin <  1e8) ? clmin : -0.6;
    cm0 = cm_at_min;
    confidence = conf_band;

    cache[rr] = Coeff{key, cd0, k, cl_max, cl_min, cm0, confidence, true};
    rr = (rr + 1) % NC;
}

Polar Surrogate::query(const std::vector<double>& shape, double cl,
                       double Re, double te) const {
    if (backend_ == Backend::NeuralFoil && nf_.loaded() && shape.size() >= 8) {
        double cd0, k, cmx, cmn, cm0, confidence = 1.0;
        nf_coeffs(shape, te, Re, cd0, k, cmx, cmn, cm0, confidence);
        Polar p;
        p.cl_max = cmx; p.cl_min = cmn; p.cm = cm0; p.confidence = confidence;
        double cl_eff = cl;
        if (cl > p.cl_max)      { cl_eff = p.cl_max; p.clamped = true; }
        else if (cl < p.cl_min) { cl_eff = p.cl_min; p.clamped = true; }
        p.cd = cd0 + k * cl_eff * cl_eff;
        if (p.clamped) {
            double over = std::fabs(cl - cl_eff);
            p.cd += 0.9 * over * over + 0.05 * over;
        }
        if (confidence < 0.5) p.clamped = true;  // out-of-distribution flag
        return p;
    }

    if (backend_ != Backend::Table || !tables_loaded_ ||
        shape.size() != hull_lo_.size())
        return analytic(cl, Re);

    const std::size_t D = shape.size();
    // inverse-distance weighting over (shape, Re), normalized per dimension
    double cd0 = 0, k = 0, cmx = 0, cmn = 0, cm0 = 0, wsum = 0;
    bool out_of_hull = false;
    for (std::size_t d = 0; d < D; ++d) {
        if (shape[d] < hull_lo_[d] - 1e-9 || shape[d] > hull_hi_[d] + 1e-9)
            out_of_hull = true;
    }
    if (Re < re_lo_ - 1.0 || Re > re_hi_ + 1.0) out_of_hull = true;

    for (const auto& s : samples_) {
        double d2 = 0.0;
        for (std::size_t d = 0; d < D; ++d) {
            double range = std::max(hull_hi_[d] - hull_lo_[d], 1e-6);
            double dd = (shape[d] - s.shape[d]) / range;
            d2 += dd * dd;
        }
        double re_range = std::max(re_hi_ - re_lo_, 1.0);
        double dre = (Re - s.Re) / re_range;
        d2 += dre * dre;
        double w = 1.0 / (d2 + 1e-9);
        cd0 += w * s.cd0; k += w * s.k; cmx += w * s.cl_max;
        cmn += w * s.cl_min; cm0 += w * s.cm0; wsum += w;
    }
    Polar p;
    if (wsum <= 0) return analytic(cl, Re);
    cd0 /= wsum; k /= wsum; cmx /= wsum; cmn /= wsum; cm0 /= wsum;
    p.cl_max = cmx; p.cl_min = cmn; p.cm = cm0;

    double cl_eff = cl;
    if (cl > p.cl_max)      { cl_eff = p.cl_max; p.clamped = true; }
    else if (cl < p.cl_min) { cl_eff = p.cl_min; p.clamped = true; }
    p.cd = cd0 + k * cl_eff * cl_eff;
    if (p.clamped) {
        double over = std::fabs(cl - cl_eff);
        p.cd += 0.9 * over * over + 0.05 * over;
    }
    if (out_of_hull) {
        p.clamped = true;   // flag shape/Re extrapolation
        p.confidence = 0.0;
    }
    return p;
}

}  // namespace viscous
}  // namespace aero
