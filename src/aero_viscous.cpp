#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/engine_core.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <limits>
#include <iostream>

namespace aero {
namespace viscous {

Polar Surrogate::analytic(double cl, double Re) const {
    Polar p;
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
    samples_.clear();
    tables_loaded_ = false;

    std::ifstream f(dir + "/polar_coeffs.csv");
    if (!f) return false;

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
    if (samples_.empty()) return false;

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

Polar Surrogate::query(const std::vector<double>& shape, double cl,
                       double Re) const {
    if (!tables_loaded_ || shape.size() != hull_lo_.size())
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
    if (out_of_hull) p.clamped = true;   // flag shape/Re extrapolation
    return p;
}

}  // namespace viscous
}  // namespace aero
