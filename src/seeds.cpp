#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/airfoil_io.h"
#include <Eigen/Dense>
#include <filesystem>
#include <random>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace aero {
namespace seeds {

namespace {
std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        // trim
        std::size_t a = tok.find_first_not_of(" \t");
        std::size_t b = tok.find_last_not_of(" \t");
        if (a != std::string::npos) out.push_back(tok.substr(a, b - a + 1));
    }
    return out;
}
}  // namespace

SeedSet load_seeds(const Config& cfg) {
    SeedSet set;
    const int order = 3;   // must match genome: 4 upper + 4 lower weights

    // (1) uploaded coordinate files (.dat or .txt; Selig/Lednicer auto-detected)
    std::string dir = cfg.gets("seed_airfoils_dir", "data/airfoils");
    std::error_code ec;
    if (std::filesystem::exists(dir, ec)) {
        for (const auto& e : std::filesystem::directory_iterator(dir, ec)) {
            if (ec) break;
            std::string ext = e.path().extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(),
                           [](unsigned char ch) { return std::tolower(ch); });
            if (ext != ".dat" && ext != ".txt") continue;
            bool ok = false;
            airfoil_io::Coords c = airfoil_io::load_dat(e.path().string(), ok);
            if (!ok) {
                std::cerr << "[seeds] skipped unreadable " << e.path() << "\n";
                continue;
            }
            double te = airfoil_io::estimate_te(c);
            set.airfoils.push_back(airfoil_io::to_airfoil(c, order, te));
            set.names.push_back(c.name.empty() ? e.path().stem().string() : c.name);
        }
    }

    // (2) generated NACA 4-digit
    for (const auto& code : split_csv(cfg.gets("seed_naca", "0012,2412,4412"))) {
        if (code.size() != 4) {
            std::cerr << "[seeds] only NACA 4-digit supported, skipping " << code
                      << "\n";
            continue;
        }
        airfoil_io::Coords c = airfoil_io::naca4(code, 100);
        double te = airfoil_io::estimate_te(c);
        set.airfoils.push_back(airfoil_io::to_airfoil(c, order, te));
        set.names.push_back("NACA " + code);
    }
    return set;
}

void widen_cst_bounds(geom::GenomeSpec& spec, const SeedSet& s, double margin) {
    if (s.airfoils.empty()) return;
    auto widen = [&](int gi, double v) {
        spec.lo[gi] = std::min(spec.lo[gi], v - margin);
        spec.hi[gi] = std::max(spec.hi[gi], v + margin);
    };
    using namespace geom;
    for (const auto& f : s.airfoils) {
        for (int k = 0; k < N_SECTIONS; ++k) {
            for (int i = 0; i < 4 && i < (int)f.wu.size(); ++i)
                widen(G_SEC(k, 0, i), f.wu[i]);
            for (int i = 0; i < 4 && i < (int)f.wl.size(); ++i)
                widen(G_SEC(k, 1, i), f.wl[i]);
        }
    }
    // G_TE gene removed; no te_thick widening needed
}

std::vector<std::vector<double>> build_seed_genomes(
    const geom::GenomeSpec& spec, const SeedSet& s, unsigned rng_seed, int count,
    double cst_jitter) {
    std::vector<std::vector<double>> out;
    if (s.airfoils.empty() || count <= 0) return out;
    std::mt19937 g(rng_seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const int n = static_cast<int>(spec.size());
    const int n_elite = std::max(1, count / 5);
    using namespace geom;
    const int planform[] = {G_ROOT, G_TAPER, G_SEMISPAN, G_SWEEP, G_WASHOUT,
                            G_BATTERY, G_CS_CHORD, G_AIL_SPAN,
                            G_CHORD_EXP, G_SWEEP_EXP, G_GULL_A, G_GULL_B, G_GULL_C,
                            G_WINGLET_CANT, G_WINGLET_ETA};

    auto clampg = [&](int gi, double v) {
        return std::min(spec.hi[gi], std::max(spec.lo[gi], v));
    };

    // ---- PCA basis from seed CST vectors (4 wu + 4 wl = 8 dims) ----
    // Jitter along principal modes keeps hybrids on the realistic-airfoil manifold
    // instead of wandering into figure-8 shapes that gate 12 then culls.
    const int D = 8;
    const int ns = (int)s.airfoils.size();
    Eigen::MatrixXd M(ns, D);
    for (int i = 0; i < ns; ++i) {
        const Airfoil& a = s.airfoils[i];
        for (int j = 0; j < 4; ++j) M(i,   j) = j < (int)a.wu.size() ? a.wu[j] : 0.0;
        for (int j = 0; j < 4; ++j) M(i, 4+j) = j < (int)a.wl.size() ? a.wl[j] : 0.0;
    }
    Eigen::MatrixXd Mc = M.rowwise() - M.colwise().mean();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mc, Eigen::ComputeThinV);
    const Eigen::MatrixXd& V_modes = svd.matrixV();
    const Eigen::VectorXd& S_vals  = svd.singularValues();
    const int n_modes = (int)S_vals.size();

    std::normal_distribution<double> nrm(0.0, 1.0);

    // Curated feasibility anchor for elite seeds: forward battery + max washout
    // ensures NP aft of CG → SM lands near 6% without burning ~40 gen of spin-up.
    static const double g_anchor[N_PLANFORM] = {
        0.26, 0.55, 0.62, 20.0, -5.0,  // root, taper, semi_span, sweep, washout
        0.03, 0.25, 0.60,  1.2,  1.2,  // battery, cs_chord, ail_span, chord_exp, sweep_exp
        0.02, 0.00, 0.00,              // gull_a, gull_b, gull_c
        45.0, 0.85,                     // winglet_cant_deg, winglet_eta
    };

    for (int k = 0; k < count; ++k) {
        std::vector<double> x(n);
        for (int gi = 0; gi < n; ++gi)
            x[gi] = spec.lo[gi] + u(g) * (spec.hi[gi] - spec.lo[gi]);

        const Airfoil& af = s.airfoils[k % s.airfoils.size()];
        bool elite = (k < n_elite);

        // PCA jitter: independent per section; each section delta lies on the
        // realistic-airfoil manifold spanned by the seed library.
        // Falls back to per-weight uniform jitter when ns=1 (all S_vals=0).
        for (int sec = 0; sec < N_SECTIONS; ++sec) {
            Eigen::VectorXd delta = Eigen::VectorXd::Zero(D);
            if (!elite && cst_jitter > 0.0) {
                if (n_modes > 0 && S_vals(0) > 1e-8) {
                    for (int m = 0; m < n_modes; ++m)
                        delta += (nrm(g) * S_vals(m) * cst_jitter) * V_modes.col(m);
                } else {
                    for (int d = 0; d < D; ++d) {
                        int gi = (d < 4) ? G_SEC(0, 0, d) : G_SEC(0, 1, d - 4);
                        delta(d) = cst_jitter * (spec.hi[gi] - spec.lo[gi]) * (u(g) - 0.5);
                    }
                }
            }
            for (int i = 0; i < 4 && i < (int)af.wu.size(); ++i)
                x[G_SEC(sec, 0, i)] = clampg(G_SEC(sec, 0, i), af.wu[i] + delta(i));
            for (int i = 0; i < 4 && i < (int)af.wl.size(); ++i)
                x[G_SEC(sec, 1, i)] = clampg(G_SEC(sec, 1, i), af.wl[i] + delta(4+i));
        }

        if (elite)
            for (int gi : planform) x[gi] = clampg(gi, g_anchor[gi]);

        out.push_back(std::move(x));
    }
    return out;
}

}  // namespace seeds
}  // namespace aero
