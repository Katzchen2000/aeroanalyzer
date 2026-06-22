#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/airfoil_io.h"
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
        for (int i = 0; i < 4 && i < (int)f.wu.size(); ++i) widen(G_WU0 + i, f.wu[i]);
        for (int i = 0; i < 4 && i < (int)f.wl.size(); ++i) widen(G_WL0 + i, f.wl[i]);
        widen(G_TE, f.te_thick);
    }
    // keep TE strictly positive
    spec.lo[G_TE] = std::max(spec.lo[G_TE], 1e-4);
}

std::vector<std::vector<double>> build_seed_genomes(
    const geom::GenomeSpec& spec, const SeedSet& s, unsigned rng_seed, int count) {
    std::vector<std::vector<double>> out;
    if (s.airfoils.empty() || count <= 0) return out;
    std::mt19937 g(rng_seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const int n = static_cast<int>(spec.size());
    const int n_elite = std::max(1, count / 5);
    using namespace geom;
    const int planform[] = {G_ROOT, G_TAPER, G_SEMISPAN, G_SWEEP, G_WASHOUT, G_BATTERY};

    auto clampg = [&](int gi, double v) {
        return std::min(spec.hi[gi], std::max(spec.lo[gi], v));
    };

    for (int k = 0; k < count; ++k) {
        std::vector<double> x(n);
        for (int gi = 0; gi < n; ++gi)
            x[gi] = spec.lo[gi] + u(g) * (spec.hi[gi] - spec.lo[gi]);  // explorer base

        // overlay the seed airfoil's CST shape
        const Airfoil& af = s.airfoils[k % s.airfoils.size()];
        for (int i = 0; i < 4 && i < (int)af.wu.size(); ++i) x[G_WU0 + i] = clampg(G_WU0 + i, af.wu[i]);
        for (int i = 0; i < 4 && i < (int)af.wl.size(); ++i) x[G_WL0 + i] = clampg(G_WL0 + i, af.wl[i]);
        x[G_TE] = clampg(G_TE, af.te_thick);

        // elites: settle planform at mid-box; hybrids keep the random planform
        if (k < n_elite)
            for (int gi : planform) x[gi] = 0.5 * (spec.lo[gi] + spec.hi[gi]);

        out.push_back(std::move(x));
    }
    return out;
}

}  // namespace seeds
}  // namespace aero
