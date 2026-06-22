// build_surrogate.cpp — OFFLINE viscous-surrogate generator (separate exe).
//
// Samples CST airfoil shapes (seeds + a DoE around them) x a Reynolds grid,
// calibrates a compact polar (cd0, k, cl_max, cl_min, cm0) for each, and writes
// data/surrogates/polar_coeffs.csv for the runtime to interpolate.
//
//   surrogate_mode = synthetic   -> analytic calibration (no XFOIL; runs now)
//   surrogate_mode = xfoil       -> drive xfoil.exe per (shape, Re)
//
// The runtime binary never depends on XFOIL — only this offline tool does.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>

#include "aeroanalyzer/config.h"
#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/airfoil_io.h"
#include <Eigen/Dense>

using namespace aero;

namespace {

Airfoil shape_to_airfoil(const std::vector<double>& shape, double te) {
    Airfoil f;
    f.wu.assign(shape.begin(), shape.begin() + 4);
    f.wl.assign(shape.begin() + 4, shape.begin() + 8);
    f.te_thick = te;
    return f;
}

double max_thickness(const Airfoil& f) {
    double tmax = 0.0;
    for (int i = 0; i <= 100; ++i) {
        double x = 0.5 * (1.0 - std::cos(PI * i / 100));
        tmax = std::max(tmax, geom::cst_upper(f, x) - geom::cst_lower(f, x));
    }
    return tmax;
}

struct Coeffs { double cd0, k, cl_max, cl_min, cm0; };

// Shape/Re-dependent analytic calibration (stand-in for XFOIL).
Coeffs synthetic_coeffs(const Airfoil& f, double Re, double ncrit) {
    geom::ThinAirfoil ta = geom::thin_airfoil(f);
    double tmax = max_thickness(f);
    double camber = -ta.alpha_L0;   // positive for cambered sections
    double rough = std::max(0.0, (9.0 - ncrit)) * 0.0006;
    Coeffs c;
    c.cd0 = 0.0065 + 0.012 * std::pow(1.0e5 / std::max(Re, 1e4), 0.5) + rough
            + 0.02 * std::max(0.0, tmax - 0.10);
    c.k = 0.012 + 0.02 * tmax;
    c.cl_max = 1.15 - 0.10 * std::max(0.0, (9.0 - ncrit)) / 5.0
               + 2.0 * camber - 0.6 * std::max(0.0, tmax - 0.15);
    c.cl_max = std::max(0.6, std::min(1.8, c.cl_max));
    c.cl_min = -0.65;
    c.cm0 = ta.cm_ac;
    return c;
}

void write_dat(const std::string& path, const Airfoil& f) {
    std::ofstream o(path);
    o << "sample\n";
    for (int i = 80; i >= 0; --i) {
        double x = 0.5 * (1.0 - std::cos(PI * i / 80));
        o << x << " " << geom::cst_upper(f, x) << "\n";
    }
    for (int i = 1; i <= 80; ++i) {
        double x = 0.5 * (1.0 - std::cos(PI * i / 80));
        o << x << " " << geom::cst_lower(f, x) << "\n";
    }
}

// Parse an XFOIL accumulated polar dump: rows of alpha CL CD CDp CM ...
bool parse_polar(const std::string& path, std::vector<double>& cl,
                 std::vector<double>& cd, std::vector<double>& cm) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        std::vector<double> v;
        double x;
        while (ss >> x) v.push_back(x);
        if (v.size() >= 5 && v[0] >= -30 && v[0] <= 40)  // alpha, CL, CD, CDp, CM
            { cl.push_back(v[1]); cd.push_back(v[2]); cm.push_back(v[4]); }
    }
    return cl.size() >= 4;
}

bool fit_polar(const std::vector<double>& cl, const std::vector<double>& cd,
               const std::vector<double>& cm, Coeffs& out) {
    // cd = cd0 + k cl^2 over the unstalled band (|cl| < 1.0)
    Eigen::MatrixXd A(static_cast<Eigen::Index>(cl.size()), 2);
    Eigen::VectorXd b(static_cast<Eigen::Index>(cl.size()));
    Eigen::Index n = 0;
    for (std::size_t i = 0; i < cl.size(); ++i) {
        if (std::fabs(cl[i]) > 1.0) continue;
        A(n, 0) = 1.0; A(n, 1) = cl[i] * cl[i]; b(n) = cd[i]; ++n;
    }
    if (n < 3) return false;
    Eigen::VectorXd x =
        A.topRows(n).colPivHouseholderQr().solve(b.head(n));
    out.cd0 = x(0);
    out.k = std::max(0.0, x(1));
    out.cl_max = *std::max_element(cl.begin(), cl.end());
    out.cl_min = *std::min_element(cl.begin(), cl.end());
    // cm0 at the smallest |cl|
    std::size_t imin = 0;
    for (std::size_t i = 1; i < cl.size(); ++i)
        if (std::fabs(cl[i]) < std::fabs(cl[imin])) imin = i;
    out.cm0 = cm[imin];
    return true;
}

bool xfoil_coeffs(const Airfoil& f, double Re, double ncrit,
                  const std::string& xfoil_exe, const std::string& tmp,
                  Coeffs& out) {
    std::string dat = tmp + "/s.dat";
    std::string pol = tmp + "/s.pol";
    std::string cmd = tmp + "/s.cmd";
    std::error_code ec;
    std::filesystem::remove(pol, ec);
    write_dat(dat, f);
    {
        std::ofstream o(cmd);
        o << "LOAD " << dat << "\n";
        o << "PANE\n";
        o << "OPER\n";
        o << "VPAR\nN " << ncrit << "\n\n";
        o << "VISC " << static_cast<long>(Re) << "\n";
        o << "ITER 250\n";
        o << "PACC\n" << pol << "\n\n";
        o << "ASEQ -8 14 0.5\n";
        o << "PACC\n\n";
        o << "QUIT\n";
    }
    std::string run = "\"" + xfoil_exe + "\" < \"" + cmd + "\" > \"" + tmp +
                      "/x.log\" 2>&1";
    std::system(run.c_str());
    std::vector<double> cl, cd, cm;
    if (!parse_polar(pol, cl, cd, cm)) return false;
    return fit_polar(cl, cd, cm, out);
}

}  // namespace

int main(int argc, char** argv) {
    std::string cfg_path = (argc > 1) ? argv[1] : "config/baseline.cfg";
    Config cfg;
    if (!cfg.load(cfg_path))
        std::cerr << "[warn] could not open " << cfg_path << " (defaults)\n";

    geom::GenomeSpec spec = geom::default_genome();
    seeds::SeedSet seedset = seeds::load_seeds(cfg);
    seeds::widen_cst_bounds(spec, seedset, cfg.getd("cst_bound_margin", 0.04));

    // ---- assemble sample shapes: seeds + random DoE within CST bounds ----
    std::vector<std::vector<double>> shapes;
    for (const auto& f : seedset.airfoils) {
        std::vector<double> s(f.wu.begin(), f.wu.end());
        s.insert(s.end(), f.wl.begin(), f.wl.end());
        if (s.size() == 8) shapes.push_back(s);
    }
    int n_doe = cfg.geti("surrogate_samples", 300);
    unsigned seed = static_cast<unsigned>(cfg.geti("ga_seed", 1));
    std::mt19937 g(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const int sg[8] = {geom::G_WU0, geom::G_WU1, geom::G_WU2, geom::G_WU3,
                       geom::G_WL0, geom::G_WL1, geom::G_WL2, geom::G_WL3};
    // light Latin-Hypercube: one stratified permutation per dimension
    std::vector<std::vector<int>> perm(8, std::vector<int>(n_doe));
    for (int d = 0; d < 8; ++d) {
        for (int i = 0; i < n_doe; ++i) perm[d][i] = i;
        std::shuffle(perm[d].begin(), perm[d].end(), g);
    }
    for (int i = 0; i < n_doe; ++i) {
        std::vector<double> s(8);
        for (int d = 0; d < 8; ++d) {
            double lo = spec.lo[sg[d]], hi = spec.hi[sg[d]];
            double q = (perm[d][i] + u(g)) / n_doe;   // stratified
            s[d] = lo + q * (hi - lo);
        }
        shapes.push_back(s);
    }

    // ---- Reynolds grid ----
    std::vector<double> res;
    {
        std::stringstream ss(cfg.gets("surrogate_re", "100000,150000,200000,300000"));
        std::string t;
        while (std::getline(ss, t, ',')) { try { res.push_back(std::stod(t)); } catch (...) {} }
        if (res.empty()) res = {1e5, 1.5e5, 2e5, 3e5};
    }

    std::string mode = cfg.gets("surrogate_mode", "synthetic");
    double ncrit = cfg.getd("ncrit", 4.0);
    std::string xfoil_exe = cfg.gets("xfoil_exe", "xfoil.exe");
    std::filesystem::create_directories("data/surrogates");
    std::filesystem::create_directories("build/xfoil_tmp");

    bool use_xfoil = (mode == "xfoil");
    if (use_xfoil) {
        // probe: does xfoil run at all?
        if (std::system(("\"" + xfoil_exe + "\" < nul > nul 2>&1").c_str()) != 0)
            std::cerr << "[warn] '" << xfoil_exe
                      << "' did not run; rows that fail fall back to synthetic.\n";
    }

    std::ofstream csv("data/surrogates/polar_coeffs.csv");
    csv << "# wu0,wu1,wu2,wu3,wl0,wl1,wl2,wl3,Re,cd0,k,cl_max,cl_min,cm0\n";
    csv << "# mode=" << mode << " ncrit=" << ncrit << " shapes=" << shapes.size()
        << " Re=" << res.size() << "\n";

    int rows = 0, xf_ok = 0, xf_fail = 0;
    for (const auto& s : shapes) {
        Airfoil f = shape_to_airfoil(s, 0.002);
        for (double Re : res) {
            Coeffs c;
            bool ok = false;
            if (use_xfoil)
                ok = xfoil_coeffs(f, Re, ncrit, xfoil_exe, "build/xfoil_tmp", c);
            if (ok) ++xf_ok;
            else {
                if (use_xfoil) ++xf_fail;
                c = synthetic_coeffs(f, Re, ncrit);   // fallback / synthetic
            }
            csv << s[0] << "," << s[1] << "," << s[2] << "," << s[3] << ","
                << s[4] << "," << s[5] << "," << s[6] << "," << s[7] << ","
                << static_cast<long>(Re) << "," << c.cd0 << "," << c.k << ","
                << c.cl_max << "," << c.cl_min << "," << c.cm0 << "\n";
            ++rows;
        }
    }
    csv.close();
    std::cout << "wrote data/surrogates/polar_coeffs.csv: " << rows << " rows ("
              << shapes.size() << " shapes x " << res.size() << " Re)\n";
    if (use_xfoil)
        std::cout << "  xfoil: " << xf_ok << " ok, " << xf_fail
                  << " fell back to synthetic\n";
    return 0;
}
