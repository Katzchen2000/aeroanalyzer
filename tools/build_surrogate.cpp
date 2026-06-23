// build_surrogate.cpp — OFFLINE viscous-surrogate generator (separate exe).
//
// Samples CST airfoil shapes (seeds + a DoE around them) x a Reynolds grid,
// calibrates a compact polar (cd0, k, cl_max, cl_min, cm0) for each, and writes
// data/surrogates/polar_coeffs.csv for the runtime to interpolate.
//
//   surrogate_mode = synthetic   -> analytic calibration (no XFOIL; runs now)
//   surrogate_mode = native      -> in-process aero::xfoil solver (no subprocess,
//                                   no disk I/O; the M4 default)
//   surrogate_mode = xfoil       -> drive the legacy xfoil.exe per (shape, Re)
//   surrogate_mode = xfoil_lib   -> in-process REAL XFOIL 6.99 Fortran core
//                                   (libxfoil.a), zero disk I/O, xfoil.exe-exact.
//                                   XFOIL's COMMON state is not thread-safe, so
//                                   this mode fans out across PROCESSES via a
//                                   self-spawning worker pool (--worker). Falls
//                                   back per-row to native, then synthetic.
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
#include <thread>

#include "aeroanalyzer/config.h"
#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/airfoil_io.h"
#include "aeroanalyzer/aero_xfoil.h"
#include <Eigen/Dense>

#ifdef HAVE_XFOIL_LIB
#include "aeroanalyzer/xfoil_lib.h"
#include <cstdio>
#if defined(_WIN32)
#include <windows.h>
#include <process.h>
#else
#include <sys/wait.h>
#include <unistd.h>
#endif
#endif

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

// Native in-process polar: sweep the aero::xfoil solver and fit the compact
// polar from the converged points. No subprocess, no temp files.
bool native_coeffs(const Airfoil& f, double Re, double ncrit, Coeffs& out) {
    xfoil::Options opt;
    opt.Ncrit = ncrit;
    std::vector<xfoil::Result> sweep = xfoil::sweep(f, -6.0, 12.0, 1.0, Re, opt);
    std::vector<double> cl, cd, cm;
    for (const auto& r : sweep)
        if (r.converged && std::isfinite(r.cd) && std::isfinite(r.cl)) {
            cl.push_back(r.cl); cd.push_back(r.cd); cm.push_back(r.cm);
        }
    if (cl.size() < 4) return false;
    return fit_polar(cl, cd, cm, out);
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

#ifdef HAVE_XFOIL_LIB
// ---- accurate engine: real XFOIL 6.99 in-process (libxfoil.a) -------------
// One global Session per process (XFOIL COMMON state is process-global). A
// worker process owns exactly one of these.
aero::xfoil_lib::Session& xfl_session() {
    static aero::xfoil_lib::Session s;   // ctor runs xfl_init() once
    return s;
}

// Compute the compact polar for one (shape, Re) with the accurate engine.
// Returns true and fills `out` on success; false if too few points converged.
bool xfoil_lib_coeffs(const Airfoil& f, double Re, double ncrit, Coeffs& out) {
    std::vector<aero::xfoil_lib::Point> pol =
        xfl_session().polar(f, Re, ncrit, -6.0, 14.0, 1.0);
    std::vector<double> cl, cd, cm;
    for (const auto& p : pol)
        if (p.converged) { cl.push_back(p.cl); cd.push_back(p.cd); cm.push_back(p.cm); }
    if (cl.size() < 4) return false;
    return fit_polar(cl, cd, cm, out);
}

// Per-row fallback chain for the accurate mode: xfoil_lib -> native -> synthetic.
// Writes the source tag used into `src`.
Coeffs row_coeffs_accurate(const Airfoil& f, double Re, double ncrit,
                           std::string& src) {
    Coeffs c;
    if (xfoil_lib_coeffs(f, Re, ncrit, c)) { src = "xfoil_lib"; return c; }
    if (native_coeffs(f, Re, ncrit, c))    { src = "native";    return c; }
    src = "synthetic";
    return synthetic_coeffs(f, Re, ncrit);
}

// ---- self-spawning worker process ----------------------------------------
// Invoked as:  build_surrogate(.exe) --worker <jobfile> <resultfile>
// jobfile lines:   idx wu0 wu1 wu2 wu3 wl0 wl1 wl2 wl3 Re ncrit
// resultfile lines: idx wu0..wl3 Re cd0 k cl_max cl_min cm0 src
// Each worker is single-threaded; XFOIL's COMMON blocks are safe within it.
int run_worker(const std::string& jobfile, const std::string& resultfile) {
    std::ifstream in(jobfile);
    if (!in) { std::cerr << "[worker] cannot open " << jobfile << "\n"; return 2; }
    std::ofstream res(resultfile);
    if (!res) { std::cerr << "[worker] cannot open " << resultfile << "\n"; return 2; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        long idx; double s[8], Re, ncrit;
        if (!(ss >> idx)) continue;
        bool ok = true;
        for (int d = 0; d < 8; ++d) ok = ok && static_cast<bool>(ss >> s[d]);
        ok = ok && static_cast<bool>(ss >> Re) && static_cast<bool>(ss >> ncrit);
        if (!ok) continue;
        Airfoil f = shape_to_airfoil(std::vector<double>(s, s + 8), 0.002);
        std::string src;
        Coeffs c = row_coeffs_accurate(f, Re, ncrit, src);
        res << idx;
        for (int d = 0; d < 8; ++d) res << ' ' << s[d];
        res << ' ' << static_cast<long>(Re) << ' ' << c.cd0 << ' ' << c.k << ' '
            << c.cl_max << ' ' << c.cl_min << ' ' << c.cm0 << ' ' << src << "\n";
        res.flush();
    }
    return 0;
}

// Absolute path to this very executable (so workers re-launch the same binary).
std::string self_exe_path(const char* argv0) {
#if defined(_WIN32)
    char buf[MAX_PATH];
    DWORD n = GetModuleFileNameA(nullptr, buf, MAX_PATH);
    if (n > 0 && n < MAX_PATH) return std::string(buf, n);
#endif
    return std::string(argv0);
}

// Launch a worker process; returns an opaque handle to wait on.
#if defined(_WIN32)
intptr_t spawn_worker(const std::string& self, const std::string& jobf,
                      const std::string& resf) {
    const char* args[] = {self.c_str(), "--worker", jobf.c_str(),
                          resf.c_str(), nullptr};
    return _spawnv(_P_NOWAIT, self.c_str(), args);
}
void wait_worker(intptr_t h) { int st = 0; _cwait(&st, h, 0); }
#else
intptr_t spawn_worker(const std::string& self, const std::string& jobf,
                      const std::string& resf) {
    pid_t pid = fork();
    if (pid == 0) {
        execl(self.c_str(), self.c_str(), "--worker", jobf.c_str(),
              resf.c_str(), (char*)nullptr);
        _exit(127);
    }
    return static_cast<intptr_t>(pid);
}
void wait_worker(intptr_t h) { int st = 0; waitpid(static_cast<pid_t>(h), &st, 0); }
#endif
#endif  // HAVE_XFOIL_LIB

}  // namespace

int main(int argc, char** argv) {
#ifdef HAVE_XFOIL_LIB
    // Worker entry point (re-launched by the master for surrogate_mode=xfoil_lib).
    // Bypasses all GA/DoE setup: just crunch a job file and exit.
    if (argc >= 4 && std::string(argv[1]) == "--worker")
        return run_worker(argv[2], argv[3]);
#endif

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

    bool use_xfoil  = (mode == "xfoil");
    bool use_native = (mode == "native");
    if (use_xfoil) {
        // probe: does xfoil run at all?
        if (std::system(("\"" + xfoil_exe + "\" < nul > nul 2>&1").c_str()) != 0)
            std::cerr << "[warn] '" << xfoil_exe
                      << "' did not run; rows that fail fall back to synthetic.\n";
    }

#ifdef HAVE_XFOIL_LIB
    // ---- accurate mode: real XFOIL via a self-spawning worker-process pool ----
    // XFOIL's COMMON state is not thread-safe, so we parallelize across
    // processes: partition the (shape, Re) jobs round-robin into N files, launch
    // N copies of THIS exe with --worker, then merge the per-row results.
    if (mode == "xfoil_lib") {
        const std::string tmp = "build/xfoil_tmp";
        std::filesystem::create_directories(tmp);

        struct Job { long idx; std::vector<double> s; double Re; };
        std::vector<Job> jobs;
        { long idx = 0;
          for (const auto& s : shapes)
              for (double Re : res) jobs.push_back({idx++, s, Re}); }

        int K = static_cast<int>(std::thread::hardware_concurrency());
        if (K < 1) K = 1;
        if (K > static_cast<int>(jobs.size())) K = static_cast<int>(jobs.size());
        if (K > 32) K = 32;

        std::vector<std::string> jpath(K), rpath(K);
        std::vector<std::ofstream> jf;
        jf.reserve(K);
        for (int k = 0; k < K; ++k) {
            jpath[k] = tmp + "/job_" + std::to_string(k) + ".txt";
            rpath[k] = tmp + "/res_" + std::to_string(k) + ".txt";
            jf.emplace_back(jpath[k]);
        }
        for (std::size_t j = 0; j < jobs.size(); ++j) {
            std::ofstream& o = jf[j % K];
            o << jobs[j].idx;
            for (double v : jobs[j].s) o << ' ' << v;
            o << ' ' << static_cast<long>(jobs[j].Re) << ' ' << ncrit << "\n";
        }
        for (auto& o : jf) o.close();

        std::string self = self_exe_path(argv[0]);
        std::cout << "xfoil_lib: spawning " << K << " worker process(es) for "
                  << jobs.size() << " jobs (" << shapes.size() << " shapes x "
                  << res.size() << " Re)...\n";
        std::vector<intptr_t> handles;
        for (int k = 0; k < K; ++k)
            handles.push_back(spawn_worker(self, jpath[k], rpath[k]));
        for (auto h : handles) wait_worker(h);

        // gather results keyed by job index (deterministic CSV order)
        std::vector<std::string> body(jobs.size());
        std::vector<std::string> src(jobs.size());
        for (int k = 0; k < K; ++k) {
            std::ifstream rf(rpath[k]);
            std::string ln;
            while (std::getline(rf, ln)) {
                std::istringstream ss(ln);
                long idx; if (!(ss >> idx)) continue;
                double s8[8], Re, cd0, kk, clmax, clmin, cm0; std::string tag;
                bool ok = true;
                for (int d = 0; d < 8; ++d) ok = ok && static_cast<bool>(ss >> s8[d]);
                ok = ok && static_cast<bool>(ss >> Re)  && static_cast<bool>(ss >> cd0)
                        && static_cast<bool>(ss >> kk)  && static_cast<bool>(ss >> clmax)
                        && static_cast<bool>(ss >> clmin) && static_cast<bool>(ss >> cm0)
                        && static_cast<bool>(ss >> tag);
                if (!ok || idx < 0 || idx >= static_cast<long>(jobs.size())) continue;
                std::ostringstream row;
                row << s8[0];
                for (int d = 1; d < 8; ++d) row << ',' << s8[d];
                row << ',' << static_cast<long>(Re) << ',' << cd0 << ',' << kk
                    << ',' << clmax << ',' << clmin << ',' << cm0;
                body[idx] = row.str();
                src[idx] = tag;
            }
        }

        std::ofstream csv("data/surrogates/polar_coeffs.csv");
        csv << "# wu0,wu1,wu2,wu3,wl0,wl1,wl2,wl3,Re,cd0,k,cl_max,cl_min,cm0\n";
        csv << "# mode=xfoil_lib ncrit=" << ncrit << " shapes=" << shapes.size()
            << " Re=" << res.size() << "\n";
        int n_lib = 0, n_nat = 0, n_syn = 0;
        for (std::size_t i = 0; i < jobs.size(); ++i) {
            if (body[i].empty()) {
                // a worker died or skipped this row: synthesize in the master
                Airfoil f = shape_to_airfoil(jobs[i].s, 0.002);
                Coeffs c = synthetic_coeffs(f, jobs[i].Re, ncrit);
                std::ostringstream row;
                row << jobs[i].s[0];
                for (int d = 1; d < 8; ++d) row << ',' << jobs[i].s[d];
                row << ',' << static_cast<long>(jobs[i].Re) << ',' << c.cd0 << ','
                    << c.k << ',' << c.cl_max << ',' << c.cl_min << ',' << c.cm0;
                body[i] = row.str();
                src[i] = "synthetic";
            }
            csv << body[i] << "\n";
            if (src[i] == "xfoil_lib") ++n_lib;
            else if (src[i] == "native") ++n_nat;
            else ++n_syn;
        }
        csv.close();
        std::cout << "wrote data/surrogates/polar_coeffs.csv: " << jobs.size()
                  << " rows\n  xfoil_lib: " << n_lib << " exact, " << n_nat
                  << " native fallback, " << n_syn << " synthetic fallback\n";
        return 0;
    }
#else
    if (mode == "xfoil_lib") {
        std::cerr << "[warn] surrogate_mode=xfoil_lib requested but this build "
                     "was compiled without HAVE_XFOIL_LIB (no gfortran/libxfoil "
                     "at build time); falling back to native.\n";
        use_native = true;
    }
#endif

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
            else if (use_native)
                ok = native_coeffs(f, Re, ncrit, c);
            if (ok) ++xf_ok;
            else {
                if (use_xfoil || use_native) ++xf_fail;
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
    if (use_xfoil || use_native)
        std::cout << "  " << mode << ": " << xf_ok << " ok, " << xf_fail
                  << " fell back to synthetic\n";
    return 0;
}
