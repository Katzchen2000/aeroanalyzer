// xfoil_xcheck.cpp -- OFFLINE oracle cross-check: tabulate the native aero::xfoil
// solver against the real xfoil.exe for a few standard sections across the Re
// grid, and report the max cd/cl deviation in the unstalled band. One-off; mirrors
// the AVL cross-check probes. Build like the other scratch probes, then run from
// the project root (needs tools/bin/xfoil.exe). Drives xfoil via a temp command
// file -- this is the ONLY path here that touches disk, and only for the oracle.
#include "aeroanalyzer/aero_xfoil.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace aero;

namespace {
const double D2R = 3.14159265358979 / 180.0;

void naca(std::vector<double>& X, std::vector<double>& Y, int nside,
          double t, double m, double p) {
    const double PI = 3.14159265358979;
    auto camber = [&](double x, double& yc, double& dy) {
        if (m == 0.0 || p == 0.0) { yc = 0; dy = 0; return; }
        if (x < p) { yc = m / (p * p) * (2 * p * x - x * x);
                     dy = 2 * m / (p * p) * (p - x); }
        else { yc = m / ((1 - p) * (1 - p)) * (1 - 2 * p + 2 * p * x - x * x);
               dy = 2 * m / ((1 - p) * (1 - p)) * (p - x); }
    };
    auto yt = [&](double x) {
        return 5 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                        0.2843 * x * x * x - 0.1015 * x * x * x * x);
    };
    X.clear(); Y.clear();
    auto pt = [&](double x, int up) {
        double yc, dy, th = std::atan(dy); (void)th;
        camber(x, yc, dy); double ang = std::atan(dy);
        double yth = yt(x);
        double xx = x - up * yth * std::sin(ang);
        double yy = yc + up * yth * std::cos(ang);
        X.push_back(xx); Y.push_back(yy);
    };
    for (int j = nside; j >= 1; --j) pt(0.5 * (1 - std::cos(PI * j / nside)), +1);
    X.push_back(0); Y.push_back(0);
    for (int j = 1; j <= nside; ++j) pt(0.5 * (1 - std::cos(PI * j / nside)), -1);
}

// Drive real xfoil.exe for one (section file, Re); parse alpha,CL,CD.
bool run_xfoil(const std::string& dat, double Re, double ncrit,
               std::vector<double>& al, std::vector<double>& cl,
               std::vector<double>& cd) {
    std::string pol = "build\\xcheck.pol", cmd = "build\\xcheck.cmd";
    std::remove(pol.c_str());
    { std::ofstream o(cmd);
      o << "LOAD " << dat << "\nPANE\nOPER\nVPAR\nN " << ncrit
        << "\n\nVISC " << (long)Re << "\nITER 220\nPACC\n" << pol
        << "\n\nASEQ 0 8 1\nPACC\n\nQUIT\n"; }
    // cmd.exe strips one outer quote pair, so wrap the whole line in quotes.
    std::string run = "\"tools\\bin\\xfoil.exe < " + cmd +
                      " > build\\xcheck.log 2>&1\"";
    std::system(run.c_str());
    std::ifstream f(pol); if (!f) return false;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line); std::vector<double> v; double x;
        while (ss >> x) v.push_back(x);
        if (v.size() >= 5 && v[0] >= -20 && v[0] <= 30)
            { al.push_back(v[0]); cl.push_back(v[1]); cd.push_back(v[2]); }
    }
    return al.size() >= 3;
}

void write_dat(const std::string& path, const std::vector<double>& X,
               const std::vector<double>& Y) {
    std::ofstream o(path); o << "section\n";
    for (std::size_t i = 0; i < X.size(); ++i) o << X[i] << " " << Y[i] << "\n";
}
}  // namespace

int main() {
    struct Sec { const char* name; double t, m, p; };
    Sec secs[] = {{"NACA0012", 0.12, 0.0, 0.0},
                  {"NACA2412", 0.12, 0.02, 0.4},
                  {"NACA4412", 0.12, 0.04, 0.4}};
    double Re = 2.0e5, ncrit = 9.0;
    std::system("mkdir build 2>nul");

    for (auto& sc : secs) {
        std::vector<double> X, Y;
        naca(X, Y, 100, sc.t, sc.m, sc.p);
        std::string dat = std::string("build\\") + sc.name + ".dat";
        write_dat(dat, X, Y);

        xfoil::Options opt; opt.Ncrit = ncrit;
        xfoil::Solver s(opt);
        if (!s.set_coords(X, Y)) { std::printf("%s: native geom FAIL\n", sc.name); continue; }

        std::vector<double> xal, xcl, xcd;
        bool xok = run_xfoil(dat, Re, ncrit, xal, xcl, xcd);

        std::printf("\n%s  Re=%.0f  Ncrit=%.1f\n", sc.name, Re, ncrit);
        std::printf("  a     cl(nat) cl(xf)   cd(nat)  cd(xf)   dcd%%\n");
        double maxdcd = 0.0;
        for (double a = 0; a <= 8.001; a += 1.0) {
            xfoil::Result r = s.solve(a * D2R, Re);
            double cxf = -1, lxf = -1;
            if (xok) for (std::size_t i = 0; i < xal.size(); ++i)
                if (std::fabs(xal[i] - a) < 0.4) { cxf = xcd[i]; lxf = xcl[i]; }
            double dcd = (cxf > 0 && r.converged) ? 100.0 * (r.cd - cxf) / cxf : 0.0;
            if (cxf > 0 && r.converged) maxdcd = std::max(maxdcd, std::fabs(dcd));
            std::printf("  %4.1f  %6.3f  %6.3f   %7.5f  %7.5f  %+6.1f%s\n",
                        a, r.cl, lxf, r.cd, cxf, dcd, r.converged ? "" : "  (nat n/c)");
        }
        std::printf("  max |dcd| in band: %.1f%%\n", maxdcd);
    }
    return 0;
}
