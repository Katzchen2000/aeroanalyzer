// xfoil_visc_probe.cpp -- exercise the native viscous BL (Squire-Young cd) on
// NACA 0012 and print a small polar to eyeball against XFOIL.
#include "aeroanalyzer/aero_xfoil.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace aero;

static void naca0012(std::vector<double>& X, std::vector<double>& Y, int nside) {
    const double t = 0.12, PI = 3.14159265358979;
    auto yt = [&](double x) {
        return 5 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                        0.2843 * x * x * x - 0.1036 * x * x * x * x);
    };
    X.clear(); Y.clear();
    for (int j = nside; j >= 1; --j) {
        double x = 0.5 * (1 - std::cos(PI * j / nside));
        X.push_back(x); Y.push_back(yt(x));
    }
    X.push_back(0.0); Y.push_back(0.0);
    for (int j = 1; j <= nside; ++j) {
        double x = 0.5 * (1 - std::cos(PI * j / nside));
        X.push_back(x); Y.push_back(-yt(x));
    }
}

int main(int argc, char** argv) {
    double Re = (argc > 1) ? std::atof(argv[1]) : 2.0e5;
    std::vector<double> X, Y;
    naca0012(X, Y, 100);
    xfoil::Options opt;
    opt.Ncrit = 9.0;
    xfoil::Solver s(opt);
    if (!s.set_coords(X, Y)) { std::printf("set_coords FAILED\n"); return 1; }

    std::printf("NACA0012  Re=%.0f  Ncrit=%.1f\n", Re, opt.Ncrit);
    std::printf("  alpha     cl       cd       cdf      cdp     xtr_t  xtr_b  conv sep\n");
    const double d2r = 3.14159265358979 / 180.0;
    for (double adeg = 0.0; adeg <= 8.001; adeg += 1.0) {
        xfoil::Result r = s.solve(adeg * d2r, Re);
        std::printf("  %5.1f  %7.4f  %7.5f  %7.5f  %7.5f  %5.3f  %5.3f   %d   %d\n",
                    adeg, r.cl, r.cd, r.cdf, r.cdp, r.xtr_top, r.xtr_bot,
                    r.converged ? 1 : 0, r.separated ? 1 : 0);
    }
    return 0;
}
