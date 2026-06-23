// xfoil_inviscid_probe.cpp -- validate the inviscid Hess-Smith stage of
// aero::xfoil against thin-airfoil theory. Build:
//   g++ -std=c++17 -O2 -I include -I <eigen> scratch/xfoil_inviscid_probe.cpp \
//       src/aero_xfoil.cpp src/geom.cpp src/config.cpp src/airfoil_io.cpp \
//       src/massprops.cpp -o build/xfoil_inviscid_probe.exe
#include "aeroanalyzer/aero_xfoil.h"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace aero;

// NACA 0012 ordered TE -> upper -> LE -> lower -> TE (closed-ish TE).
static void naca0012(std::vector<double>& X, std::vector<double>& Y, int nside) {
    const double t = 0.12, PI = 3.14159265358979;
    auto yt = [&](double x) {   // closed-TE coefficient (-0.1036) -> sharp TE
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

int main() {
    std::vector<double> X, Y;
    naca0012(X, Y, 80);
    xfoil::Options opt;
    xfoil::Solver s(opt);
    if (!s.set_coords(X, Y)) { std::printf("set_coords FAILED\n"); return 1; }

    std::printf("NACA0012 inviscid (Hess-Smith), %zu nodes\n", X.size());
    std::printf("  alpha    cl        cm        dcl/dalpha[/deg]\n");
    double prev_cl = 0, prev_a = 0;
    const double d2r = 3.14159265358979 / 180.0;
    for (double adeg = 0; adeg <= 8.001; adeg += 2.0) {
        xfoil::Result r = s.solve(adeg * d2r, 1e6);
        double slope = (adeg > 0) ? (r.cl - prev_cl) / (adeg - prev_a) : 0.0;
        std::printf("  %5.1f  %8.4f  %8.4f   %8.4f\n", adeg, r.cl, r.cm, slope);
        prev_cl = r.cl; prev_a = adeg;
    }
    std::printf("expect: cl slope ~0.10-0.12 /deg (2*pi/rad), cm ~ 0\n");
    return 0;
}
