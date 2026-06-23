// xfoil_cyl_probe.cpp -- validate Hess-Smith core on a cylinder (Cp=1-4sin^2).
#include "aeroanalyzer/aero_xfoil.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace aero;

int main() {
    const double PI = 3.14159265358979;
    int n = 100;
    std::vector<double> X, Y;
    // CCW circle, radius 0.5, centered origin. Start at theta=0.
    for (int i = 0; i <= n; ++i) {
        double th = 2 * PI * i / n;
        X.push_back(0.5 * std::cos(th));
        Y.push_back(0.5 * std::sin(th));
    }
    xfoil::Solver s;
    if (!s.set_coords(X, Y)) { std::printf("set_coords FAILED\n"); return 1; }
    xfoil::Result r = s.solve(0.0, 1e6);
    std::printf("cylinder: cl=%.5f (expect 0), cm=%.5f (expect 0)\n", r.cl, r.cm);
    return 0;
}
