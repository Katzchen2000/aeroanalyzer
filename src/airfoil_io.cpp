#include "aeroanalyzer/airfoil_io.h"
#include "aeroanalyzer/geom.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace aero {
namespace airfoil_io {

namespace {
bool parse_two(const std::string& line, double& a, double& b) {
    std::istringstream ss(line);
    return static_cast<bool>(ss >> a >> b);
}

// Normalize so x spans [0,1] (chord-normalize defensively).
void normalize(std::vector<std::pair<double, double>>& pts, double xmin,
               double scale) {
    for (auto& p : pts) {
        p.first = (p.first - xmin) * scale;
        p.second = p.second * scale;
    }
}
}  // namespace

Coords load_dat(const std::string& path, bool& ok) {
    ok = false;
    Coords c;
    std::ifstream f(path);
    if (!f) return c;

    std::string line;
    std::vector<std::pair<double, double>> raw;
    bool first = true;
    while (std::getline(f, line)) {
        double a, b;
        if (!parse_two(line, a, b)) {
            if (first) c.name = line;   // title line
            first = false;
            continue;
        }
        first = false;
        raw.push_back({a, b});
    }
    if (raw.size() < 5) return c;

    // Lednicer header: first numeric pair is a pair of point counts (> 1.5).
    bool lednicer = (raw[0].first > 1.5 && raw[0].second > 1.5);
    if (lednicer) {
        int nu = static_cast<int>(std::lround(raw[0].first));
        int nl = static_cast<int>(std::lround(raw[0].second));
        if (nu < 2 || nl < 2 || 1 + nu + nl > static_cast<int>(raw.size()))
            return c;
        c.upper.assign(raw.begin() + 1, raw.begin() + 1 + nu);
        c.lower.assign(raw.begin() + 1 + nu, raw.begin() + 1 + nu + nl);
    } else {
        // Selig: TE(upper) -> LE -> TE(lower). Split at min-x.
        std::size_t imin = 0;
        for (std::size_t i = 1; i < raw.size(); ++i)
            if (raw[i].first < raw[imin].first) imin = i;
        c.upper.assign(raw.begin(), raw.begin() + imin + 1);
        std::reverse(c.upper.begin(), c.upper.end());  // -> LE..TE
        c.lower.assign(raw.begin() + imin, raw.end());
    }

    // Ensure ascending x and chord-normalize to [0,1].
    auto byx = [](const std::pair<double, double>& a,
                  const std::pair<double, double>& b) { return a.first < b.first; };
    std::sort(c.upper.begin(), c.upper.end(), byx);
    std::sort(c.lower.begin(), c.lower.end(), byx);
    double xmin = std::min(c.upper.front().first, c.lower.front().first);
    double xmax = std::max(c.upper.back().first, c.lower.back().first);
    if (xmax - xmin < 1e-9) return c;
    double scale = 1.0 / (xmax - xmin);
    normalize(c.upper, xmin, scale);
    normalize(c.lower, xmin, scale);
    ok = true;
    return c;
}

Coords naca4(const std::string& code, int n) {
    Coords c;
    c.name = "NACA " + code;
    double m = 0, p = 0, t = 0.12;
    if (code.size() >= 4) {
        m = (code[0] - '0') / 100.0;
        p = (code[1] - '0') / 10.0;
        t = ((code[2] - '0') * 10 + (code[3] - '0')) / 100.0;
    }
    for (int i = 0; i < n; ++i) {
        double x = 0.5 * (1.0 - std::cos(PI * i / (n - 1)));
        double yt = 5.0 * t *
                    (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                     0.2843 * x * x * x - 0.1015 * x * x * x * x);
        double yc = 0.0;
        if (p > 1e-6 && m > 1e-9) {
            if (x < p) yc = m / (p * p) * (2.0 * p * x - x * x);
            else yc = m / ((1.0 - p) * (1.0 - p)) *
                      ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
        }
        c.upper.push_back({x, yc + yt});
        c.lower.push_back({x, yc - yt});
    }
    return c;
}

double estimate_te(const Coords& c) {
    if (c.upper.empty() || c.lower.empty()) return 0.0;
    double te = c.upper.back().second - c.lower.back().second;
    return std::max(0.0, te);
}

Airfoil to_airfoil(const Coords& c, int order, double te) {
    Airfoil f = geom::fit_cst(c.upper, c.lower, order, te);
    f.te_thick = te;
    return f;
}

}  // namespace airfoil_io
}  // namespace aero
