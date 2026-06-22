#include "test_harness.h"
#include "aeroanalyzer/airfoil_io.h"
#include "aeroanalyzer/geom.h"
#include <filesystem>
#include <fstream>

using namespace aero;

// NACA 0012: ~12% max thickness, symmetric (zero camber).
TEST(naca0012_thickness) {
    airfoil_io::Coords c = airfoil_io::naca4("0012", 200);
    double tmax = 0.0, cam = 0.0;
    for (std::size_t i = 0; i < c.upper.size(); ++i) {
        tmax = std::max(tmax, c.upper[i].second - c.lower[i].second);
        cam = std::max(cam, std::fabs(c.upper[i].second + c.lower[i].second));
    }
    CHECK_NEAR(tmax, 0.12, 3e-3);
    CHECK_NEAR(cam, 0.0, 1e-9);    // symmetric
}

// NACA 2412: positive camber => negative zero-lift angle after CST fit.
TEST(naca2412_cambered) {
    airfoil_io::Coords c = airfoil_io::naca4("2412", 120);
    Airfoil f = airfoil_io::to_airfoil(c, 3, 0.0);
    geom::ThinAirfoil t = geom::thin_airfoil(f);
    CHECK(t.alpha_L0 < 0.0);
}

// to_airfoil round-trips the surface within fit tolerance.
TEST(naca_to_cst_round_trip) {
    airfoil_io::Coords c = airfoil_io::naca4("2412", 120);
    Airfoil f = airfoil_io::to_airfoil(c, 4, airfoil_io::estimate_te(c));
    for (double x = 0.1; x < 0.95; x += 0.2) {
        // find nearest sampled point
        double zu = 0, best = 1e9;
        for (auto& p : c.upper)
            if (std::fabs(p.first - x) < best) { best = std::fabs(p.first - x); zu = p.second; }
        CHECK_NEAR(geom::cst_upper(f, x), zu, 5e-3);
    }
}

static std::string tmp(const std::string& name) {
    return (std::filesystem::temp_directory_path() / name).string();
}

// Selig format: TE(upper)->LE->TE(lower) single loop.
TEST(load_selig) {
    airfoil_io::Coords n = airfoil_io::naca4("0012", 30);
    std::string path = tmp("aero_selig.dat");
    {
        std::ofstream o(path);
        o << "TEST SELIG\n";
        for (auto it = n.upper.rbegin(); it != n.upper.rend(); ++it)
            o << it->first << " " << it->second << "\n";          // TE->LE
        for (std::size_t i = 1; i < n.lower.size(); ++i)
            o << n.lower[i].first << " " << n.lower[i].second << "\n";  // LE->TE
    }
    bool ok = false;
    airfoil_io::Coords c = airfoil_io::load_dat(path, ok);
    CHECK(ok);
    CHECK_NEAR(c.upper.front().first, 0.0, 1e-6);
    CHECK_NEAR(c.upper.back().first, 1.0, 1e-6);
    CHECK(c.lower.size() >= 2);
}

// Lednicer format: "nu nl" header, then upper LE->TE, then lower LE->TE.
TEST(load_lednicer) {
    airfoil_io::Coords n = airfoil_io::naca4("0012", 20);
    std::string path = tmp("aero_lednicer.dat");
    {
        std::ofstream o(path);
        o << "TEST LEDNICER\n";
        o << n.upper.size() << " " << n.lower.size() << "\n";
        for (auto& p : n.upper) o << p.first << " " << p.second << "\n";
        for (auto& p : n.lower) o << p.first << " " << p.second << "\n";
    }
    bool ok = false;
    airfoil_io::Coords c = airfoil_io::load_dat(path, ok);
    CHECK(ok);
    CHECK(c.upper.size() == n.upper.size());
    CHECK(c.lower.size() == n.lower.size());
    CHECK_NEAR(c.upper.back().first, 1.0, 1e-6);
}
