#include "test_harness.h"
#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/config.h"

using namespace aero;

static Config naca_only(const std::string& list) {
    Config c;
    c.set("seed_airfoils_dir", "___no_such_dir___");
    c.set("seed_naca", list);
    return c;
}

TEST(seeds_naca_loaded) {
    seeds::SeedSet s = seeds::load_seeds(naca_only("0012,2412"));
    CHECK(s.airfoils.size() == 2);
    CHECK(s.airfoils[0].wu.size() == 4);
    CHECK(s.airfoils[0].wl.size() == 4);
}

TEST(seeds_widen_bounds_contains_seed) {
    seeds::SeedSet s = seeds::load_seeds(naca_only("4412"));
    geom::GenomeSpec spec = geom::default_genome();
    double lo0 = spec.lo[geom::G_WL3], hi0 = spec.hi[geom::G_WL3];
    seeds::widen_cst_bounds(spec, s, 0.04);
    CHECK(spec.lo[geom::G_WU0] <= s.airfoils[0].wu[0]);
    CHECK(spec.hi[geom::G_WU0] >= s.airfoils[0].wu[0]);
    CHECK(spec.lo[geom::G_WL3] <= lo0 + 1e-12);   // never narrower
    CHECK(spec.hi[geom::G_WL3] >= hi0 - 1e-12);
}

TEST(seed_genome_decodes_to_seed_shape) {
    seeds::SeedSet s = seeds::load_seeds(naca_only("2412"));
    geom::GenomeSpec spec = geom::default_genome();
    seeds::widen_cst_bounds(spec, s, 0.04);
    auto g = seeds::build_seed_genomes(spec, s, 1, 6);
    CHECK(g.size() == 6);
    WingGeometry w = geom::decode(g[0], spec);
    for (int i = 0; i < 4; ++i)
        CHECK_NEAR(w.section.wu[i], s.airfoils[0].wu[i], 1e-6);
}

TEST(seeds_empty_when_none) {
    seeds::SeedSet s = seeds::load_seeds(naca_only(""));
    geom::GenomeSpec spec = geom::default_genome();
    auto g = seeds::build_seed_genomes(spec, s, 1, 10);
    CHECK(s.airfoils.empty());
    CHECK(g.empty());
}
