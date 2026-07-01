#include "test_harness.h"
#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/config.h"
#include <cmath>

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
    int wu0 = geom::G_SEC(0,0,0), wl3 = geom::G_SEC(0,1,3);
    double lo0 = spec.lo[wl3], hi0 = spec.hi[wl3];
    seeds::widen_cst_bounds(spec, s, 0.04);
    CHECK(spec.lo[wu0] <= s.airfoils[0].wu[0]);
    CHECK(spec.hi[wu0] >= s.airfoils[0].wu[0]);
    CHECK(spec.lo[wl3] <= lo0 + 1e-12);   // never narrower
    CHECK(spec.hi[wl3] >= hi0 - 1e-12);
}

TEST(seed_genome_decodes_to_seed_shape) {
    seeds::SeedSet s = seeds::load_seeds(naca_only("2412"));
    geom::GenomeSpec spec = geom::default_genome();
    seeds::widen_cst_bounds(spec, s, 0.04);
    auto g = seeds::build_seed_genomes(spec, s, 1, 6);
    CHECK(g.size() == 6);
    CHECK((int)g[0].size() == geom::N_GENES);  // 67 genes
    WingGeometry w = geom::decode(g[0], spec);
    for (int i = 0; i < 4; ++i)
        CHECK_NEAR(w.sections[0].wu[i], s.airfoils[0].wu[i], 1e-6);
}

// Elite seeds put identical CST into all K sections (uniform loft); the GA
// discovers spanwise variation from there.
TEST(seed_all_sections_identical_for_elite) {
    seeds::SeedSet s = seeds::load_seeds(naca_only("2412"));
    geom::GenomeSpec spec = geom::default_genome();
    seeds::widen_cst_bounds(spec, s, 0.04);
    auto g = seeds::build_seed_genomes(spec, s, 1, 6);   // index 0 = elite
    WingGeometry w = geom::decode(g[0], spec);
    for (int k = 1; k < geom::N_SECTIONS; ++k) {
        for (int i = 0; i < 4; ++i) {
            CHECK_NEAR(w.sections[k].wu[i], w.sections[0].wu[i], 1e-9);
            CHECK_NEAR(w.sections[k].wl[i], w.sections[0].wl[i], 1e-9);
        }
    }
}

// Hybrid seeds jitter sections independently -> tip section differs from root.
TEST(seed_sections_jitter_for_hybrid) {
    seeds::SeedSet s = seeds::load_seeds(naca_only("2412"));
    geom::GenomeSpec spec = geom::default_genome();
    seeds::widen_cst_bounds(spec, s, 0.04);
    auto g = seeds::build_seed_genomes(spec, s, 1, 10, 0.3);  // jitter on
    WingGeometry w = geom::decode(g.back(), spec);            // a hybrid (k >= n_elite)
    bool differs = false;
    for (int i = 0; i < 4; ++i)
        if (std::fabs(w.sections[4].wu[i] - w.sections[0].wu[i]) > 1e-6) differs = true;
    CHECK(differs);
}

TEST(seeds_empty_when_none) {
    seeds::SeedSet s = seeds::load_seeds(naca_only(""));
    geom::GenomeSpec spec = geom::default_genome();
    auto g = seeds::build_seed_genomes(spec, s, 1, 10);
    CHECK(s.airfoils.empty());
    CHECK(g.empty());
}
