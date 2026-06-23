#include "test_harness.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <filesystem>
#include <fstream>

using namespace aero;

TEST(surrogate_tables_load_and_query) {
    std::string dir =
        (std::filesystem::temp_directory_path() / "aero_surr_test").string();
    std::filesystem::create_directories(dir);
    {
        std::ofstream o(dir + "/polar_coeffs.csv");
        o << "# wu0..wl3,Re,cd0,k,cl_max,cl_min,cm0\n";
        o << "0.20,0.17,0.14,0.11,-0.12,-0.09,-0.02,0.06,200000,"
             "0.010,0.013,1.20,-0.60,0.00\n";
        o << "0.25,0.20,0.16,0.12,-0.10,-0.07,0.00,0.08,200000,"
             "0.012,0.014,1.25,-0.60,0.02\n";
    }
    Config cfg;
    cfg.set("viscous_backend", "table");   // exercise the legacy IDW path
    viscous::Surrogate surr;
    bool ok = surr.load(dir, cfg);
    CHECK(ok);
    CHECK(surr.using_tables());
    CHECK(surr.sample_count() == 2);

    std::vector<double> shape = {0.20, 0.17, 0.14, 0.11, -0.12, -0.09, -0.02, 0.06};
    viscous::Polar p = surr.query(shape, 0.4, 200000);
    CHECK(!p.clamped);
    CHECK(p.cd > 0.0);
    CHECK_NEAR(p.cd, 0.010 + 0.013 * 0.16, 2e-3);   // ~ sample 1's polar

    // shape outside the trained hull -> flagged
    std::vector<double> far(8, 1.0);
    CHECK(surr.query(far, 0.3, 200000).clamped);
    // absurd cl beyond cl_max -> flagged
    CHECK(surr.query(shape, 3.0, 200000).clamped);
}

TEST(surrogate_missing_table_falls_back) {
    Config cfg;
    cfg.set("viscous_backend", "table");   // no weights, no CSV -> analytic
    viscous::Surrogate surr;
    bool ok = surr.load("___no_such_surrogate_dir___", cfg);
    CHECK(!ok);
    CHECK(!surr.using_tables());
    // still usable via analytic fallback
    std::vector<double> shape = {0.20, 0.17, 0.14, 0.11, -0.12, -0.09, -0.02, 0.06};
    CHECK(surr.query(shape, 0.3, 200000).cd > 0.0);
}
