// panel_fuzz.cpp - robustness sweep of the Morino panel across the gene space.
// The GA reaches corners (high sweep, extreme taper, strong washout/reflex, thin
// TE, low tip Re); the panel must NEVER NaN/Inf/diverge there -- one poisoned
// candidate silently corrupts the Pareto sort (a NaN objective compares false in
// every direction and wrecks non-dominated ranking). This draws thousands of
// random in-bounds genomes, runs the full evaluation chain (decode -> loft ->
// massprops -> panel::solve at several alpha/delta_e -> trim), and asserts every
// output is finite and within sane bounds.
//
// Parallel + deterministic: each genome it gets its own mt19937(seed ^ it), so
// results are reproducible regardless of OpenMP scheduling, and the genomes run
// concurrently (the panel cache is thread_local). Mesh uses nc=6 -- finiteness
// coverage does not need fine chordwise resolution, and the smaller AIC keeps
// the trim Newton (which rebuilds+refactorizes per finite-difference step) fast.
//
// The thing that can actually emit NaN is the panel KERNEL/assembly/solve;
// stability::trim is just a Newton wrapper that calls solve repeatedly, and its
// finite-difference Jacobian rebuilds+refactorizes the AIC on EVERY alpha
// perturbation (the cache key includes wake_alpha) -- a ~50-100x cost multiplier.
// So we sweep thousands of genomes through solve() (the broad finiteness gate)
// and run trim() on only a small subset (to catch the rarer singular-Jacobian /
// NaN-step case). The wake-rebuild-per-FD-step is itself the headline step-6
// perf finding; fixing it (freeze the wake across the Jacobian) would speed up
// both this and the production GA, but is a deliberate solver change, not a fuzz
// shortcut.
//
// Usage: panel_fuzz.exe [N] [margin] [trim_n]
//   N      = genomes through solve() (default 4000)
//   margin = extra widening of the CST/TE bounds to mimic the startup seed-driven
//            widening the GA morphs through (default 0.06).
//   trim_n = how many of the first genomes also go through trim() (default 300).
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/stability.h"
#include "aeroanalyzer/config.h"
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <vector>
#include <string>

using namespace aero;

// Per-genome tally (thread-local; no shared-state contention in the hot path).
struct Tally {
    long checks = 0;
    int fails = 0;
    char first[160] = {0};   // first failing check on this genome
};

static void chk(Tally& t, const char* what, double v, double lo, double hi) {
    ++t.checks;
    if (!(std::isfinite(v) && v >= lo && v <= hi)) {
        if (!t.fails) std::snprintf(t.first, sizeof t.first,
                                    "%s=%.6g (allowed [%.3g,%.3g])", what, v, lo, hi);
        ++t.fails;
    }
}

static void check_state(Tally& t, const AeroState& st) {
    chk(t, "CL", st.CL, -10, 10);
    chk(t, "CDi", st.CDi, 0, 50);
    chk(t, "CDp", st.CDp, 0, 50);
    chk(t, "CD", st.CD, 0, 50);
    chk(t, "e", st.e, 0.0, 1.0001);
    chk(t, "x_np", st.x_np, -1, 1);
    chk(t, "static_margin", st.static_margin, -50, 50);
    chk(t, "CM", st.CM, -100, 100);
    chk(t, "hinge_moment", st.hinge_moment, 0, 1e4);
    chk(t, "x_np_high", st.x_np_high, -1, 1);
    for (double cl : st.cl_local) chk(t, "cl_local", cl, -20, 20);
}

int main(int argc, char** argv) {
    int N = (argc > 1) ? std::atoi(argv[1]) : 4000;
    double margin = (argc > 2) ? std::atof(argv[2]) : 0.06;
    int TRIM_N = (argc > 3) ? std::atoi(argv[3]) : 300;
    int NC = (argc > 4) ? std::atoi(argv[4]) : 6;   // chordwise panels (6=coarse fast, 10=production)

    Config cfg; cfg.set("aero_model", "panel");
    cfg.set("panel_chordwise", std::to_string(NC)); cfg.set("panel_wake_chords", "20");
    viscous::Surrogate surr; surr.load("data/surrogates/polar_coeffs.csv", cfg);

    geom::GenomeSpec spec = geom::default_genome();
    for (int i = geom::G_WU0; i <= geom::G_WL3; ++i) {
        spec.lo[i] -= margin; spec.hi[i] += margin;
    }
    spec.lo[geom::G_TE] = std::max(0.0005, spec.lo[geom::G_TE] - 0.5 * margin);
    spec.hi[geom::G_TE] += 0.5 * margin;

    const double alphas[] = {-4 * DEG2RAD, 0.0, 4 * DEG2RAD, 9 * DEG2RAD};

    long tot_checks = 0;
    int tot_fails = 0, n_trim_fail = 0, n_geom_skip = 0;
    std::vector<double> bad_genes;
    std::string bad_msg;

    #pragma omp parallel for schedule(dynamic, 16) \
        reduction(+:tot_checks, tot_fails, n_trim_fail, n_geom_skip)
    for (int it = 0; it < N; ++it) {
        std::mt19937 rng(0x9e3779b9u ^ static_cast<unsigned>(it));
        std::uniform_real_distribution<double> U(0.0, 1.0);
        std::vector<double> genes(spec.size());
        for (std::size_t k = 0; k < spec.size(); ++k)
            genes[k] = spec.lo[k] + U(rng) * (spec.hi[k] - spec.lo[k]);

        WingGeometry w = geom::decode(genes, spec);
        geom::loft(w, 20);
        MassProps mp = massprops::compute(w, cfg);

        Tally t;
        chk(t, "mass", mp.mass, 1e-4, 50);
        chk(t, "S_ref", mp.S_ref, 1e-6, 10);
        chk(t, "mac", mp.mac, 1e-4, 5);
        chk(t, "AR", mp.AR, 0.1, 100);
        chk(t, "x_cg", mp.x_cg, -1, 1);
        if (t.fails) { ++n_geom_skip; }   // record but still report below

        if (!t.fails) {
            for (double a : alphas)
                check_state(t, potential::solve(w, mp, surr, cfg, a, 0.0));
            check_state(t, potential::solve(w, mp, surr, cfg, 3 * DEG2RAD,  0.15));
            check_state(t, potential::solve(w, mp, surr, cfg, 3 * DEG2RAD, -0.15));
            // trim() only on the first TRIM_N genomes (it is the cost driver).
            if (it < TRIM_N) {
                AeroState tr = stability::trim(w, mp, surr, cfg);
                check_state(t, tr);
                if (!tr.trimmed) ++n_trim_fail;
            }
        }

        tot_checks += t.checks;
        tot_fails  += t.fails;
        if (t.fails) {
            #pragma omp critical
            if (bad_genes.empty()) { bad_genes = genes; bad_msg = t.first; }
        }
    }

    printf("\n=== panel_fuzz: %d genomes, margin=%.3f, nc=%d, %d threads, %ld checks ===\n",
           N, margin, NC, omp_get_max_threads(), tot_checks);
    printf("non-finite / out-of-bounds failures : %d\n", tot_fails);
    printf("geometry-rejected (pre-aero)        : %d\n", n_geom_skip);
    printf("trim did not converge (subset)      : %d / %d\n", n_trim_fail, TRIM_N);
    if (tot_fails && !bad_genes.empty()) {
        printf("first failure: %s\n  genome:", bad_msg.c_str());
        for (double v : bad_genes) printf(" %.4f", v);
        printf("\n");
    }
    printf("%s\n", tot_fails == 0 ? "RESULT: PASS (all finite/bounded)"
                                  : "RESULT: FAIL");
    return tot_fails ? 1 : 0;
}
