// panel_timing.cpp - step-6 microbenchmark + accuracy + determinism.
// Does NOT run a GA for timing (that conflates GA overhead with solver cost).
// Each GA candidate has fresh geometry, so the realistic per-candidate aero cost
// is ONE cold stability::trim() (the evaluator's only aero call). We therefore
// time trim() on FRESH wings (each a cold cache miss, no preceding solve to warm
// it), for three wake-freeze strategies, and project to pop120 x gen150 = 18000
// candidates arithmetically:
//   legacy : rebuild the dense N*N AIC on every residual solve AND every FD probe
//   jac    : freeze the wake across the FD Jacobian probes only
//   full   : freeze the wake for the whole Newton solve (1 build / candidate)
// It also reports the trim-point DIFFERENCE (full vs legacy) so the accuracy cost
// of freezing is visible, and checks panel determinism (1 vs N threads) with a
// tiny GA (the panel cache is thread_local keyed on geom_sig_panel).
//
// Usage: panel_timing.exe [n_wings] [nc]
#include "aeroanalyzer/config.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/stability.h"
#include "aeroanalyzer/evaluate.h"
#include "aeroanalyzer/ga.h"
#include "aeroanalyzer/seeds.h"
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

using namespace aero;
using Clock = std::chrono::steady_clock;

static double ms_since(Clock::time_point t0) {
    return 1e3 * std::chrono::duration<double>(Clock::now() - t0).count();
}
static double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// A spread of distinct in-bounds wings so every solve is a cold cache miss.
static std::vector<WingGeometry> make_wings(int n) {
    geom::GenomeSpec spec = geom::default_genome();
    std::vector<WingGeometry> out;
    std::mt19937 rng(12345u);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        std::vector<double> g(spec.size());
        for (std::size_t k = 0; k < spec.size(); ++k)
            g[k] = spec.lo[k] + U(rng) * (spec.hi[k] - spec.lo[k]);
        WingGeometry w = geom::decode(g, spec);
        geom::loft(w, 20);
        out.push_back(w);
    }
    return out;
}

struct TrimRun { double median_ms; std::vector<AeroState> st; };

// Time trim() on each FRESH wing (no preceding solve -> realistic cold candidate).
static TrimRun trim_all(const std::vector<WingGeometry>& wings,
                        const viscous::Surrogate& surr, const Config& cfg) {
    TrimRun r;
    std::vector<double> t;
    for (const auto& w : wings) {
        MassProps mp = massprops::compute(w, cfg);
        auto t0 = Clock::now();
        AeroState s = stability::trim(w, mp, surr, cfg);
        t.push_back(ms_since(t0));
        r.st.push_back(s);
    }
    r.median_ms = median(t);
    return r;
}

static double cold_solve_ms(const std::vector<WingGeometry>& wings,
                            const viscous::Surrogate& surr, const Config& cfg) {
    std::vector<double> t;
    for (const auto& w : wings) {
        MassProps mp = massprops::compute(w, cfg);
        auto t0 = Clock::now();
        AeroState s = potential::solve(w, mp, surr, cfg, 3.0 * DEG2RAD, 0.0);
        t.push_back(ms_since(t0)); (void)s;
    }
    return median(t);
}

// Largest |difference| in the trim outputs that feed GA objectives/constraints.
static void state_diff(const std::vector<AeroState>& a, const std::vector<AeroState>& b,
                       double& dCL, double& dCM, double& dxnp, double& dsm, int& dtrim) {
    dCL = dCM = dxnp = dsm = 0.0; dtrim = 0;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) {
        dCL  = std::max(dCL,  std::fabs(a[i].CL - b[i].CL));
        dCM  = std::max(dCM,  std::fabs(a[i].CM - b[i].CM));
        dxnp = std::max(dxnp, std::fabs(a[i].x_np - b[i].x_np));
        dsm  = std::max(dsm,  std::fabs(a[i].static_margin - b[i].static_margin));
        if (a[i].trimmed != b[i].trimmed) ++dtrim;
    }
}

static int front0(const std::vector<Candidate>& p) {
    int n = 0; for (const auto& c : p) if (c.rank == 0 && c.cv <= 0.0) ++n; return n;
}

int main(int argc, char** argv) {
    int NW = (argc > 1) ? std::atoi(argv[1]) : 24;
    int nc = (argc > 2) ? std::atoi(argv[2]) : 10;

    Config base;
    if (!base.load("config/baseline.cfg")) { printf("no baseline.cfg\n"); return 1; }
    base.set("panel_chordwise", std::to_string(nc));
    base.set("aero_model", "panel");
    viscous::Surrogate surr; surr.load("data/surrogates/polar_coeffs.csv", base);
    int maxth = omp_get_max_threads();
    omp_set_num_threads(1);   // time the SOLVER, not OpenMP fan-out

    std::vector<WingGeometry> wings = make_wings(NW);

    // VLM reference (no wake cache; freeze knobs inert).
    Config vlm = base; vlm.set("aero_model", "vlm");
    double vlm_trim = trim_all(wings, surr, vlm).median_ms;

    // Three panel freeze strategies.
    Config legacy = base; legacy.set("panel_trim_freeze_jac", "0");
                          legacy.set("panel_trim_freeze_wake", "0");
    Config jac    = base; jac.set("panel_trim_freeze_jac", "1");
                          jac.set("panel_trim_freeze_wake", "0");
    Config full   = base; full.set("panel_trim_freeze_jac", "1");
                          full.set("panel_trim_freeze_wake", "1");

    double cold = cold_solve_ms(wings, surr, base);
    TrimRun rL = trim_all(wings, surr, legacy);
    TrimRun rJ = trim_all(wings, surr, jac);
    TrimRun rF = trim_all(wings, surr, full);

    printf("=== panel microbenchmark (nc=%d, %d fresh wings, single-thread) ===\n", nc, NW);
    printf("  cold panel solve (1 AIC build) : %7.2f ms\n", cold);
    printf("  vlm   trim                     : %7.2f ms\n", vlm_trim);
    printf("  panel trim  legacy             : %7.2f ms\n", rL.median_ms);
    printf("  panel trim  jac-freeze         : %7.2f ms  (%.2fx vs legacy)\n",
           rJ.median_ms, rL.median_ms / std::max(1e-9, rJ.median_ms));
    printf("  panel trim  full-freeze        : %7.2f ms  (%.2fx vs legacy)\n",
           rF.median_ms, rL.median_ms / std::max(1e-9, rF.median_ms));

    double dCL, dCM, dxnp, dsm; int dtrim;
    state_diff(rL.st, rF.st, dCL, dCM, dxnp, dsm, dtrim);
    printf("  accuracy  full-freeze vs legacy: dCL=%.2e dCM=%.2e dx_np=%.2e dSM=%.2e trim_flag_diffs=%d\n",
           dCL, dCM, dxnp, dsm, dtrim);

    double cand = 18000.0;
    auto proj = [&](double ms) { return ms * cand / 1e3 / maxth; };
    printf("  extrapolated pop120 x gen150 (=18000 trims) at %d threads:\n", maxth);
    printf("    vlm                 ~ %7.1f s (%.1f min)\n", proj(vlm_trim), proj(vlm_trim)/60);
    printf("    panel legacy        ~ %7.1f s (%.1f min)\n", proj(rL.median_ms), proj(rL.median_ms)/60);
    printf("    panel full-freeze   ~ %7.1f s (%.1f min)\n", proj(rF.median_ms), proj(rF.median_ms)/60);

    // ---- determinism: tiny panel GA, serial vs OpenMP, byte-identical ----
    auto run_ga = [&](int threads) {
        Config cfg = base;
        omp_set_num_threads(threads);
        Evaluator eval(cfg);
        ga::Params p; p.pop = 12; p.generations = 4; p.seed = 1;
        ga::NSGA2 opt(eval, p);
        opt.set_seeds(seeds::build_seed_genomes(eval.spec(), eval.seeds(), p.seed, 6));
        return opt.run();
    };
    auto a = run_ga(1);
    auto b = run_ga(maxth);
    double md = 0.0; bool same = (a.size() == b.size());
    for (std::size_t i = 0; same && i < a.size(); ++i) {
        for (std::size_t k = 0; k < a[i].genes.size(); ++k)
            md = std::max(md, std::fabs(a[i].genes[k] - b[i].genes[k]));
        for (int o = 0; o < N_OBJ; ++o)
            md = std::max(md, std::fabs(a[i].objectives[o] - b[i].objectives[o]));
        md = std::max(md, std::fabs(a[i].cv - b[i].cv));
    }
    printf("\n=== determinism (panel, seed=1, 1-thread vs %d-thread) ===\n", maxth);
    printf("  front0(1t)=%d front0(%dt)=%d  max|diff|=%.3g  -> %s\n",
           front0(a), maxth, front0(b), md,
           (same && md == 0.0) ? "IDENTICAL" : "DIVERGENT");
    return (same && md == 0.0) ? 0 : 2;
}
