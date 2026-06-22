// seeds.h — ancestry injection (plan §3): build the GA's initial DNA from real
// airfoils. Loads uploaded .dat seeds + generated NACA, fits them to CST, widens
// the CST gene bounds to contain them, and emits seed genomes for NSGA-II.
#pragma once
#include <vector>
#include <string>
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/config.h"

namespace aero {
namespace seeds {

struct SeedSet {
    std::vector<Airfoil> airfoils;     // fitted at order 3 (4+4 weights)
    std::vector<std::string> names;
};

// Load .dat from cfg "seed_airfoils_dir" + NACA from cfg "seed_naca".
SeedSet load_seeds(const Config& cfg);

// Expand CST gene bounds (G_WU0..3, G_WL0..3, G_TE) in place so every seed fits,
// with the given absolute margin for morph room. No-op if no seeds.
void widen_cst_bounds(geom::GenomeSpec& spec, const SeedSet& s, double margin);

// Build up to `count` seed genomes: the first ~20% are "elites" (seed shape +
// mid-box planform), the rest "hybrids" (seed shape + randomized planform). The
// GA fills the remaining population with random "explorers".
std::vector<std::vector<double>> build_seed_genomes(
    const geom::GenomeSpec& spec, const SeedSet& s, unsigned rng_seed, int count);

}  // namespace seeds
}  // namespace aero
