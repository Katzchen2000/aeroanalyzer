// evaluate.h — the candidate fitness pipeline (phases 1→5 -> objectives+cv).
#pragma once
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/seeds.h"
#include "aeroanalyzer/config.h"

namespace aero {

// Full per-candidate detail, kept for reporting / AVL export of the incumbents.
struct EvalResult {
    WingGeometry geom;
    MassProps mp;
    AeroState aero;
    std::array<double, N_OBJ> objectives{{0, 0, 0}};
    double cv = 0.0;
};

class Evaluator {
public:
    explicit Evaluator(const Config& cfg);

    // Fills c.objectives and c.cv (thread-safe: const, surrogate is read-only).
    void evaluate(Candidate& c) const;

    EvalResult detail(const std::vector<double>& genes) const;

    const geom::GenomeSpec& spec() const { return spec_; }
    const seeds::SeedSet& seeds() const { return seeds_; }

private:
    EvalResult run(const std::vector<double>& genes) const;

    const Config& cfg_;
    geom::GenomeSpec spec_;
    seeds::SeedSet seeds_;
    viscous::Surrogate surr_;
    int n_stations_;
};

}  // namespace aero
