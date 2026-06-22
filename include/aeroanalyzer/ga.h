// ga.h — NSGA-II with Deb's constraint-domination (drag / mass / SM front).
#pragma once
#include <vector>
#include <functional>
#include "aeroanalyzer/engine_core.h"
#include "aeroanalyzer/evaluate.h"

namespace aero {
namespace ga {

struct Params {
    int pop = 120;
    int generations = 150;
    unsigned seed = 1;
    double eta_cx = 15.0;
    double eta_mut = 20.0;
    double pcx = 0.9;
    double pmut = -1.0;   // <0 => 1/n_genes
};

// Constraint-domination (Deb 2002): feasible beats infeasible; among infeasible
// the smaller violation wins; among feasible, Pareto dominance on objectives.
bool constrained_dominates(const Candidate& a, const Candidate& b);

// Ranks the population into fronts (sets rank on each candidate).
void fast_nondominated_sort(std::vector<Candidate>& pop,
                            std::vector<std::vector<int>>& fronts);

// Sets crowding distance on the candidates listed in `front`.
void crowding_distance(std::vector<Candidate>& pop, const std::vector<int>& front);

class NSGA2 {
public:
    NSGA2(const Evaluator& eval, const Params& p);

    // Inject seed genomes (ancestry injection). The first seeds.size() members
    // of the initial population are taken from these (clamped to bounds); the
    // rest are random explorers.
    void set_seeds(std::vector<std::vector<double>> s) { seed_genes_ = std::move(s); }

    // Runs the optimization; on_gen (optional) receives a dashboard each gen.
    // Returns the final population (front 0 = the Pareto set).
    std::vector<Candidate> run(
        const std::function<void(const EngineDashboard&)>& on_gen = nullptr);

private:
    const Evaluator& eval_;
    Params p_;
    std::vector<std::vector<double>> seed_genes_;
};

}  // namespace ga
}  // namespace aero
