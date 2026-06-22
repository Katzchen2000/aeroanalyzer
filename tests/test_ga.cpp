#include "test_harness.h"
#include "aeroanalyzer/ga.h"

using namespace aero;

static Candidate mk(double d, double m, double s, double cv) {
    Candidate c;
    c.objectives = {d, m, s};
    c.cv = cv;
    return c;
}

// Deb's constraint-domination rules.
TEST(constraint_domination) {
    Candidate feas = mk(2, 2, 2, 0.0);
    Candidate infeas = mk(1, 1, 1, 0.5);     // "better" objectives but infeasible
    CHECK(ga::constrained_dominates(feas, infeas));
    CHECK(!ga::constrained_dominates(infeas, feas));

    Candidate low = mk(5, 5, 5, 0.1);
    Candidate high = mk(0, 0, 0, 0.9);
    CHECK(ga::constrained_dominates(low, high));   // both infeasible: lower cv

    Candidate a = mk(1, 1, 1, 0.0);
    Candidate b = mk(2, 2, 2, 0.0);
    CHECK(ga::constrained_dominates(a, b));        // both feasible: Pareto
    Candidate c = mk(1, 2, 1, 0.0);
    CHECK(ga::constrained_dominates(a, c));        // a <= c all, < in one
    CHECK(!ga::constrained_dominates(c, a));
}

// Fast non-dominated sort produces the known fronts.
// A=(1,1,0) dominates all; C=(1,2,0),D=(2,1,0) are front 1; B=(2,2,0) front 2.
TEST(nondominated_sort_fronts) {
    std::vector<Candidate> pop = {
        mk(1, 1, 0, 0), // A -> rank 0
        mk(2, 2, 0, 0), // B -> rank 2
        mk(1, 2, 0, 0), // C -> rank 1
        mk(2, 1, 0, 0), // D -> rank 1
    };
    std::vector<std::vector<int>> fronts;
    ga::fast_nondominated_sort(pop, fronts);
    CHECK(pop[0].rank == 0);
    CHECK(pop[2].rank == 1);
    CHECK(pop[3].rank == 1);
    CHECK(pop[1].rank == 2);
    CHECK(fronts.size() == 3);
    CHECK(fronts[0].size() == 1);
}

// Crowding distance: boundary points get infinity.
TEST(crowding_boundaries_infinite) {
    std::vector<Candidate> pop = {
        mk(0, 3, 0, 0), mk(1, 2, 0, 0), mk(2, 1, 0, 0), mk(3, 0, 0, 0),
    };
    std::vector<int> front = {0, 1, 2, 3};
    for (auto& c : pop) c.rank = 0;
    ga::crowding_distance(pop, front);
    int inf_count = 0;
    for (auto& c : pop)
        if (std::isinf(c.crowding)) ++inf_count;
    CHECK(inf_count >= 2);            // at least the two extremes per objective
}
