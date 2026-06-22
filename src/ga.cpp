#include "aeroanalyzer/ga.h"
#include <random>
#include <algorithm>
#include <limits>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace aero {
namespace ga {

bool constrained_dominates(const Candidate& a, const Candidate& b) {
    bool a_feas = (a.cv <= 0.0), b_feas = (b.cv <= 0.0);
    if (a_feas != b_feas) return a_feas;          // feasible beats infeasible
    if (!a_feas) return a.cv < b.cv;              // both infeasible: less viol.
    bool not_worse = true, strictly_better = false;
    for (int k = 0; k < N_OBJ; ++k) {
        if (a.objectives[k] > b.objectives[k]) not_worse = false;
        if (a.objectives[k] < b.objectives[k]) strictly_better = true;
    }
    return not_worse && strictly_better;
}

void fast_nondominated_sort(std::vector<Candidate>& pop,
                            std::vector<std::vector<int>>& fronts) {
    const int n = static_cast<int>(pop.size());
    fronts.clear();
    std::vector<std::vector<int>> S(n);
    std::vector<int> ndom(n, 0);
    std::vector<int> f0;
    for (int p = 0; p < n; ++p) {
        for (int qi = 0; qi < n; ++qi) {
            if (qi == p) continue;
            if (constrained_dominates(pop[p], pop[qi])) S[p].push_back(qi);
            else if (constrained_dominates(pop[qi], pop[p])) ndom[p]++;
        }
        if (ndom[p] == 0) { pop[p].rank = 0; f0.push_back(p); }
    }
    fronts.push_back(f0);
    int fi = 0;
    while (!fronts[fi].empty()) {
        std::vector<int> next;
        for (int p : fronts[fi]) {
            for (int qi : S[p]) {
                if (--ndom[qi] == 0) { pop[qi].rank = fi + 1; next.push_back(qi); }
            }
        }
        ++fi;
        fronts.push_back(next);
    }
    if (!fronts.empty() && fronts.back().empty()) fronts.pop_back();
}

void crowding_distance(std::vector<Candidate>& pop,
                       const std::vector<int>& front) {
    const int m = static_cast<int>(front.size());
    if (m == 0) return;
    for (int idx : front) pop[idx].crowding = 0.0;
    const double INF = std::numeric_limits<double>::infinity();
    for (int k = 0; k < N_OBJ; ++k) {
        std::vector<int> order(front);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return pop[a].objectives[k] < pop[b].objectives[k];
        });
        pop[order.front()].crowding = INF;
        pop[order.back()].crowding = INF;
        double fmin = pop[order.front()].objectives[k];
        double fmax = pop[order.back()].objectives[k];
        double range = fmax - fmin;
        if (range <= 0.0) continue;
        for (int i = 1; i < m - 1; ++i) {
            if (pop[order[i]].crowding == INF) continue;
            pop[order[i]].crowding +=
                (pop[order[i + 1]].objectives[k] -
                 pop[order[i - 1]].objectives[k]) / range;
        }
    }
}

namespace {
bool crowded_better(const Candidate& a, const Candidate& b) {
    if (a.rank != b.rank) return a.rank < b.rank;
    return a.crowding > b.crowding;
}

struct Rng {
    std::mt19937 gen;
    std::uniform_real_distribution<double> u{0.0, 1.0};
    explicit Rng(unsigned s) : gen(s) {}
    double operator()() { return u(gen); }
};

// SBX crossover for one variable, bounded.
void sbx(double& c1, double& c2, double p1, double p2, double lo, double hi,
         double eta, Rng& rng) {
    if (std::fabs(p1 - p2) < 1e-14) { c1 = p1; c2 = p2; return; }
    if (p1 > p2) std::swap(p1, p2);
    double rand = rng();
    auto spread = [&](double beta) {
        double a = 2.0 - std::pow(beta, -(eta + 1.0));
        return (rand <= 1.0 / a) ? std::pow(rand * a, 1.0 / (eta + 1.0))
                                 : std::pow(1.0 / (2.0 - rand * a),
                                            1.0 / (eta + 1.0));
    };
    double beta1 = 1.0 + 2.0 * (p1 - lo) / (p2 - p1);
    double beta2 = 1.0 + 2.0 * (hi - p2) / (p2 - p1);
    double bq1 = spread(beta1), bq2 = spread(beta2);
    c1 = 0.5 * ((p1 + p2) - bq1 * (p2 - p1));
    c2 = 0.5 * ((p1 + p2) + bq2 * (p2 - p1));
    c1 = std::max(lo, std::min(hi, c1));
    c2 = std::max(lo, std::min(hi, c2));
}

// Polynomial mutation for one variable, bounded.
void poly_mut(double& x, double lo, double hi, double eta, Rng& rng) {
    if (hi <= lo) return;
    double u = rng();
    double dx = hi - lo;
    double delta;
    if (u < 0.5) {
        double b = 2.0 * u + (1.0 - 2.0 * u) *
                   std::pow((hi - x) / dx, eta + 1.0);
        delta = std::pow(b, 1.0 / (eta + 1.0)) - 1.0;
    } else {
        double b = 2.0 * (1.0 - u) + 2.0 * (u - 0.5) *
                   std::pow((x - lo) / dx, eta + 1.0);
        delta = 1.0 - std::pow(b, 1.0 / (eta + 1.0));
    }
    x = std::max(lo, std::min(hi, x + delta * dx));
}
}  // namespace

NSGA2::NSGA2(const Evaluator& eval, const Params& p) : eval_(eval), p_(p) {}

std::vector<Candidate> NSGA2::run(
    const std::function<void(const EngineDashboard&)>& on_gen) {
    const geom::GenomeSpec& spec = eval_.spec();
    const int n_genes = static_cast<int>(spec.size());
    const double pmut = (p_.pmut < 0.0) ? 1.0 / n_genes : p_.pmut;
    Rng rng(p_.seed);

    auto eval_pop = [&](std::vector<Candidate>& pop) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < static_cast<int>(pop.size()); ++i)
            eval_.evaluate(pop[i]);
    };

    // ---- initial population: seeds first, then random explorers ----
    std::vector<Candidate> P(p_.pop);
    for (int i = 0; i < p_.pop; ++i) {
        Candidate& c = P[i];
        c.genes.resize(n_genes);
        bool seeded = (i < static_cast<int>(seed_genes_.size()));
        for (int g = 0; g < n_genes; ++g) {
            double v = (seeded && g < static_cast<int>(seed_genes_[i].size()))
                           ? seed_genes_[i][g]
                           : spec.lo[g] + rng() * (spec.hi[g] - spec.lo[g]);
            c.genes[g] = std::min(spec.hi[g], std::max(spec.lo[g], v));
        }
    }
    eval_pop(P);
    std::vector<std::vector<int>> fronts;
    fast_nondominated_sort(P, fronts);
    for (auto& f : fronts) crowding_distance(P, f);

    auto tournament = [&](const std::vector<Candidate>& pop) -> const Candidate& {
        int a = static_cast<int>(rng() * pop.size());
        int b = static_cast<int>(rng() * pop.size());
        if (a >= (int)pop.size()) a = (int)pop.size() - 1;
        if (b >= (int)pop.size()) b = (int)pop.size() - 1;
        return crowded_better(pop[a], pop[b]) ? pop[a] : pop[b];
    };

    for (int gen = 0; gen < p_.generations; ++gen) {
        // ---- offspring ----
        std::vector<Candidate> Q(p_.pop);
        for (int i = 0; i < p_.pop; i += 2) {
            const Candidate& p1 = tournament(P);
            const Candidate& p2 = tournament(P);
            Candidate c1, c2;
            c1.genes.resize(n_genes);
            c2.genes.resize(n_genes);
            for (int g = 0; g < n_genes; ++g) {
                double a = p1.genes[g], b = p2.genes[g];
                if (rng() < p_.pcx)
                    sbx(c1.genes[g], c2.genes[g], a, b, spec.lo[g], spec.hi[g],
                        p_.eta_cx, rng);
                else { c1.genes[g] = a; c2.genes[g] = b; }
                if (rng() < pmut)
                    poly_mut(c1.genes[g], spec.lo[g], spec.hi[g], p_.eta_mut, rng);
                if (rng() < pmut)
                    poly_mut(c2.genes[g], spec.lo[g], spec.hi[g], p_.eta_mut, rng);
            }
            Q[i] = c1;
            if (i + 1 < p_.pop) Q[i + 1] = c2;
        }
        eval_pop(Q);

        // ---- elitist merge + truncation ----
        std::vector<Candidate> R;
        R.reserve(P.size() + Q.size());
        R.insert(R.end(), P.begin(), P.end());
        R.insert(R.end(), Q.begin(), Q.end());
        fast_nondominated_sort(R, fronts);

        std::vector<Candidate> next;
        next.reserve(p_.pop);
        for (auto& f : fronts) {
            crowding_distance(R, f);
            if ((int)(next.size() + f.size()) <= p_.pop) {
                for (int idx : f) next.push_back(R[idx]);
            } else {
                std::vector<int> rem(f);
                std::sort(rem.begin(), rem.end(), [&](int a, int b) {
                    return R[a].crowding > R[b].crowding;
                });
                for (int idx : rem) {
                    if ((int)next.size() >= p_.pop) break;
                    next.push_back(R[idx]);
                }
                break;
            }
        }
        P.swap(next);
        fast_nondominated_sort(P, fronts);
        for (auto& f : fronts) crowding_distance(P, f);

        // ---- dashboard ----
        if (on_gen) {
            EngineDashboard d;
            d.generation = gen;
            d.pop_size = (int)P.size();
            d.front0_size = fronts.empty() ? 0 : (int)fronts[0].size();
            int feas = 0;
            int bd = -1, bm = -1, bs = -1;
            for (int i = 0; i < (int)P.size(); ++i) {
                if (P[i].cv <= 0.0) {
                    feas++;
                    if (bd < 0 || P[i].objectives[OBJ_DRAG] < P[bd].objectives[OBJ_DRAG]) bd = i;
                    if (bm < 0 || P[i].objectives[OBJ_MASS] < P[bm].objectives[OBJ_MASS]) bm = i;
                    if (bs < 0 || P[i].objectives[OBJ_SM] < P[bs].objectives[OBJ_SM]) bs = i;
                }
            }
            d.feasible_count = feas;
            if (bd >= 0) d.best_drag = P[bd];
            if (bm >= 0) d.best_mass = P[bm];
            if (bs >= 0) d.best_sm = P[bs];
            on_gen(d);
        }
    }
    return P;
}

}  // namespace ga
}  // namespace aero
