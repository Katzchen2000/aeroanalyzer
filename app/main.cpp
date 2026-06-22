// main.cpp — AeroAnalyzer Pro driver. Thin: load config, run NSGA-II, report.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <cmath>

#include "aeroanalyzer/config.h"
#include "aeroanalyzer/evaluate.h"
#include "aeroanalyzer/ga.h"
#include "aeroanalyzer/avl_export.h"

using namespace aero;

namespace {
void print_dash(const EngineDashboard& d) {
    if (d.generation % 10 != 0 && d.generation != 0) return;
    std::cout << "gen " << std::setw(4) << d.generation
              << " | feasible " << std::setw(4) << d.feasible_count
              << "/" << d.pop_size
              << " | front0 " << std::setw(4) << d.front0_size;
    if (d.feasible_count > 0) {
        std::cout << " | min drag " << std::fixed << std::setprecision(3)
                  << d.best_drag.objectives[OBJ_DRAG] << " N"
                  << " | min mass " << d.best_mass.objectives[OBJ_MASS] << " kg"
                  << " | best SMdev " << std::setprecision(4)
                  << d.best_sm.objectives[OBJ_SM];
    }
    std::cout << "\n";
}
}  // namespace

int main(int argc, char** argv) {
    std::string cfg_path = (argc > 1) ? argv[1] : "config/baseline.cfg";
    Config cfg;
    if (!cfg.load(cfg_path))
        std::cerr << "[warn] could not open " << cfg_path
                  << " — using built-in defaults\n";

    std::cout << "AeroAnalyzer Pro — NSGA-II (drag / mass / static-margin)\n";
    std::cout << "NOTE: aerodynamics = analytic REFERENCE MODEL (Milestone 3 "
                 "swaps in the Morino panel solver).\n\n";

    Evaluator eval(cfg);

    ga::Params p;
    p.pop = cfg.geti("ga_pop", 120);
    p.generations = cfg.geti("ga_generations", 150);
    p.seed = static_cast<unsigned>(cfg.geti("ga_seed", 1));
    p.eta_cx = cfg.getd("ga_eta_cx", 15.0);
    p.eta_mut = cfg.getd("ga_eta_mut", 20.0);
    p.pcx = cfg.getd("ga_pcx", 0.9);
    p.pmut = cfg.getd("ga_pmut", -1.0);

    ga::NSGA2 opt(eval, p);

    // Ancestry injection: ~50% of the initial population from seed airfoils
    // (elites + hybrids), the rest random explorers (plan §3).
    int n_seed = static_cast<int>(std::lround(0.5 * p.pop));
    auto seed_genes = seeds::build_seed_genomes(eval.spec(), eval.seeds(),
                                                p.seed, n_seed);
    opt.set_seeds(seed_genes);
    std::cout << "seeded " << seed_genes.size() << " genomes from "
              << eval.seeds().airfoils.size() << " airfoil(s)";
    if (!eval.seeds().names.empty()) {
        std::cout << " [";
        for (std::size_t i = 0; i < eval.seeds().names.size(); ++i)
            std::cout << (i ? ", " : "") << eval.seeds().names[i];
        std::cout << "]";
    }
    std::cout << "\n\n";

    auto pop = opt.run(&print_dash);

    // ---- collect feasible Pareto set (front 0, cv == 0) ----
    std::vector<int> front;
    for (int i = 0; i < (int)pop.size(); ++i)
        if (pop[i].rank == 0 && pop[i].cv <= 0.0) front.push_back(i);

    std::filesystem::create_directories("out");
    std::ofstream csv("out/pareto.csv");
    csv << "idx,drag_N,mass_kg,sm_dev,static_margin,span_m,AR,root_c,tip_c,"
           "sweep_deg,washout_deg,CL,CD,hinge_kgcm,mode\n";
    csv << std::fixed << std::setprecision(5);
    for (std::size_t k = 0; k < front.size(); ++k) {
        const Candidate& c = pop[front[k]];
        EvalResult r = eval.detail(c.genes);
        csv << k << "," << c.objectives[OBJ_DRAG] << "," << c.objectives[OBJ_MASS]
            << "," << c.objectives[OBJ_SM] << "," << r.aero.static_margin << ","
            << r.mp.b_full << "," << r.mp.AR << "," << r.geom.root_chord << ","
            << r.geom.tip_chord << "," << (r.geom.le_sweep * RAD2DEG) << ","
            << (r.geom.washout * RAD2DEG) << "," << r.aero.CL << "," << r.aero.CD
            << "," << r.aero.hinge_moment << ","
            << (r.geom.mode == ControlMode::Elevon ? "elevon" : "split") << "\n";
    }
    csv.close();

    std::cout << "\nPareto set (feasible front 0): " << front.size()
              << " designs -> out/pareto.csv\n";

    if (!front.empty()) {
        auto pick = [&](int obj) {
            int best = front[0];
            for (int idx : front)
                if (pop[idx].objectives[obj] < pop[best].objectives[obj]) best = idx;
            return best;
        };
        int i_drag = pick(OBJ_DRAG), i_mass = pick(OBJ_MASS);
        // knee: min normalized L2 of the three objectives
        double dN = 1e-9, mN = 1e-9, sN = 1e-9;
        for (int idx : front) {
            dN = std::max(dN, pop[idx].objectives[OBJ_DRAG]);
            mN = std::max(mN, pop[idx].objectives[OBJ_MASS]);
            sN = std::max(sN, pop[idx].objectives[OBJ_SM]);
        }
        int i_knee = front[0];
        double best_l2 = 1e18;
        for (int idx : front) {
            double a = pop[idx].objectives[OBJ_DRAG] / dN;
            double b = pop[idx].objectives[OBJ_MASS] / mN;
            double cc = pop[idx].objectives[OBJ_SM] / sN;
            double l2 = a * a + b * b + cc * cc;
            if (l2 < best_l2) { best_l2 = l2; i_knee = idx; }
        }

        struct Pick { const char* name; int idx; };
        Pick picks[] = {{"min_drag", i_drag}, {"min_mass", i_mass},
                        {"knee", i_knee}};
        std::cout << "\nincumbents:\n";
        for (auto& pk : picks) {
            EvalResult r = eval.detail(pop[pk.idx].genes);
            std::string stem = std::string("out/") + pk.name;
            avl::write_case(stem, r.geom, r.mp, cfg);
            std::cout << "  " << std::setw(9) << pk.name
                      << " : drag " << std::setprecision(3)
                      << r.objectives[OBJ_DRAG] << " N, mass "
                      << r.objectives[OBJ_MASS] << " kg, SM "
                      << std::setprecision(3) << (r.aero.static_margin * 100.0)
                      << " %, AR " << r.mp.AR
                      << "  -> " << stem << ".avl\n";
        }
    } else {
        std::cout << "No feasible designs found — relax constraints or check "
                     "config.\n";
    }
    return 0;
}
