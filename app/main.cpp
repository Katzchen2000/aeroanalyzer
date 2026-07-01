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
void print_aero_report(const std::string& stem, const EvalResult& r,
                       const std::string& aero_model) {
    using std::cout; using std::fixed; using std::setprecision; using std::setw;
    constexpr double RAD2DEG_L = 180.0 / 3.14159265358979323846;
    auto w = [&](const char* label, double val, int dec, const char* unit="") {
        cout << "  " << std::left << setw(34) << label
             << std::right << setw(8) << fixed << setprecision(dec) << val
             << "  " << unit << "\n";
    };
    cout << "\n================================================================\n"
         << "  AERODYNAMIC REPORT \xe2\x80\x94  " << stem << "\n"
         << "================================================================\n";

    std::string solver_desc = (aero_model == "vlm")
        ? "analytic VLM" : "panel solver, relaxed wake";
    cout << "\n--- Trim condition (" << solver_desc << ") ---\n";
    w("Trim alpha",          r.aero.alpha * RAD2DEG_L, 3, "deg");
    w("Lift coefficient CL", r.aero.CL,    3);
    w("Induced drag CDi",    r.aero.CDi,   3);
    w("Oswald efficiency e", r.aero.e,     3);
    w("Mean aero chord",     r.mp.mac,     3, "m");
    w("CG position",         r.mp.x_cg,   3, "m");
    w("Neutral point",       r.aero.x_np,  3, "m");
    w("Static margin",       r.aero.static_margin * 100.0, 3, "% MAC");

    cout << "\n--- Performance (NSGA-II Pareto front) ---\n";
    w("Inverse L/D (CD/CL)", r.objectives[OBJ_DRAG], 3);
    w("Total mass",          r.objectives[OBJ_MASS], 3, "kg");
    w("Full span",           r.mp.b_full,  3, "m");
    w("Aspect ratio",        r.mp.AR,      3);
    w("Root chord",          r.geom.root_chord, 3, "m");
    w("Tip chord",           r.geom.tip_chord,  3, "m");
    w("LE sweep",            r.geom.le_sweep * RAD2DEG_L, 3, "deg");
    w("Washout",             r.geom.washout * RAD2DEG_L,  3, "deg");
    w("Total drag coeff CD", r.aero.CD,    3);
    w("Hinge moment",        r.aero.hinge_moment, 3, "kg*cm");
    cout << "  " << std::left << setw(34) << "Control mode"
         << (r.geom.mode == ControlMode::Elevon ? "elevon" : "split") << "\n";

    cout << "\n--- Lateral / directional dynamics ---\n";
    w("Roll authority pb/2V", r.aero.roll_helix, 3);

    cout << "\n--- Dynamic stability ---\n";
    w("Dutch-roll damp. zeta", r.aero.dutch_roll_zeta, 3);
    w("Phugoid damp. zeta",    r.aero.phugoid_zeta,    3);
    cout << "    Dutch-roll: ("
         << (r.aero.dutch_roll_zeta >= 0.05 ? "stable"
             : r.aero.dutch_roll_zeta >= 0.0 ? "lightly damped" : "divergent")
         << ")\n";
    cout << "    Phugoid:    ("
         << (r.aero.phugoid_zeta >= 0.05 ? "stable"
             : r.aero.phugoid_zeta >= 0.0 ? "lightly damped" : "divergent")
         << ")\n";

    cout << "\n--- Banked-turn stability (2g / 60deg bank, no extra solve) ---\n";
    w("Cn_beta at turn CL",    r.aero.cn_beta_turn,        3, "/rad");
    w("Dutch-roll zeta (turn)", r.aero.dutch_roll_zeta_turn, 3);
    cout << "  " << std::left << setw(34) << "Tip-stall in turn"
         << ": " << (r.aero.tip_stall_turn ? "YES" : "no") << "\n";
    cout << "    Cn_beta(turn): ("
         << (r.aero.cn_beta_turn > 0.0 ? "stable" : "UNSTABLE") << ")   "
         << "Dutch-roll: ("
         << (r.aero.dutch_roll_zeta_turn >= 0.0 ? "stable" : "divergent") << ")\n";
    cout << "\n";
}

void print_dash(const EngineDashboard& d) {
    if (d.generation % 10 != 0 && d.generation != 0) return;
    std::cout << "gen " << std::setw(4) << d.generation
              << " | feasible " << std::setw(4) << d.feasible_count
              << "/" << d.pop_size
              << " | front0 " << std::setw(4) << d.front0_size;
    if (d.feasible_count > 0) {
        std::cout << " | best L/D " << std::fixed << std::setprecision(3)
                  << (d.best_drag.objectives[OBJ_DRAG] > 0
                      ? 1.0 / d.best_drag.objectives[OBJ_DRAG] : 0.0) << " "
                  << " | min mass " << d.best_mass.objectives[OBJ_MASS] << " kg"
                  << " | best SMdev " << std::setprecision(4)
                  << d.best_sm.objectives[OBJ_SM];
    }
    std::cout << "\n" << std::flush;
}
}  // namespace

int main(int argc, char** argv) {
    std::string cfg_path = (argc > 1) ? argv[1] : "config/baseline.cfg";
    Config cfg;
    if (!cfg.load(cfg_path))
        std::cerr << "[warn] could not open " << cfg_path
                  << " — using built-in defaults\n";

    std::string aero_model = cfg.gets("aero_model", "panel");
    std::cout << "AeroAnalyzer Pro — NSGA-II (drag / mass / static-margin)\n";
    if (aero_model == "vlm")
        std::cout << "aerodynamics: analytic vortex-lattice FALLBACK "
                     "(aero_model = vlm).\n\n";
    else
        std::cout << "aerodynamics: Morino panel solver (default, AVL-validated; "
                     "set aero_model = vlm for the analytic fallback).\n\n";

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
    double cst_jitter = cfg.getd("cst_seed_jitter", 0.25);
    auto seed_genes = seeds::build_seed_genomes(eval.spec(), eval.seeds(),
                                                p.seed, n_seed, cst_jitter);
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
    csv << "idx,inv_LD,mass_kg,sm_dev,static_margin,span_m,AR,root_c,tip_c,"
           "sweep_deg,washout_deg,CL,CD,hinge_kgcm,roll_helix,mode,"
           "dutch_roll_zeta,phugoid_zeta,"
           "cn_beta_turn,dutch_roll_zeta_turn,tip_stall_turn";
    for (std::size_t g = 0; g < eval.spec().size(); ++g)
        csv << "," << eval.spec().names[g];
    csv << "\n";
    csv << std::fixed << std::setprecision(5);

    // Accumulate for summary (reuses the eval.detail calls made for the CSV)
    std::vector<double> s_drag, s_mass, s_sm, s_ar, s_span;

    for (std::size_t k = 0; k < front.size(); ++k) {
        const Candidate& c = pop[front[k]];
        EvalResult r = eval.detail(c.genes);
        csv << k << "," << c.objectives[OBJ_DRAG] << "," << c.objectives[OBJ_MASS]
            << "," << c.objectives[OBJ_SM] << "," << r.aero.static_margin << ","
            << r.mp.b_full << "," << r.mp.AR << "," << r.geom.root_chord << ","
            << r.geom.tip_chord << "," << (r.geom.le_sweep * RAD2DEG) << ","
            << (r.geom.washout * RAD2DEG) << "," << r.aero.CL << "," << r.aero.CD
            << "," << r.aero.hinge_moment << "," << r.aero.roll_helix << ","
            << (r.geom.mode == ControlMode::Elevon ? "elevon" : "split")
            << "," << r.aero.dutch_roll_zeta << "," << r.aero.phugoid_zeta
            << "," << r.aero.cn_beta_turn
            << "," << r.aero.dutch_roll_zeta_turn
            << "," << (r.aero.tip_stall_turn ? 1 : 0);
        for (std::size_t g = 0; g < c.genes.size(); ++g)
            csv << "," << c.genes[g];
        csv << "\n";
        s_drag.push_back(c.objectives[OBJ_DRAG]);
        s_mass.push_back(c.objectives[OBJ_MASS]);
        s_sm.push_back(r.aero.static_margin * 100.0);
        s_ar.push_back(r.mp.AR);
        s_span.push_back(r.mp.b_full);
    }
    csv.close();

    std::cout << "\nPareto set (feasible front 0): " << front.size()
              << " designs -> out/pareto.csv\n";

    if (!s_drag.empty()) {
        // Print min/median/max summary
        auto med = [](std::vector<double> v) -> double {
            std::sort(v.begin(), v.end());
            std::size_t n = v.size();
            return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
        };
        std::cout << std::fixed;
        std::cout << "\n  front summary         min       median       max\n";
        std::cout << "  drag   N :       " << std::setprecision(3)
                  << *std::min_element(s_drag.begin(),s_drag.end()) << "      "
                  << med(s_drag) << "      "
                  << *std::max_element(s_drag.begin(),s_drag.end()) << "\n";
        std::cout << "  mass   kg:       " << std::setprecision(3)
                  << *std::min_element(s_mass.begin(),s_mass.end()) << "      "
                  << med(s_mass) << "      "
                  << *std::max_element(s_mass.begin(),s_mass.end()) << "\n";
        std::cout << "  SM     % :       " << std::setprecision(2)
                  << *std::min_element(s_sm.begin(),s_sm.end()) << "       "
                  << med(s_sm) << "       "
                  << *std::max_element(s_sm.begin(),s_sm.end()) << "\n";
        std::cout << "  AR       :       " << std::setprecision(2)
                  << *std::min_element(s_ar.begin(),s_ar.end()) << "        "
                  << med(s_ar) << "        "
                  << *std::max_element(s_ar.begin(),s_ar.end()) << "\n";
        std::cout << "  span   m :       " << std::setprecision(3)
                  << *std::min_element(s_span.begin(),s_span.end()) << "      "
                  << med(s_span) << "      "
                  << *std::max_element(s_span.begin(),s_span.end()) << "\n";
    }

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
            // Re-loft at export resolution; GA uses n_stations for speed.
            WingGeometry rg = r.geom;
            geom::loft(rg, cfg.geti("n_stations_export", 60));
            avl::write_case(stem, rg, r.mp, cfg);
            avl::write_3d_csv(stem, rg, cfg);
            avl::write_stl(stem, rg, cfg);
            // Sidecar: panel reference numbers for validate_avl.ps1 cross-check.
            {
                std::ofstream sc(stem + "_panel.txt");
                sc << std::fixed << std::setprecision(6);
                sc << "alpha_deg = " << (r.aero.alpha * RAD2DEG) << "\n";
                sc << "CL        = " << r.aero.CL                << "\n";
                sc << "CDi       = " << r.aero.CDi               << "\n";
                sc << "e         = " << r.aero.e                 << "\n";
                sc << "x_np      = " << r.aero.x_np              << "\n";
                sc << "x_cg      = " << r.mp.x_cg                << "\n";
                sc << "mac       = " << r.mp.mac                 << "\n";
                sc << "sm        = " << r.aero.static_margin     << "\n";
                sc << "dutch_roll_zeta      = " << r.aero.dutch_roll_zeta       << "\n";
                sc << "phugoid_zeta         = " << r.aero.phugoid_zeta          << "\n";
                sc << "cn_beta              = " << r.aero.cn_beta               << "\n";
                sc << "dutch_roll_zeta_turn = " << r.aero.dutch_roll_zeta_turn  << "\n";
                sc << "cn_beta_turn         = " << r.aero.cn_beta_turn          << "\n";
                sc << std::setprecision(4);
                sc << "batt_box_center = " << r.mp.batt_cx << ", "
                   << r.mp.batt_cy << ", " << r.mp.batt_cz << "\n";
                sc << "batt_box_dims   = " << r.mp.batt_lx << " x "
                   << r.mp.batt_ly << " x " << r.mp.batt_lz << "  (L x W x H, m)\n";
            }
            print_aero_report(stem, r, aero_model);
            std::cout << "  " << std::setw(9) << pk.name
                      << " : L/D " << std::setprecision(3)
                      << (r.objectives[OBJ_DRAG] > 0 ? 1.0/r.objectives[OBJ_DRAG] : 0.0)
                      << ", mass "
                      << r.objectives[OBJ_MASS] << " kg, SM "
                      << std::setprecision(3) << (r.aero.static_margin * 100.0)
                      << " %, AR " << r.mp.AR
                      << ", DRz " << std::setprecision(3) << r.aero.dutch_roll_zeta
                      << "  -> " << stem << ".avl/.stl"
                      << "\n             CAD: mirror Y=0";
            {
                double bed_mm = 250.0, b_semi = r.geom.semi_span;
                if (r.mp.b_full * 1000.0 > bed_mm) {
                    int ns = (int)std::ceil(r.mp.b_full * 1000.0 / bed_mm) - 1;
                    for (int sp = 1; sp <= ns; ++sp)
                        std::cout << "  print-split Y="
                                  << std::setprecision(3) << (b_semi * sp / (ns+1)) << "m";
                }
            }
            std::cout << "\n             batt ctr (x,y,z)=("
                      << std::setprecision(4) << r.mp.batt_cx << ","
                      << r.mp.batt_cy << "," << r.mp.batt_cz
                      << ") box " << r.mp.batt_lx << "x"
                      << r.mp.batt_ly << "x" << r.mp.batt_lz << " m\n";
        }
    } else {
        std::cout << "No feasible designs found — relax constraints or check "
                     "config.\n";
    }
    return 0;
}
