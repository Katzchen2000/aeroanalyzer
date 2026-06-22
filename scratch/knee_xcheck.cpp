// knee_xcheck.cpp - reconstruct the three AVL decks (knee / min_drag / min_mass)
// from the CURRENT out/*.avl + out/*.dat and run the panel solver on them, so
// their CLa / e / x_np / Cm / SM can be diffed against a headless AVL run on the
// SAME decks. Parses each .avl directly (planform: root/tip chord, semi-span,
// LE sweep, washout=tip Ainc; CG=Xref; mac=Cref) so it can never drift out of
// sync with out/ again -- the prior version hardcoded a stale planform + baked
// AVL numbers, which became a Frankenstein once out/ was regenerated.
//
// AVL columns are left blank for the user to fill from `avl352.exe` (OPER -> x):
//   reads CLa from two ST runs (or dCL/da), e and Xnp from the ST dump, Cm@a=2.
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/airfoil_io.h"
#include "aeroanalyzer/aero_potential.h"
#include "aeroanalyzer/massprops.h"
#include "aeroanalyzer/aero_viscous.h"
#include "aeroanalyzer/config.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace aero;

struct Deck {
    double Sref = 0, Cref = 0, Bref = 0, Xref = 0;
    double root_c = 0, tip_xle = 0, semi = 0, tip_c = 0, tip_ainc = 0;
    std::string dat;
    bool ok = false;
};

// Minimal AVL-deck parser: grabs the Sref/Cref/Bref line, the Xref (CG) line,
// and the two SECTION coordinate lines (root then tip).
static Deck parse_avl(const std::string& path) {
    Deck d;
    std::ifstream f(path);
    if (!f) return d;
    std::string line;
    int sect = 0;
    bool want_coords = false;
    while (std::getline(f, line)) {
        if (line.find("Sref") != std::string::npos) {
            std::istringstream s(line); s >> d.Sref >> d.Cref >> d.Bref; continue;
        }
        if (line.find("Xref") != std::string::npos) {
            std::istringstream s(line); s >> d.Xref; continue;
        }
        if (line.find("SECTION") != std::string::npos) { want_coords = true; continue; }
        if (want_coords) {
            // "Xle Yle Zle Chord Ainc"
            std::istringstream s(line);
            double xle, yle, zle, chord, ainc;
            if (s >> xle >> yle >> zle >> chord >> ainc) {
                if (sect == 0) { d.root_c = chord; }
                else { d.tip_xle = xle; d.semi = yle; d.tip_c = chord; d.tip_ainc = ainc; }
                ++sect; want_coords = false;
            }
            continue;
        }
        if (line.find("AFILE") != std::string::npos) {
            if (std::getline(f, line)) {
                // trim
                std::size_t a = line.find_first_not_of(" \t\r\n");
                std::size_t b = line.find_last_not_of(" \t\r\n");
                if (a != std::string::npos) d.dat = line.substr(a, b - a + 1);
            }
        }
    }
    d.ok = (d.root_c > 0 && d.tip_c > 0 && d.semi > 0 && !d.dat.empty());
    return d;
}

static void run(const std::string& name) {
    Deck d = parse_avl("out/" + name + ".avl");
    if (!d.ok) { printf("[%s] deck parse FAILED\n", name.c_str()); return; }

    bool ok = false;
    auto coords = airfoil_io::load_dat(d.dat, ok);
    if (!ok) { printf("[%s] failed to load %s\n", name.c_str(), d.dat.c_str()); return; }
    double te = airfoil_io::estimate_te(coords);
    Airfoil af = airfoil_io::to_airfoil(coords, 3, te);

    WingGeometry w;
    w.root_chord = d.root_c;
    w.tip_chord  = d.tip_c;
    w.semi_span  = d.semi;
    w.le_sweep   = std::atan2(d.tip_xle, d.semi);
    w.washout    = d.tip_ainc * DEG2RAD;   // deck Ainc == w.washout*RAD2DEG
    w.section    = af;
    geom::loft(w, 20);

    Config cfg; viscous::Surrogate s;
    s.load("data/surrogates/polar_coeffs.csv", cfg);
    cfg.set("panel_chordwise", "10"); cfg.set("panel_wake_chords", "20");
    MassProps mp = massprops::compute(w, cfg);

    // Use the DECK's mac and CG for an apples-to-apples SM vs AVL (we lack the
    // genome's battery_x, so massprops' CG would not match the deck's Xref).
    double mac = (d.Cref > 0) ? d.Cref : mp.mac;
    double cg  = d.Xref;

    auto solve = [&](const char* model, double a) {
        Config c = cfg; c.set("aero_model", model);
        return potential::solve(w, mp, s, c, a * DEG2RAD, 0.0);
    };
    printf("[%s] sweep=%.2f taper=%.3f wash=%.3f  deck: Sref=%.5f mac=%.5f b=%.5f CG=%.5f | recon mac=%.5f AR=%.3f\n",
           name.c_str(), w.le_sweep * RAD2DEG, w.tip_chord / w.root_chord,
           w.washout * RAD2DEG, d.Sref, d.Cref, d.Bref, d.Xref, mp.mac, mp.AR);
    for (const char* model : {"vlm", "panel"}) {
        AeroState s1 = solve(model, 1.0), s2 = solve(model, 5.0), s2d = solve(model, 2.0);
        double CLa = (s2.CL - s1.CL) / (4.0 * DEG2RAD);
        double sm_deck = (s2d.x_np - cg) / mac;
        printf("   %-5s CLa=%.4f  e=%.4f  x_np=%.5f  Cm@2=%+.4f  SM(deckCG)=%.4f\n",
               model, CLa, s2d.e, s2d.x_np, s2d.CM, sm_deck);
    }
    printf("   AVL   CLa=______  e=______  x_np=______  Cm@2=______  (fill from avl352.exe)\n");
}

int main() {
    printf("=== panel/VLM vs AVL cross-check on CURRENT out/ decks ===\n");
    run("knee");
    run("min_drag");
    run("min_mass");
    return 0;
}
