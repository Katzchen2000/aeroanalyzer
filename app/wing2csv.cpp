// wing2csv — AeroAnalyzer wing geometry -> 3D CSV + aerodynamic report
// Usage: wing2csv [stem]   e.g. wing2csv out/min_drag   (default: out/min_drag)
//
// Reads <stem>.avl + <stem>_s*_*.dat + <stem>_panel.txt + out/pareto.csv
// Writes <stem>_3d.csv  (Station_ID, Side, X, Y, Z — import into Fusion360 etc.)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ─── types ──────────────────────────────────────────────────────────────────

struct Pt2 { double x, y; };   // normalised airfoil coord

struct Section {
    std::string label;          // e.g. "s0_root"
    double le_x, le_y, le_z;
    double chord, twist_deg;
    std::string dat_path;
};

// ─── helpers ────────────────────────────────────────────────────────────────

static std::string trim_ws(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

// Simple key = value parser (our _panel.txt format)
static std::map<std::string, double> read_kv(const std::string& path) {
    std::map<std::string, double> m;
    std::ifstream f(path);
    if (!f) return m;
    std::string line;
    while (std::getline(f, line)) {
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = trim_ws(line.substr(0, eq));
        try { m[k] = std::stod(trim_ws(line.substr(eq + 1))); } catch (...) {}
    }
    return m;
}

// ─── AVL parser ─────────────────────────────────────────────────────────────
// Avoids seekg (unreliable on Windows text-mode CRLF files); uses a one-line
// push-back buffer instead.

static std::vector<Section> parse_avl(const std::string& path) {
    std::vector<Section> secs;
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; return secs; }

    std::string pending;
    auto nextline = [&](std::string& out) -> bool {
        if (!pending.empty()) { out = pending; pending.clear(); return true; }
        return (bool)std::getline(f, out);
    };

    std::string line;
    int idx = 0;
    while (nextline(line)) {
        if (trim_ws(line) != "SECTION") continue;
        std::string data;
        if (!nextline(data)) continue;
        Section s;
        std::istringstream ss(data);
        if (!(ss >> s.le_x >> s.le_y >> s.le_z >> s.chord >> s.twist_deg)) continue;
        // Look ahead for AFILE; push back any SECTION/SURFACE we hit
        while (nextline(line)) {
            std::string t = trim_ws(line);
            if (t == "SECTION" || t == "SURFACE") { pending = line; break; }
            if (t == "AFILE") {
                std::string p;
                if (nextline(p)) s.dat_path = trim_ws(p);
            }
        }
        s.label = "S" + std::to_string(idx++);
        secs.push_back(s);
    }
    return secs;
}

// ─── dat loader ─────────────────────────────────────────────────────────────
// Our write_dat: header line, then 161 (x y) pairs, upper TE->LE then lower LE->TE

static std::vector<Pt2> load_dat(const std::string& path) {
    std::vector<Pt2> pts;
    std::ifstream f(path);
    if (!f) { std::cerr << "[warn] cannot open dat: " << path << "\n"; return pts; }
    std::string line;
    std::getline(f, line);  // skip header
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        Pt2 p;
        if (ss >> p.x >> p.y) pts.push_back(p);
    }
    return pts;
}

// ─── 3D transform ────────────────────────────────────────────────────────────
// Airfoil: x_norm=0 at LE, x_norm=1 at TE; y_norm = thickness (up = positive)
// Wing:    X = chordwise (aft), Y = spanwise (right), Z = up

static void write_section_3d(std::ofstream& csv,
                              const std::string& label, const std::string& side,
                              const std::vector<Pt2>& pts,
                              double le_x, double le_y, double le_z,
                              double chord, double twist_deg, double y_sign) {
    double rad  = -twist_deg * M_PI / 180.0;
    double cosT = std::cos(rad);
    double sinT = std::sin(rad);

    for (const auto& p : pts) {
        double xc = p.x * chord;
        double zc = p.y * chord;
        double xr = xc * cosT - zc * sinT;
        double zr = xc * sinT + zc * cosT;
        double X  = le_x + xr;
        double Y  = y_sign * le_y;
        double Z  = le_z + zr;
        csv << label << "," << side << ","
            << std::fixed << std::setprecision(6)
            << X << "," << Y << "," << Z << "\n";
    }
}

// ─── pareto CSV matcher ──────────────────────────────────────────────────────
// Finds the row whose CL matches `cl_target` within tol.

static std::map<std::string, std::string> find_pareto_row(const std::string& path,
                                                           double cl_target) {
    std::ifstream f(path);
    if (!f) return {};
    std::string header_line;
    if (!std::getline(f, header_line)) return {};

    // Parse header
    std::vector<std::string> cols;
    {
        std::istringstream ss(header_line);
        std::string tok;
        while (std::getline(ss, tok, ',')) cols.push_back(trim_ws(tok));
    }
    int cl_col = -1;
    for (int i = 0; i < (int)cols.size(); ++i)
        if (cols[i] == "CL") { cl_col = i; break; }
    if (cl_col < 0) return {};

    // Find best-matching row
    std::string best_row;
    double best_err = 1e9;
    std::string row_line;
    while (std::getline(f, row_line)) {
        std::istringstream ss(row_line);
        std::string tok;
        std::vector<std::string> vals;
        while (std::getline(ss, tok, ',')) vals.push_back(tok);
        if ((int)vals.size() <= cl_col) continue;
        try {
            double err = std::fabs(std::stod(vals[cl_col]) - cl_target);
            if (err < best_err) { best_err = err; best_row = row_line; }
        } catch (...) {}
    }
    if (best_err > 0.01 || best_row.empty()) return {};

    // Build column map for winner; first occurrence wins (pareto.csv has
    // duplicate "mode" and "washout_deg" columns — metric col beats gene col)
    std::map<std::string, std::string> m;
    std::istringstream ss(best_row);
    std::string tok;
    std::vector<std::string> vals;
    while (std::getline(ss, tok, ',')) vals.push_back(tok);
    for (int i = 0; i < (int)std::min(cols.size(), vals.size()); ++i)
        if (!m.count(cols[i])) m[cols[i]] = vals[i];
    return m;
}

// ─── main ────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    std::string stem = (argc > 1) ? argv[1] : "out/min_drag";
    std::string avl_path    = stem + ".avl";
    std::string panel_path  = stem + "_panel.txt";
    std::string csv_out     = stem + "_3d.csv";
    std::string pareto_path = "out/pareto.csv";

    std::cout << "Wing2CSV  —  AeroAnalyzer geometry + aero report\n";
    std::cout << "Stem: " << stem << "\n\n";

    // 1. Parse AVL ─────────────────────────────────────────────────────────
    auto sections = parse_avl(avl_path);
    if (sections.empty()) {
        std::cerr << "No SECTION entries found in " << avl_path << "\n";
        return 1;
    }
    std::cout << "Parsed " << sections.size() << " wing sections from " << avl_path << "\n";

    // 2. Load dat files + write 3D CSV ─────────────────────────────────────
    std::ofstream csv(csv_out);
    if (!csv) { std::cerr << "Cannot write " << csv_out << "\n"; return 1; }
    csv << "Station_ID,Side,X,Y,Z\n";

    for (const auto& sec : sections) {
        auto pts = load_dat(sec.dat_path);
        if (pts.empty()) {
            std::cerr << "[warn] no points for " << sec.label << "\n";
            continue;
        }
        // Extract span-label from dat path (e.g. "s0_root") for readability
        std::string sp_label = sec.label;
        {
            size_t last = sec.dat_path.find_last_of("/\\");
            std::string fname = (last == std::string::npos) ? sec.dat_path
                                                            : sec.dat_path.substr(last+1);
            // "min_drag_s0_root.dat" -> find last underscore group before .dat
            size_t dot = fname.rfind(".dat");
            if (dot != std::string::npos) {
                std::string noext = fname.substr(0, dot);
                // find second-to-last '_' to get "s0_root" suffix
                size_t p1 = noext.rfind('_');
                if (p1 != std::string::npos) {
                    size_t p2 = noext.rfind('_', p1-1);
                    if (p2 != std::string::npos)
                        sp_label = noext.substr(p2+1);
                }
            }
        }
        write_section_3d(csv, sp_label, "R",
                         pts, sec.le_x, sec.le_y, sec.le_z,
                         sec.chord, sec.twist_deg, +1.0);
        write_section_3d(csv, sp_label, "L",
                         pts, sec.le_x, sec.le_y, sec.le_z,
                         sec.chord, sec.twist_deg, -1.0);
        std::cout << "  " << sp_label
                  << "  y=" << std::fixed << std::setprecision(3) << sec.le_y
                  << " m  chord=" << sec.chord
                  << " m  twist=" << sec.twist_deg << " deg"
                  << "  (" << pts.size() << " pts)\n";
    }
    csv.close();
    std::cout << "\nCSV written: " << csv_out << "\n";
    std::cout << "(import as XY point cloud in Fusion360 / SolidWorks — one loft per Station_ID)\n";

    // 3. Panel metrics ─────────────────────────────────────────────────────
    auto pv = read_kv(panel_path);

    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "  AERODYNAMIC REPORT  —  " << stem << "\n";
    std::cout << "================================================================\n";

    auto pf = [](const char* label, double v, const char* unit) {
        std::cout << "  " << std::left << std::setw(28) << label
                  << std::right << std::fixed << std::setw(10)
                  << v << "  " << unit << "\n";
    };
    auto ps = [](const char* label, const std::string& v) {
        std::cout << "  " << std::left << std::setw(28) << label
                  << std::right << std::setw(10) << v << "\n";
    };

    if (!pv.empty()) {
        std::cout << "\n--- Trim condition (panel solver, relaxed wake) ---\n";
        if (pv.count("alpha_deg")) pf("Trim alpha",      pv["alpha_deg"], "deg");
        if (pv.count("CL"))        pf("Lift coefficient CL",  pv["CL"],        "");
        if (pv.count("CDi"))       pf("Induced drag CDi",     pv["CDi"],       "");
        if (pv.count("e"))         pf("Oswald efficiency e",  pv["e"],         "");
        if (pv.count("mac"))       pf("Mean aero chord",  pv["mac"],           "m");
        if (pv.count("x_cg"))      pf("CG position",     pv["x_cg"],          "m");
        if (pv.count("x_np"))      pf("Neutral point",   pv["x_np"],          "m");
        if (pv.count("sm"))        pf("Static margin",   pv["sm"] * 100.0,    "% MAC");
    } else {
        std::cout << "[warn] " << panel_path << " not found or empty\n";
    }

    // 4. Full pareto metrics ───────────────────────────────────────────────
    double cl_key = pv.count("CL") ? pv["CL"] : -1.0;
    auto row = find_pareto_row(pareto_path, cl_key);

    auto get = [&](const std::string& k, double fallback = 0.0) -> double {
        if (!row.count(k)) return fallback;
        try { return std::stod(row.at(k)); } catch (...) { return fallback; }
    };

    if (!row.empty()) {
        std::cout << "\n--- Performance (NSGA-II Pareto front) ---\n";
        pf("Drag force",           get("drag_N"),       "N");
        pf("Total mass",           get("mass_kg"),      "kg");
        pf("Full span",            get("span_m"),       "m");
        pf("Aspect ratio",         get("AR"),           "");
        pf("Root chord",           get("root_c"),       "m");
        pf("Tip chord",            get("tip_c"),        "m");
        pf("LE sweep",             get("sweep_deg"),    "deg");
        pf("Washout",              get("washout_deg"),  "deg");
        pf("Total drag coeff CD",  get("CD"),           "");
        pf("Hinge moment",         get("hinge_kgcm"),   "kg*cm");
        if (row.count("mode")) ps("Control mode", row.at("mode"));

        std::cout << "\n--- Lateral / directional dynamics ---\n";
        pf("Roll authority pb/2V", get("roll_helix"), "");
        std::cout << "\n--- Dynamic stability ---\n";
        double dr_z = get("dutch_roll_zeta");
        double ph_z = get("phugoid_zeta");
        pf("Dutch-roll damp. zeta",  dr_z,  "");
        pf("Phugoid damp. zeta",     ph_z,  "");
        // Qualitative assessment
        auto flag = [](double z, double floor) -> const char* {
            if (z >= floor) return "(stable)";
            if (z >= 0.0)   return "(lightly damped)";
            return "(UNSTABLE)";
        };
        std::cout << "    Dutch-roll: " << flag(dr_z, 0.08) << "\n";
        std::cout << "    Phugoid:    " << flag(ph_z, 0.04) << "\n";
    } else {
        std::cout << "\n[info] pareto.csv not found or no CL match — "
                     "run the optimizer first.\n";
    }

    // 5. Wing geometry summary ─────────────────────────────────────────────
    std::cout << "\n--- Wing planform (from AVL sections) ---\n";
    if (sections.size() >= 2) {
        double b_semi = sections.back().le_y;
        std::cout << "  Semi-span              : " << std::fixed << std::setprecision(3)
                  << b_semi << " m   (full span " << 2.0*b_semi << " m)\n";
        std::cout << "  Control sections       : " << sections.size() << "\n";
        std::cout << "  Span stations:\n";
        for (const auto& sec : sections) {
            double eta = (b_semi > 0) ? sec.le_y / b_semi : 0.0;
            std::cout << "    eta=" << std::setprecision(3) << eta
                      << "  y=" << sec.le_y << " m"
                      << "  chord=" << sec.chord << " m"
                      << "  twist=" << sec.twist_deg << " deg\n";
        }
    }

    std::cout << "================================================================\n";
    return 0;
}
