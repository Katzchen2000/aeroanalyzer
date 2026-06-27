#include "aeroanalyzer/avl_export.h"
#include "aeroanalyzer/geom.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>

namespace aero {
namespace avl {

static bool write_dat(const std::string& path, const Airfoil& f) {
    std::ofstream o(path);
    if (!o) return false;
    o << "AeroAnalyzer section\n";
    const int N = 80;
    // upper surface TE -> LE
    for (int i = N; i >= 0; --i) {
        double th = PI * i / N;
        double x = 0.5 * (1.0 - std::cos(th));
        o << x << " " << geom::cst_upper(f, x) << "\n";
    }
    // lower surface LE -> TE
    for (int i = 1; i <= N; ++i) {
        double th = PI * i / N;
        double x = 0.5 * (1.0 - std::cos(th));
        o << x << " " << geom::cst_lower(f, x) << "\n";
    }
    return true;
}

bool write_case(const std::string& stem, const WingGeometry& w,
                const MassProps& mp, const Config& cfg) {
    const int K = (int)w.sections.size();
    if (K == 0) return false;

    // Write per-section .dat files
    const char* sec_label[] = {"root","eta50","eta75","eta875","tip"};
    for (int k = 0; k < K && k < 5; ++k) {
        std::string dat = stem + "_s" + std::to_string(k) + "_" + sec_label[k] + ".dat";
        if (!write_dat(dat, w.sections[k])) return false;
    }

    std::ofstream o(stem + ".avl");
    if (!o) return false;
    double mach = cfg.getd("v_cruise", V_CRUISE) / 340.0;
    o << "AeroAnalyzer flying wing\n";
    o << mach << "                      | Mach\n";
    o << "0     0     0.0          | iYsym  iZsym  Zsym\n";
    o << mp.S_ref << " " << mp.mac << " " << mp.b_full
      << "   | Sref Cref Bref\n";
    o << mp.x_cg << " 0.0 0.0      | Xref Yref Zref (CG)\n";
    o << "#\n";
    o << "SURFACE\n";
    o << "Wing\n";
    o << "12  1.0   24  -1.5       | Nchord Cspace Nspan Sspace\n";
    o << "YDUPLICATE\n0.0\n";
    o << "SCALE\n1.0 1.0 1.0\n";
    o << "TRANSLATE\n0.0 0.0 0.0\n";
    o << "ANGLE\n0.0\n#\n";

    // K SECTION entries using canonical η breakpoints; power-law + gull + winglet fold.
    const double b_h     = w.semi_span;
    const double eta_wl  = w.winglet_eta;
    const double cant    = w.winglet_cant;
    for (int k = 0; k < K && k < 5; ++k) {
        double eta = geom::SECTION_ETA[k];
        double se  = std::pow(eta, w.sweep_exp);
        double te  = std::pow(eta, w.chord_exp);
        double xle   = b_h * std::tan(w.le_sweep) * se;
        double chord = w.root_chord - (w.root_chord - w.tip_chord) * te;
        double twist = w.washout * eta * RAD2DEG;
        double phys_y = 0.0;
        double phys_z = 0.0;
        if (w.stations.empty()) {
            phys_y = eta * b_h;
            phys_z = 0.0;
        } else if (eta <= w.stations.front().eta) {
            phys_y = w.stations.front().y;
            phys_z = w.stations.front().z;
        } else if (eta >= w.stations.back().eta) {
            phys_y = w.stations.back().y;
            phys_z = w.stations.back().z;
        } else {
            int idx = 0;
            while (idx < (int)w.stations.size() && w.stations[idx].eta < eta) {
                idx++;
            }
            const auto& s0 = w.stations[idx - 1];
            const auto& s1 = w.stations[idx];
            double denom = s1.eta - s0.eta;
            double f = (denom > 1e-9) ? (eta - s0.eta) / denom : 0.0;
            phys_y = s0.y + f * (s1.y - s0.y);
            phys_z = s0.z + f * (s1.z - s0.z);
        }
        std::string dat = stem + "_s" + std::to_string(k) + "_" + sec_label[k] + ".dat";
        o << "SECTION\n";
        o << xle << " " << phys_y << " " << phys_z << " " << chord << " " << twist << "\n";
        o << "AFILE\n" << dat << "\n#\n";
    }
    return true;
}

bool write_3d_csv(const std::string& stem, const WingGeometry& w) {
    if (w.stations.empty()) return false;
    std::ofstream csv(stem + "_3d.csv");
    if (!csv) return false;
    csv << "Station_ID,Side,X,Y,Z\n";
    csv << std::fixed << std::setprecision(6);
    const int N = 80;
    for (int si = 0; si < (int)w.stations.size(); ++si) {
        const Station& s = w.stations[si];
        double rad  = -s.twist;
        double cosT = std::cos(rad);
        double sinT = std::sin(rad);
        std::string label = "st" + std::to_string(si);
        for (int side = 0; side < 2; ++side) {
            double y_sign      = (side == 0) ? +1.0 : -1.0;
            const char* sname  = (side == 0) ? "R"   : "L";
            // upper surface TE→LE
            for (int i = N; i >= 0; --i) {
                double th = PI * i / N;
                double x  = 0.5 * (1.0 - std::cos(th));
                double xc = x * s.chord;
                double zc = geom::cst_upper(s.af, x) * s.chord;
                double zr = xc * sinT + zc * cosT;
                double Y  = y_sign * (s.y - zr * std::sin(s.dihedral));
                double Z  = s.z + zr * std::cos(s.dihedral);
                csv << label << "," << sname << ","
                    << (s.x_le + xc * cosT - zc * sinT) << ","
                    << Y << ","
                    << Z << "\n";
            }
            // lower surface LE→TE
            for (int i = 1; i <= N; ++i) {
                double th = PI * i / N;
                double x  = 0.5 * (1.0 - std::cos(th));
                double xc = x * s.chord;
                double zc = geom::cst_lower(s.af, x) * s.chord;
                double zr = xc * sinT + zc * cosT;
                double Y  = y_sign * (s.y - zr * std::sin(s.dihedral));
                double Z  = s.z + zr * std::cos(s.dihedral);
                csv << label << "," << sname << ","
                    << (s.x_le + xc * cosT - zc * sinT) << ","
                    << Y << ","
                    << Z << "\n";
            }
        }
    }
    return true;
}

bool write_stl(const std::string& stem, const WingGeometry& w, const Config& cfg) {
    if (w.stations.empty()) return false;
    std::ofstream stl(stem + ".stl");
    if (!stl) return false;
    stl << std::fixed << std::setprecision(6);
    stl << "solid wing\n";

    const int N = 80;                // chordwise half-count (matches write_3d_csv)
    const int M = 2 * N;             // closed-loop vertex count (TE not repeated)

    using Pt3 = std::array<double, 3>;

    // Build closed contour for one station / side. Returns 2N+1 pts;
    // [0]=[2N]=TE (sharp TE → same point), so use first M indices as closed loop.
    auto make_contour = [&](const Station& s, double ys) {
        std::vector<Pt3> pts;
        pts.reserve(2 * N + 1);
        double rad = -s.twist, cosT = std::cos(rad), sinT = std::sin(rad);
        auto push = [&](double x, double zc_hat) {
            double xc = x * s.chord, zc = zc_hat * s.chord;
            double zr = xc * sinT + zc * cosT;
            pts.push_back({s.x_le + xc * cosT - zc * sinT,
                           ys * (s.y - zr * std::sin(s.dihedral)),
                           s.z  + zr * std::cos(s.dihedral)});
        };
        for (int i = N; i >= 0; --i) {   // upper: TE → LE
            double th = PI * i / N, x = 0.5 * (1.0 - std::cos(th));
            push(x, geom::cst_upper(s.af, x));
        }
        for (int i = 1; i <= N; ++i) {   // lower: LE → TE
            double th = PI * i / N, x = 0.5 * (1.0 - std::cos(th));
            push(x, geom::cst_lower(s.af, x));
        }
        return pts;
    };

    // Emit one STL facet; normal via cross product (always geometrically correct).
    auto tri = [&](const Pt3& A, const Pt3& B, const Pt3& C) {
        double ux=B[0]-A[0], uy=B[1]-A[1], uz=B[2]-A[2];
        double vx=C[0]-A[0], vy=C[1]-A[1], vz=C[2]-A[2];
        double nx=uy*vz-uz*vy, ny=uz*vx-ux*vz, nz=ux*vy-uy*vx;
        double len = std::sqrt(nx*nx+ny*ny+nz*nz);
        if (len > 1e-12) { nx/=len; ny/=len; nz/=len; }
        stl << "facet normal " << nx << " " << ny << " " << nz << "\n"
               "  outer loop\n"
               "    vertex " << A[0] << " " << A[1] << " " << A[2] << "\n"
               "    vertex " << B[0] << " " << B[1] << " " << B[2] << "\n"
               "    vertex " << C[0] << " " << C[1] << " " << C[2] << "\n"
               "  endloop\nendfacet\n";
    };

    const int n = (int)w.stations.size();

    for (int side = 0; side < 2; ++side) {
        double ys = (side == 0) ? +1.0 : -1.0;

        // Build all station contours for this side
        std::vector<std::vector<Pt3>> C(n);
        for (int si = 0; si < n; ++si) C[si] = make_contour(w.stations[si], ys);

        // Lateral surface: quad-strip between adjacent stations → 2 triangles each
        for (int si = 0; si < n - 1; ++si) {
            for (int j = 0; j < M; ++j) {
                int j1 = (j + 1) % M;
                tri(C[si][j],   C[si+1][j],  C[si+1][j1]);
                tri(C[si][j],   C[si+1][j1], C[si][j1]);
            }
        }
        // Root cap: fan from LE point (contour index N = end of upper surface)
        for (int j = 0; j < M; ++j)
            tri(C[0][N], C[0][(j+1)%M], C[0][j]);
        // Tip cap: fan from LE (reversed winding vs root so both normals face outward)
        for (int j = 0; j < M; ++j)
            tri(C[n-1][N], C[n-1][j], C[n-1][(j+1)%M]);
    }

    // Prop-disk visual marker: triangle fan in the Y-Z plane at x = root_chord + hub_gap.
    // ponytail: visual only, not structurally part of the wing mesh.
    double prop_r  = 0.5 * cfg.getd("prop_diameter", 0.203);
    double hub_gap = cfg.getd("prop_hub_gap", 0.010);
    double x_disk  = w.root_chord + hub_gap;
    const int ND   = 24;
    Pt3 ctr = {x_disk, 0.0, 0.0};
    for (int i = 0; i < ND; ++i) {
        double a0 = 2.0*PI*i/ND, a1 = 2.0*PI*(i+1)/ND;
        Pt3 p0 = {x_disk, prop_r*std::cos(a0), prop_r*std::sin(a0)};
        Pt3 p1 = {x_disk, prop_r*std::cos(a1), prop_r*std::sin(a1)};
        tri(ctr, p0, p1);
    }

    stl << "endsolid wing\n";
    return true;
}

}  // namespace avl
}  // namespace aero
