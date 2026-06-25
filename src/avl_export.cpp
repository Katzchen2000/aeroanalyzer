#include "aeroanalyzer/avl_export.h"
#include "aeroanalyzer/geom.h"
#include <fstream>
#include <cmath>

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

    // K SECTION entries using canonical η breakpoints; bow formula for x_le/chord.
    for (int k = 0; k < K && k < 5; ++k) {
        double eta = geom::SECTION_ETA[k];
        double y   = eta * w.semi_span;
        double bow = 4.0 * eta * (1.0 - eta);
        double xle = y * std::tan(w.le_sweep) + w.le_bow * bow;
        double chord = w.root_chord + (w.tip_chord - w.root_chord) * eta + w.te_bow * bow;
        double twist = w.washout * eta * RAD2DEG;
        std::string dat = stem + "_s" + std::to_string(k) + "_" + sec_label[k] + ".dat";
        o << "SECTION\n";
        o << xle << " " << y << " 0.0 " << chord << " " << twist << "\n";
        o << "AFILE\n" << dat << "\n#\n";
    }
    return true;
}

}  // namespace avl
}  // namespace aero
