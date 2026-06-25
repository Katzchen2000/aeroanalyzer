// aero_panel.cpp — Morino (Dirichlet) constant-strength source/doublet panel
// method. STAGE 1: influence kernels (this file grows over the later stages;
// solve() is a stub until the mesh/assembly/forces stages land).
//
// Kernel conventions (calibrated against brute-force surface quadrature in
// scratch/panel_probe.cpp):
//   doublet_potential(q,P) = -Omega(P) / (4*pi)      Omega = signed solid angle
//   source_potential(q,P)  = -(1/4pi) * INT_q (1/|P-Q|) dS
// with Omega computed by the Van Oosterom-Strackee formula (robust, no slope
// singularities) and the source surface integral by the edge-log closed form
// minus the off-plane solid-angle correction.

#include "aeroanalyzer/aero_panel.h"
#include "aeroanalyzer/control.h"
#include "aeroanalyzer/geom.h"
#include "aeroanalyzer/aero_potential.h"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <limits>

namespace aero {
namespace panel {

namespace {
constexpr double FOUR_PI = 4.0 * PI;
constexpr double KEPS = 1.0e-12;
}  // namespace

Vec3 quad_centroid(const Quad& q) {
    return (q.c[0] + q.c[1] + q.c[2] + q.c[3]) * 0.25;
}

Vec3 quad_normal(const Quad& q) {
    Vec3 d1 = q.c[2] - q.c[0];
    Vec3 d2 = q.c[3] - q.c[1];
    Vec3 n = d1.cross(d2);
    double L = n.norm();
    return (L > KEPS) ? n * (1.0 / L) : Vec3(0, 0, 1);
}

double quad_area(const Quad& q) {
    // Two triangles (0,1,2) and (0,2,3).
    Vec3 a = (q.c[1] - q.c[0]).cross(q.c[2] - q.c[0]);
    Vec3 b = (q.c[2] - q.c[0]).cross(q.c[3] - q.c[0]);
    return 0.5 * (a.norm() + b.norm());
}

// Signed solid angle subtended by the quad at P (fan of triangles from c0).
// Sign: positive when P sees the corners wound CCW about the panel's +normal
// (i.e. P on the +normal side); -> -2*pi for a point on the -normal side just
// behind the panel. Van Oosterom & Strackee (1983).
static double solid_angle(const Quad& q, const Vec3& P) {
    Vec3 r[4];
    double n[4];
    for (int i = 0; i < 4; ++i) { r[i] = q.c[i] - P; n[i] = r[i].norm(); }
    double omega = 0.0;
    // triangles (0,1,2) and (0,2,3)
    const int tri[2][3] = {{0, 1, 2}, {0, 2, 3}};
    for (auto& t : tri) {
        const Vec3& a = r[t[0]]; const Vec3& b = r[t[1]]; const Vec3& c = r[t[2]];
        double na = n[t[0]], nb = n[t[1]], nc = n[t[2]];
        double num = a.dot(b.cross(c));
        double den = na * nb * nc + a.dot(b) * nc + b.dot(c) * na + c.dot(a) * nb;
        omega += 2.0 * std::atan2(num, den);
    }
    return omega;
}

double doublet_potential(const Quad& q, const Vec3& P) {
    return -solid_angle(q, P) / FOUR_PI;
}

double source_potential(const Quad& q, const Vec3& P) {
    // Work in the panel local frame: origin at centroid, +z along the normal.
    Vec3 o = quad_centroid(q);
    Vec3 nz = quad_normal(q);
    Vec3 t = q.c[1] - q.c[0];
    Vec3 lx = t - nz * t.dot(nz);
    double lxn = lx.norm();
    lx = (lxn > KEPS) ? lx * (1.0 / lxn) : Vec3(1, 0, 0);
    Vec3 ly = nz.cross(lx);

    auto to_local = [&](const Vec3& v) {
        Vec3 d = v - o;
        return Vec3(d.dot(lx), d.dot(ly), d.dot(nz));
    };
    Vec3 Pl = to_local(P);
    double x = Pl.x, y = Pl.y, z = Pl.z;
    double X[4], Y[4];
    for (int i = 0; i < 4; ++i) { Vec3 cl = to_local(q.c[i]); X[i] = cl.x; Y[i] = cl.y; }

    double edge = 0.0;
    for (int k = 0; k < 4; ++k) {
        int k1 = (k + 1) & 3;
        double dx = X[k1] - X[k], dy = Y[k1] - Y[k];
        double d = std::sqrt(dx * dx + dy * dy);
        if (d < KEPS) continue;
        double rk  = std::sqrt((x - X[k]) * (x - X[k]) + (y - Y[k]) * (y - Y[k]) + z * z);
        double rk1 = std::sqrt((x - X[k1]) * (x - X[k1]) + (y - Y[k1]) * (y - Y[k1]) + z * z);
        double num = (x - X[k]) * dy - (y - Y[k]) * dx;     // perpendicular moment / d
        double denom = rk + rk1 - d;
        if (denom < KEPS) denom = KEPS;
        edge += (num / d) * std::log((rk + rk1 + d) / denom);
    }
    // INT 1/r dS = -(edge-log term) - |z| * (solid angle magnitude).
    // (sign of the edge sum calibrated against brute-force quadrature.)
    double omega = solid_angle(q, P);
    double integral = -edge - std::fabs(z) * std::fabs(omega);
    return -integral / FOUR_PI;
}

// ---- Surface mesh of the closed body -------------------------------------
namespace {

struct Panel {
    Quad q;
    Vec3 cp;       // collocation point (centroid nudged just inside the body)
    Vec3 n;        // outward unit normal
    double area = 0.0;
    int strip = 0;    // spanwise strip index
    int surf = 0;     // 0 = lower, 1 = upper
    int jc = 0;       // chordwise index (0 = LE-most)
};

struct Mesh {
    std::vector<Panel> panels;
    int n_strips = 0;
    int nc = 0;                 // chordwise panels per surface
    std::vector<int> te_up;     // per strip: index of upper-TE panel
    std::vector<int> te_lo;     // per strip: index of lower-TE panel
    std::vector<double> strip_y, strip_dy, strip_chord;   // per strip
    std::vector<double> strip_xqc;                        // per strip: quarter-chord x
    std::vector<Vec3>  te_pt;   // per full-station: camber TE point (size strips+1)
};

// Section point: chord-local (a from LE, b = thickness ordinate z/c), rotated
// by the station's geometric twist about the LE, placed into 3D.
Vec3 section_point(const Station& s, double a, double b) {
    // Geometric twist must act like a LOCAL angle-of-attack increment, matching
    // both the VLM convention (rhs = -(alpha + twist - aL0)) and the freestream
    // Vinf = (cos a, 0, sin a): a positive twist tilts the section nose-up, i.e.
    // the TE moves DOWN (-z), which presents positive incidence to the +x flow.
    // (A previous "fix" used +twist and inverted the washout loading sign.)
    double tw = s.twist;
    double ct = std::cos(tw), stt = std::sin(tw);
    double xr =  a * ct + b * stt;
    double zr = -a * stt + b * ct;
    return Vec3(s.x_le + s.chord * xr, s.y, s.z + s.chord * zr);
}

// Build the HALF-wing (y in [0,semi]) closed-body mesh; the other half is
// represented by a y=0 mirror image folded into the influence coefficients
// (see build_system). This forces a symmetric solution — excluding the spurious
// antisymmetric circulation mode that twist would otherwise excite — and halves
// the matrix dimension. nc chordwise panels per surface.
Mesh build_mesh(const WingGeometry& w, int nc, bool half_cosine) {
    Mesh m;
    m.nc = nc;

    const std::vector<Station>& full = w.stations;   // half-wing, root..tip

    // Chordwise nodes. Full cosine clusters at BOTH the LE and the TE; that TE
    // clustering makes the trailing-edge panel vanishingly thin at high nc,
    // driving the upper/lower TE collocation points toward coincidence and the
    // TE doublet-jump (from which CL is extracted) into catastrophic
    // cancellation -- the chordwise non-convergence (CLa diverges and sign-flips
    // as nc grows). Half cosine (the default) has zero slope at the LE (still
    // clustered for the suction peak) and maximum slope at the TE (kept coarse),
    // so the TE panels stay well separated and CLa converges under refinement.
    std::vector<double> xi(nc + 1);
    for (int j = 0; j <= nc; ++j)
        xi[j] = half_cosine ? (1.0 - std::cos(0.5 * PI * j / nc))
                            : 0.5 * (1.0 - std::cos(PI * j / nc));

    m.n_strips = static_cast<int>(full.size()) - 1;
    m.te_up.assign(m.n_strips, -1);
    m.te_lo.assign(m.n_strips, -1);

    // Per-strip span/chord and per-station camber TE point (for the wake).
    m.strip_y.resize(m.n_strips);
    m.strip_dy.resize(m.n_strips);
    m.strip_chord.resize(m.n_strips);
    m.strip_xqc.resize(m.n_strips);
    for (int sp = 0; sp < m.n_strips; ++sp) {
        m.strip_y[sp]     = 0.5 * (full[sp].y + full[sp + 1].y);
        m.strip_dy[sp]    = std::fabs(full[sp + 1].y - full[sp].y);
        m.strip_chord[sp] = 0.5 * (full[sp].chord + full[sp + 1].chord);
        double xle = 0.5 * (full[sp].x_le + full[sp + 1].x_le);
        m.strip_xqc[sp]   = xle + 0.25 * m.strip_chord[sp];   // quarter-chord x
    }
    m.te_pt.resize(full.size());
    for (std::size_t i = 0; i < full.size(); ++i)
        m.te_pt[i] = section_point(full[i], 1.0, 0.0);   // camber TE (ord=0)

    auto add_surface = [&](int surf) {   // surf: 1 upper, 0 lower
        for (int sp = 0; sp < m.n_strips; ++sp) {
            const Station& A = full[sp];
            const Station& B = full[sp + 1];
            for (int j = 0; j < nc; ++j) {
                double x0 = xi[j], x1 = xi[j + 1];
                auto ord = [&](const Station& S, double x) {
                    return surf ? geom::cst_upper(S.af, x)
                                : geom::cst_lower(S.af, x);
                };
                // Four corners of the panel.
                Vec3 pA0 = section_point(A, x0, ord(A, x0));
                Vec3 pA1 = section_point(A, x1, ord(A, x1));
                Vec3 pB0 = section_point(B, x0, ord(B, x0));
                Vec3 pB1 = section_point(B, x1, ord(B, x1));

                Panel p;
                // nominal winding; fixed below to make the normal outward.
                p.q.c = {pA0, pB0, pB1, pA1};
                p.area = quad_area(p.q);
                Vec3 cen = quad_centroid(p.q);

                // Outward reference: from the camber line toward this panel.
                double xm = 0.5 * (x0 + x1);
                double camb = 0.5 * (geom::cst_upper(A.af, xm) +
                                     geom::cst_lower(A.af, xm));
                Vec3 inner = section_point(A, xm, camb);
                Vec3 outward = cen - inner;

                Vec3 nrm = quad_normal(p.q);
                if (nrm.dot(outward) < 0.0) {            // flip winding -> outward
                    p.q.c = {pA0, pA1, pB1, pB0};
                    nrm = quad_normal(p.q);
                }
                p.n = nrm;
                double sz = std::sqrt(std::max(p.area, 1e-12));
                p.cp = cen - p.n * (1.0e-4 * sz);        // nudge just inside
                p.strip = sp;
                p.surf = surf;
                p.jc = j;

                int idx = static_cast<int>(m.panels.size());
                m.panels.push_back(p);
                if (j == nc - 1) {                       // TE-most panel
                    if (surf) m.te_up[sp] = idx; else m.te_lo[sp] = idx;
                }
            }
        }
    };
    add_surface(0);   // lower
    add_surface(1);   // upper
    return m;
}

}  // namespace

// Chordwise spacing selector. Default (and any unrecognised value) -> the
// convergent half-cosine; "cosine" selects the legacy LE+TE clustering.
static bool half_cosine_spacing(const Config& cfg) {
    return cfg.gets("panel_chord_spacing", "halfcosine") != "cosine";
}

MeshStats mesh_debug(const WingGeometry& w, const Config& cfg) {
    int nc = cfg.geti("panel_chordwise", 10);
    Mesh m = build_mesh(w, nc, half_cosine_spacing(cfg));
    MeshStats st;
    st.n_panels = static_cast<int>(m.panels.size());
    st.n_strips = m.n_strips;
    st.nc = m.nc;
    // Count panels whose enforced normal actually points INWARD vs a robust
    // deep-interior reference (chord-line point at 35% chord at the panel span).
    int flipped = 0;
    for (const auto& p : m.panels) {
        st.wetted_area += p.area;
        Vec3 cen = quad_centroid(p.q);
        // interior ref at the same span as the panel centroid, 35% chord, ord 0
        Station ref; ref.y = cen.y; ref.x_le = 0; ref.chord = 1; ref.twist = 0;
        // approximate: use the panel centroid projected toward the section mid.
        Vec3 interior(cen.x, cen.y, 0.0);   // chord-line height ~ 0 (thin section)
        // crude but robust: inward if normal points toward z=0 from above/below
        double zsign = (p.surf ? 1.0 : -1.0);   // upper should have n.z>0
        if (p.n.z * zsign < 0.0) ++flipped;
        (void)interior;
    }
    st.min_outward = static_cast<double>(flipped);   // reuse field: # inward normals
    // Self doublet at a representative panel (should be ~ -0.5).
    if (!m.panels.empty()) {
        const Panel& p = m.panels[m.panels.size() / 2];
        st.self_doublet = doublet_potential(p.q, p.cp);
    }
    return st;
}

// ---- Morino system: assembly, Kutta, factorisation cache -----------------
namespace {

struct PanelSystem {
    std::uint64_t sig = 0;
    double wake_alpha = -999.0;
    bool valid = false;
    Mesh mesh;
    Eigen::MatrixXd B;                          // source influence (N x N)
    Eigen::PartialPivLU<Eigen::MatrixXd> lu;    // factored doublet+Kutta matrix
    std::vector<Vec3> n;
    double S_ref = 1.0;
    int N = 0;
    double cond_C = 0.0;      // 2-norm condition number of C (debug only)
    double min_cp_gap = 0.0;  // smallest collocation-point separation (debug only)
};
thread_local PanelSystem t_sys;

std::uint64_t geom_sig_panel(const WingGeometry& w, int nc, bool half_cosine) {
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&h](double d) {
        std::uint64_t b; std::memcpy(&b, &d, sizeof(b));
        h = (h ^ b) * 1099511628211ull;
    };
    mix(w.root_chord); mix(w.tip_chord); mix(w.semi_span);
    mix(w.le_sweep); mix(w.washout); mix(w.le_bow); mix(w.te_bow);
    for (const auto& sec : w.sections) {
        mix(sec.te_thick);
        for (double v : sec.wu) mix(v);
        for (double v : sec.wl) mix(v);
    }
    mix(static_cast<double>(nc));
    mix(half_cosine ? 1.0 : 0.0);
    mix(static_cast<double>(w.stations.size()));
    for (const auto& s : w.stations) { mix(s.y); mix(s.chord); mix(s.x_le); mix(s.twist); }
    return h;
}

// A constant-strength doublet wake panel: TE edge (teA->teB) swept downstream
// along +x (body axis) for wake_len. Wound so the wake normal points +z (the
// lift sense), consistent with the upper-surface body-panel normals so the
// folded Kutta coefficient has the correct sign.
Quad wake_quad(const Vec3& teA, const Vec3& teB, const Vec3& dn) {
    Quad q;
    q.c = { teA, teA + dn, teB + dn, teB };
    return q;
}

// Mirror image of a panel across the y=0 plane. Reflecting (x,y,z)->(x,-y,z)
// reverses the winding, so we also reverse corner order to keep the image's
// outward normal physically correct (n_y flips, n_x/n_z preserved). The image
// carries the SAME source/doublet strength as the original (symmetric loading).
Quad image_quad(const Quad& q) {
    Quad r;
    for (int k = 0; k < 4; ++k) {
        const Vec3& c = q.c[3 - k];
        r.c[k] = Vec3(c.x, -c.y, c.z);
    }
    return r;
}

void build_system(const WingGeometry& w, int nc, double wake_len, const Vec3& wdir,
                  bool half_cosine) {
    t_sys.valid = false;
    Mesh m = build_mesh(w, nc, half_cosine);
    const int N = static_cast<int>(m.panels.size());
    if (N < 1) return;
    Vec3 dnw = wdir * wake_len;

    bool use_image = !getenv("PANEL_NO_IMAGE");
    std::vector<Quad> img(N);
    for (int j = 0; j < N; ++j) img[j] = image_quad(m.panels[j].q);

    // ---- Build base matrices (N x N) ----
    Eigen::MatrixXd C(N, N), B(N, N);
    for (int i = 0; i < N; ++i) {
        const Vec3& Pi = m.panels[i].cp;
        for (int j = 0; j < N; ++j) {
            // Self doublet term: use the analytic jump value (-1/2 for a
            // collocation point just inside the body), NOT the solid-angle
            // formula. On a WARPED panel the centroid-offset collocation point
            // straddles the panel's two defining triangles, so the numerical
            // solid angle flips branch (+1/2 vs -1/2) discontinuously; that
            // corrupted the diagonal whenever the wing had spanwise-varying
            // twist (washout) and collapsed the lift response. The +-1/2 jump
            // is exact and shape-independent (Katz & Plotkin, Morino).
            C(i, j) = (i == j) ? -0.5 : doublet_potential(m.panels[j].q, Pi);
            B(i, j) = source_potential(m.panels[j].q, Pi);
            if (use_image) {
                C(i, j) += doublet_potential(img[j], Pi);
                B(i, j) += source_potential(img[j], Pi);
            }
        }
    }

    // ---- Kutta: local wake doublet per strip ----
    for (int sp = 0; sp < m.n_strips; ++sp) {
        if (m.te_up[sp] < 0 || m.te_lo[sp] < 0) continue;
        Quad wq  = wake_quad(m.te_pt[sp], m.te_pt[sp + 1], dnw);
        Quad wqi = image_quad(wq);
        for (int i = 0; i < N; ++i) {
            double cw = doublet_potential(wq, m.panels[i].cp) +
                        (use_image ? doublet_potential(wqi, m.panels[i].cp) : 0.0);
            C(i, m.te_up[sp]) += cw;
            C(i, m.te_lo[sp]) -= cw;
        }
    }

    // Cache panel normals and reference area BEFORE moving the mesh out.
    t_sys.n.resize(N);
    for (int i = 0; i < N; ++i) t_sys.n[i] = m.panels[i].n;
    double S = 0.0;
    for (int sp = 0; sp < m.n_strips; ++sp) S += m.strip_chord[sp] * m.strip_dy[sp];
    t_sys.S_ref = (S > 0) ? S : 1.0;

    // ---- Optional conditioning diagnostics (R&D; guarded, GA path skips) ----
    // Tests the chordwise non-convergence hypothesis: as nc grows the cosine
    // chordwise clustering drives the upper/lower trailing-edge panels (and their
    // inside-nudged collocation points) near-coincident, producing near-identical
    // rows in C -> exploding condition number -> the CLa blow-up/sign-flip.
    if (std::getenv("PANEL_DEBUG_COND")) {
        Eigen::BDCSVD<Eigen::MatrixXd> svd(C);   // singular values only (fast, large N)
        const auto& sv = svd.singularValues();
        double smax = sv(0), smin = sv(sv.size() - 1);
        t_sys.cond_C = (smin > 0) ? smax / smin : std::numeric_limits<double>::infinity();
        double gap = std::numeric_limits<double>::infinity();
        for (int a = 0; a < N; ++a)
            for (int b = a + 1; b < N; ++b)
                gap = std::min(gap, (m.panels[a].cp - m.panels[b].cp).norm());
        t_sys.min_cp_gap = gap;
    } else {
        t_sys.cond_C = 0.0;
        t_sys.min_cp_gap = 0.0;
    }

    t_sys.B = std::move(B);
    t_sys.lu.compute(C);
    t_sys.mesh = std::move(m);       // mesh moved out LAST (nothing reads m after this)
    t_sys.N = N;
    t_sys.sig = geom_sig_panel(w, nc, half_cosine);
    t_sys.wake_alpha = std::atan2(wdir.z, wdir.x);
    t_sys.valid = true;
}

// Solve for the doublet strengths at a given freestream direction (V=1).
// Returns the per-panel mu and (out) the source flux for the watertight check.
Eigen::VectorXd solve_mu(const Vec3& Vinf, double& sigma_flux) {
    const int N = t_sys.N;
    if (N == 0) return Eigen::VectorXd();

    // Katz & Plotkin combined source/doublet, perturbation-potential interior
    // Dirichlet: sigma = -n.Vinf (freestream encoded in the source term; NO
    // separate phi_inf forcing — that double-counts).  C.mu = -B.sigma.
    Eigen::VectorXd sigma(N);
    sigma_flux = 0.0;
    for (int j = 0; j < N; ++j) {
        double s = -t_sys.n[j].dot(Vinf);
        sigma(j) = s;
        sigma_flux += s * t_sys.mesh.panels[j].area;
    }
    Eigen::VectorXd rhs = -(t_sys.B * sigma);
    return t_sys.lu.solve(rhs);
}
double trefftz_CL(const Eigen::VectorXd& mu, double& gamma_max,
                  std::vector<double>* gamma_out = nullptr) {
    const Mesh& m = t_sys.mesh;
    double sumGdy = 0.0;
    gamma_max = 0.0;
    if (gamma_out) gamma_out->assign(m.n_strips, 0.0);
    for (int sp = 0; sp < m.n_strips; ++sp) {
        double g = mu(m.te_up[sp]) - mu(m.te_lo[sp]);
        if (gamma_out) (*gamma_out)[sp] = g;
        sumGdy += g * m.strip_dy[sp];
        gamma_max = std::max(gamma_max, std::fabs(g));
    }
    // CL = 2 * sum(Gamma * dy) / (V * S_ref), V = 1, full span already summed.
    return 2.0 * sumGdy / t_sys.S_ref;
}

// DIAGNOSTIC (not the production x_np): neutral point from the chordwise
// doublet distribution. D(x) = mu_upper(x) - mu_lower(x) along a strip is the
// circulation accumulated from the LE (D(LE)=0, D(TE)=Gamma_strip), so the
// strip's chordwise centre of pressure follows from integration by parts
// (dL ~ dD acts at x):  x_cp = INT x dD / INT dD = x_te - (INT D dx)/D_te.
//
// This was evaluated as a replacement for the quarter-chord proxy that solve()
// uses, on the theory that the proxy ignored camber/reflex. It proved WORSE and
// is retained only for the xnp diagnostic: on a flat plate it sits ~1.5-3% MAC
// forward of the exact 0.25c (slow chordwise convergence of the constant-
// strength doublet), and on a reflexed flying-wing section it lands ~25% MAC AFT
// of AVL -- a reflexed airfoil's D(x) is non-monotonic (the tail unloads), so
// the small net circulation divided into large fore/aft local loading makes the
// CoP integral ill-conditioned. The quarter-chord proxy (the thin-airfoil
// aerodynamic centre) is instead within ~1.5-2% MAC of AVL on the same wing, so
// solve() keeps it. See scratch/xnp_probe.cpp for the side-by-side.
double neutral_point_load(const Eigen::VectorXd& mu, double fallback) {
    const Mesh& m = t_sys.mesh;
    const int nc = m.nc, ns = m.n_strips;
    if (nc < 1 || ns < 1) return fallback;
    // (surf, strip, jc) -> panel index, matching build_mesh's add order
    // (lower surface first, then upper; each loops strip then chordwise j).
    auto pid = [&](int surf, int sp, int jc) { return (surf * ns + sp) * nc + jc; };
    double sumL = 0.0, sumLx = 0.0;
    for (int sp = 0; sp < ns; ++sp) {
        double x_le = m.strip_xqc[sp] - 0.25 * m.strip_chord[sp];
        double x_te = x_le + m.strip_chord[sp];
        // Trapezoidal INT D dx with endpoints (x_le, 0) and (x_te, D_te); the
        // panel centroids supply the interior samples.
        double intD = 0.0, prevx = x_le, prevD = 0.0, Dte = 0.0, xc_last = x_le;
        for (int jc = 0; jc < nc; ++jc) {
            int up = pid(1, sp, jc), lo = pid(0, sp, jc);
            double D = mu(up) - mu(lo);
            double xc = 0.5 * (quad_centroid(m.panels[up].q).x +
                               quad_centroid(m.panels[lo].q).x);
            intD += 0.5 * (D + prevD) * (xc - prevx);
            prevx = xc; prevD = D; Dte = D; xc_last = xc;
        }
        intD += Dte * (x_te - xc_last);            // flat close to the camber TE
        if (std::fabs(Dte) > 1e-12) {
            double xcp = x_te - intD / Dte;
            double L = Dte * m.strip_dy[sp];       // strip lift ~ Gamma * dy
            sumL += L; sumLx += L * xcp;
        }
    }
    return (std::fabs(sumL) > 1e-12) ? sumLx / sumL : fallback;
}

void ensure_system(const WingGeometry& w, const Config& cfg, double alpha) {
    int nc = cfg.geti("panel_chordwise", 10);
    bool half = half_cosine_spacing(cfg);
    double wlen = cfg.getd("panel_wake_chords", 20.0) * w.root_chord;
    Vec3 wdir(std::cos(alpha), 0.0, std::sin(alpha));   // freestream-aligned wake
    bool geom_ok = t_sys.valid && t_sys.sig == geom_sig_panel(w, nc, half);
    // Frozen-wake mode: once a system is built for this geometry, an alpha change
    // reuses the existing factorization instead of rebuilding the dense N*N AIC.
    // The RHS in solve() is recomputed from the live alpha regardless, so a frozen
    // wake + perturbed freestream yields the exact partial derivative used by the
    // trim Jacobian -- without the per-FD-step rebuild (the GA hot-path cost
    // driver). Off by default; trim turns it on only for its Jacobian probes.
    bool freeze = cfg.geti("panel_freeze_wake", 0) != 0;
    bool stale = !geom_ok ||
                 (!freeze && std::fabs(alpha - t_sys.wake_alpha) > 1e-9);
    if (stale) build_system(w, nc, wlen, wdir, half);
}

// Span efficiency from the computed loading via a Glauert sine-series fit
// (odd modes, symmetric loading): Gamma(theta) = sum C_n sin(n theta),
// y = -s cos(theta). e = C_1^2 / sum_n (n C_n^2). e -> 1 for elliptic loading.
double span_efficiency(const std::vector<double>& gamma) {
    const Mesh& m = t_sys.mesh;
    int K = m.n_strips;
    if (K < 1) return 1.0;
    // Reconstruct the full-span loading by mirroring the half (symmetric).
    const double s = m.strip_y.back() + 0.5 * m.strip_dy.back();   // ~ semi-span
    if (s <= 0) return 1.0;
    std::vector<double> Y, G;
    for (int k = K - 1; k >= 0; --k) { Y.push_back(-m.strip_y[k]); G.push_back(gamma[k]); }
    for (int k = 0; k < K; ++k)      { Y.push_back( m.strip_y[k]); G.push_back(gamma[k]); }
    const int modes[5] = {1, 3, 5, 7, 9};
    const int M = 5;
    int KK = static_cast<int>(Y.size());
    Eigen::MatrixXd A(KK, M);
    Eigen::VectorXd b(KK);
    for (int k = 0; k < KK; ++k) {
        double yy = std::max(-0.999999, std::min(0.999999, Y[k] / s));
        double th = std::acos(-yy);
        for (int j = 0; j < M; ++j) A(k, j) = std::sin(modes[j] * th);
        b(k) = G[k];
    }
    Eigen::VectorXd C = A.colPivHouseholderQr().solve(b);
    double num = C(0) * C(0);
    double den = 0.0;
    for (int j = 0; j < M; ++j) den += modes[j] * C(j) * C(j);
    if (den <= 1e-30) return 1.0;
    double e = num / den;
    // Sanity clamp ONLY. The inviscid planar span efficiency is physically in
    // (0,1]. The lower bound is a NUMERICAL degeneracy guard -- it keeps CDi =
    // CL^2/(pi*e*AR) finite when a washout-dominated, near-zero-net-lift loading
    // makes the fundamental Glauert mode C_1 (hence the denominator's leading
    // term) vanish -- NOT a physical floor. The previous 0.30 floor was a real
    // bug: it clipped legitimately low-e designs and pinned e for ANY wing
    // probed at low alpha (where basic/twist loading dominates). At the trim
    // operating point, where the drag objective evaluates, e sits well inside
    // the band (e.g. a moderately tapered, washed-out knee trims to e ~ 0.92).
    constexpr double E_DEGEN = 0.05;
    return std::max(E_DEGEN, std::min(1.0, e));
}

// ---- Pitch-control derivatives (linear flap heuristic; mesh is undeflected,
//      same model as the VLM path so trim/delta_e behave identically) --------
}  // namespace

PanelSolveStats panel_solve_debug(const WingGeometry& w, const Config& cfg,
                                  double alpha) {
    ensure_system(w, cfg, alpha);
    PanelSolveStats st;
    st.n_panels = t_sys.N;
    if (!t_sys.valid) return st;
    Vec3 Vinf(std::cos(alpha), 0.0, std::sin(alpha));
    double flux = 0.0;
    Eigen::VectorXd mu = solve_mu(Vinf, flux);
    st.sigma_flux = flux;
    st.cl = trefftz_CL(mu, st.gamma_max);
    st.cond_C = t_sys.cond_C;
    st.min_cp_gap = t_sys.min_cp_gap;
    return st;
}

LoadingDump panel_loading_debug(const WingGeometry& w, const Config& cfg,
                                double alpha) {
    ensure_system(w, cfg, alpha);
    LoadingDump d;
    if (!t_sys.valid) return d;
    Vec3 Vinf(std::cos(alpha), 0.0, std::sin(alpha));
    double flux = 0.0;
    Eigen::VectorXd mu = solve_mu(Vinf, flux);
    double gmax = 0.0;
    std::vector<double> g;
    trefftz_CL(mu, gmax, &g);
    d.gamma = g;
    d.y = t_sys.mesh.strip_y;
    return d;
}

XnpDebug panel_xnp_debug(const WingGeometry& w, const Config& cfg, double alpha) {
    ensure_system(w, cfg, alpha);
    XnpDebug d;
    if (!t_sys.valid) return d;
    Vec3 Vinf(std::cos(alpha), 0.0, std::sin(alpha));
    double flux = 0.0;
    Eigen::VectorXd mu = solve_mu(Vinf, flux);
    double gmax = 0.0;
    std::vector<double> g;
    d.cl = trefftz_CL(mu, gmax, &g);
    // Quarter-chord proxy (the legacy x_np): Gamma-weighted strip quarter-chord.
    const Mesh& m = t_sys.mesh;
    double sg = 0.0, sgx = 0.0;
    for (int sp = 0; sp < m.n_strips; ++sp) {
        sg  += g[sp] * m.strip_dy[sp];
        sgx += g[sp] * m.strip_dy[sp] * m.strip_xqc[sp];
    }
    d.proxy = (std::fabs(sg) > 1e-12) ? sgx / sg : potential::wing_ac_x(w);
    d.load  = neutral_point_load(mu, d.proxy);
    return d;
}

// ---- Main aerodynamic solve (Morino panel + strip viscous coupling) -------
AeroState solve(const WingGeometry& w, const MassProps& mp,
                const viscous::Surrogate& surr, const Config& cfg,
                double alpha, double delta_e) {
    AeroState st;
    st.alpha = alpha;
    st.delta_e = delta_e;

    geom::ThinAirfoil ta = geom::thin_airfoil(w.sections.empty() ? Airfoil{} : w.sections[0]);
    const double AR = mp.AR;
    const double taper = (w.root_chord > 0) ? w.tip_chord / w.root_chord : 1.0;
    const double e_ref = potential::oswald_e(AR, taper);
    const double a_ref = potential::lift_curve_slope_3d(ta.cl_alpha, AR, e_ref);
    control::Derivs pc = control::compute(w, mp, a_ref, alpha, delta_e, cfg);

    // ---- Morino panel solve for the doublet distribution ----------------
    ensure_system(w, cfg, alpha);
    const Mesh& m = t_sys.mesh;
    Vec3 Vinf(std::cos(alpha), 0.0, std::sin(alpha));
    double flux = 0.0;
    std::vector<double> gamma;
    double panel_CL = 0.0, gmax = 0.0, x_np = potential::wing_ac_x(w), e_panel = e_ref;
    if (t_sys.valid && t_sys.N > 0) {
        Eigen::VectorXd mu = solve_mu(Vinf, flux);
        panel_CL = trefftz_CL(mu, gmax, &gamma);

        // Neutral point: circulation-weighted quarter-chord (sweep-aware).
        double sg = 0.0, sgx = 0.0, sgabs = 0.0;
        double xqc_lo = 1e300, xqc_hi = -1e300;
        for (int sp = 0; sp < m.n_strips; ++sp) {
            double wsp = gamma[sp] * m.strip_dy[sp];
            sg  += wsp;
            sgx += wsp * m.strip_xqc[sp];
            sgabs += std::fabs(wsp);
            xqc_lo = std::min(xqc_lo, m.strip_xqc[sp]);
            xqc_hi = std::max(xqc_hi, m.strip_xqc[sp]);
        }
        // The load-weighted centroid is only meaningful when there is net
        // circulation. At low diagnostic alpha a washed-out / reflexed wing can
        // carry near-zero NET load over large +/- local strip loads, so |sg|
        // collapses while sgx stays finite -- a bare |sg|>1e-12 guard then flings
        // x_np to +/-20 (a finite but DIVERGENT value that poisons the Pareto
        // sort). Two guards: (1) require net load to be a non-trivial fraction of
        // the gross absolute load, else keep the bounded geometric-AC fallback;
        // (2) a true load centroid of quarter-chord points must lie within their
        // own range, so clamp -- this stays exact for coherent loadings and
        // repairs the mixed-sign cancellation case unconditionally.
        if (sgabs > 1e-12 && std::fabs(sg) > 1e-3 * sgabs) {
            x_np = sgx / sg;
            if (xqc_hi >= xqc_lo)
                x_np = std::min(xqc_hi, std::max(xqc_lo, x_np));
        }
        e_panel = span_efficiency(gamma);
    }

    // Optional relaxed-wake pass: re-align the trailing wake with the local
    // downwash (induced angle) and re-solve, capturing the dominant free-wake
    // effect on CDi. Off by default; the reporting path enables it for the
    // incumbent only. ponytail: first-order wake tilt, not full Biot-Savart
    // rollup -- upgrade to multi-segment convection only if sub-0.1% CDi matters.
    if (cfg.geti("panel_relaxed_wake", 0) != 0 && t_sys.valid && t_sys.N > 0 && AR > 0) {
        int nc = cfg.geti("panel_chordwise", 10);
        bool half = half_cosine_spacing(cfg);
        double wlen = cfg.getd("panel_wake_chords", 20.0) * w.root_chord;
        int steps = std::max(1, cfg.geti("panel_wake_relax_steps", 3));
        for (int it = 0; it < steps && e_panel > 1e-6; ++it) {
            double a_ind = panel_CL / (PI * e_panel * AR);   // induced angle, rad
            Vec3 wdir(std::cos(alpha - a_ind), 0.0, std::sin(alpha - a_ind));
            build_system(w, nc, wlen, wdir, half);
            double fl = 0.0;
            Eigen::VectorXd mu = solve_mu(Vinf, fl);
            panel_CL = trefftz_CL(mu, gmax, &gamma);
            e_panel = span_efficiency(gamma);
        }
    }
    st.e = e_panel;
    st.x_np = x_np;
    st.static_margin = (mp.mac > 0) ? (x_np - mp.x_cg) / mp.mac : 0.0;

    // Total CL including the control-surface increment.
    st.CL = panel_CL + pc.CLde * delta_e;

    // Induced drag from the computed span efficiency (Trefftz-consistent).
    st.CDi = (AR > 0) ? st.CL * st.CL / (PI * e_panel * AR) : 0.0;

    // ---- Strip viscous coupling over the right half-wing ----------------
    const double V    = cfg.getd("v_cruise", V_CRUISE);
    const double cosL = std::cos(w.le_sweep);
    const double sweep_deg = w.le_sweep * RAD2DEG;
    const double xflow_deg = cfg.getd("sweep_crossflow_deg", 25.0);
    const double xflow_fac = cfg.getd("crossflow_factor", 1.15);
    const double Shalf = 0.5 * (mp.S_ref > 0 ? mp.S_ref : 1.0);

    std::vector<double> shape;  // rebuilt per strip from the lofted sections

    st.cl_local.assign(w.stations.size(), 0.0);
    double CDp_num = 0.0;
    bool tip_stall = false;
    double min_conf = 1.0;

    // Post-stall cap (M5): gate on cfg flag; accumulates capped-load centroid for x_np_high
    // ponytail: flag set by stability::trim() high-alpha pass only; zero cost on normal cruise solve
    bool post_stall_cap = cfg.geti("post_stall_cap", 0) != 0;
    double sumL_cap = 0.0, sumLx_cap = 0.0, sumLabs_cap = 0.0;
    bool root_capped = false, tip_capped = false;
    double xqc_cap_lo = 1e300, xqc_cap_hi = -1e300;

    // Strips ARE the half-wing strips (root..tip), one per [station i, i+1].
    int nstrip = m.n_strips;
    int nh = static_cast<int>(w.stations.size());
    for (int sp = 0; sp < nstrip; ++sp) {
        int i = sp;                                  // half-wing strip index
        double chord = m.strip_chord[sp];
        double dy    = m.strip_dy[sp];
        double t = (w.semi_span > 0) ? m.strip_y[sp] / w.semi_span : 0.0;
        double cl_i = (chord > 0 && !gamma.empty()) ? 2.0 * gamma[sp] / chord : 0.0;
        if (i < nh) st.cl_local[i] = cl_i;

        double Re      = RHO * V * cosL * chord / MU;
        double cl_norm = cl_i / (cosL * cosL);
        // strip section = mean of the two bounding lofted stations
        const Airfoil& afA = w.stations[sp].af;
        const Airfoil& afB = w.stations[sp + 1].af;
        shape.clear();
        for (std::size_t j = 0; j < afA.wu.size(); ++j)
            shape.push_back(0.5 * (afA.wu[j] + afB.wu[j]));
        for (std::size_t j = 0; j < afA.wl.size(); ++j)
            shape.push_back(0.5 * (afA.wl[j] + afB.wl[j]));
        double te_strip = 0.5 * (afA.te_thick + afB.te_thick);
        viscous::Polar pol = surr.query(shape, cl_norm, Re, te_strip);
        min_conf = std::min(min_conf, pol.confidence);

        // Accumulate capped-load centroid for the post-stall NP (M5 high-alpha pass)
        if (post_stall_cap && chord > 0) {
            double cl_cap = std::min(cl_norm, pol.cl_max);
            bool capped = (cl_norm > pol.cl_max && cl_i > 0.0);
            double L_cap = cl_cap * (cosL * cosL) * chord * dy;
            sumL_cap    += L_cap;
            sumLx_cap   += L_cap * m.strip_xqc[sp];
            sumLabs_cap += std::fabs(L_cap);
            xqc_cap_lo = std::min(xqc_cap_lo, m.strip_xqc[sp]);
            xqc_cap_hi = std::max(xqc_cap_hi, m.strip_xqc[sp]);
            if (capped) {
                if (t < 0.4) root_capped = true;
                if (t > 0.6) tip_capped  = true;
            }
        }

        double ramp = 0.0;
        if (sweep_deg > xflow_deg - 3.0) {
            double u = (sweep_deg - (xflow_deg - 3.0)) / 6.0;
            ramp = std::max(0.0, std::min(1.0, u));
        }
        double outer  = std::max(0.0, (t - 0.5) / 0.5);
        double factor = 1.0 + (xflow_fac - 1.0) * ramp * outer;
        CDp_num += pol.cd * factor * chord * dy;
        if (pol.clamped && cl_i > 0.0 && t > 0.6) tip_stall = true;
    }
    st.CDp = (Shalf > 0) ? CDp_num / Shalf : 0.0;
    st.CD  = st.CDi + st.CDp;
    st.tip_stall = tip_stall;
    st.polar_confidence = min_conf;
    // adverse yaw only on the cruise solve; the high-alpha pass discards it.
    if (!post_stall_cap)
        st.cn_da = control::adverse_yaw_cn_da(w, mp, surr, st.cl_local, a_ref, cfg);

    // Post-stall cap override: replace x_np and tip_stall with capped-load results (M5)
    if (post_stall_cap) {
        // Same guard as the Gamma-weighted NP block: require non-trivial net load
        if (sumLabs_cap > 1e-12 && std::fabs(sumL_cap) > 1e-3 * sumLabs_cap) {
            double xnp_cap = sumLx_cap / sumL_cap;
            if (xqc_cap_hi >= xqc_cap_lo)
                xnp_cap = std::min(xqc_cap_hi, std::max(xqc_cap_lo, xnp_cap));
            st.x_np = xnp_cap;
            st.static_margin = (mp.mac > 0) ? (xnp_cap - mp.x_cg) / mp.mac : 0.0;
        }
        // ponytail: tips-before-root is the pitch-up stall precursor; both capped = symmetric, benign
        st.tip_stall = tip_capped && !root_capped;
    }

    // ---- Pitching moment about CG ---------------------------------------
    double Cm0   = ta.cm_ac;
    double Cm_cl = st.CL * ((mp.x_cg - x_np) / (mp.mac > 0 ? mp.mac : 1.0));
    double Cm_de = pc.Cmde * delta_e;
    double q     = 0.5 * RHO * V * V;
    double T     = st.CD * q * (mp.S_ref > 0 ? mp.S_ref : 1.0);
    double zt    = cfg.getd("thrust_z_offset", 0.020);
    double denom = q * (mp.S_ref > 0 ? mp.S_ref : 1.0) * (mp.mac > 0 ? mp.mac : 1.0);
    double Cm_thr = (denom > 0) ? (T * zt) / denom : 0.0;
    st.CM = Cm0 + Cm_cl + Cm_de - Cm_thr;

    // ---- Hinge moment (Glauert thin-airfoil, worst-case pitch+roll) -----
    st.hinge_moment = pc.hinge_moment;

    // ---- Roll control / damping / helix (M6) ----------------------------
    st.cl_da      = pc.Cl_da;
    st.cl_p       = pc.Cl_p;
    st.roll_helix = pc.roll_helix;

    // x_np_high: cruise-NP fallback; stability::trim() overwrites this via a
    // high-alpha capped solve (M5) for the panel model.
    st.x_np_high = st.x_np;

    return st;
}

}  // namespace panel
}  // namespace aero
