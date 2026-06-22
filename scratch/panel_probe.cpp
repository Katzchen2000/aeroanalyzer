// panel_probe.cpp — validate the Morino panel influence kernels against
// brute-force surface quadrature of the underlying integrals.
//   source_potential = -(1/4pi) * INT 1/r dS
//   solid_angle      =  INT (Q-P).n / r^3 dS   (used by doublet_potential)
#include "aeroanalyzer/aero_panel.h"
#include <cstdio>
#include <cmath>

using namespace aero;
using namespace aero::panel;

static const double PI4 = 4.0 * 3.14159265358979323846;

// Bilinear point on the quad at parameters (u,v) in [0,1]^2.
static Vec3 bilinear(const Quad& q, double u, double v) {
    Vec3 a = q.c[0] * ((1 - u) * (1 - v));
    Vec3 b = q.c[1] * (u * (1 - v));
    Vec3 c = q.c[2] * (u * v);
    Vec3 d = q.c[3] * ((1 - u) * v);
    return a + b + c + d;
}

// Brute-force INT 1/r dS and INT (Q-P).n/r^3 dS over the quad (fine midpoint).
static void brute(const Quad& q, const Vec3& P, double& src_int, double& sa) {
    const int M = 400;
    Vec3 n = quad_normal(q);
    src_int = 0.0; sa = 0.0;
    double h = 1.0 / M;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) {
            double u = (i + 0.5) * h, v = (j + 0.5) * h;
            Vec3 Q = bilinear(q, u, v);
            // local metric for dS
            Vec3 Qu = bilinear(q, u + 1e-6, v) - bilinear(q, u - 1e-6, v);
            Vec3 Qv = bilinear(q, u, v + 1e-6) - bilinear(q, u, v - 1e-6);
            double dS = (Qu * (1.0 / 2e-6)).cross(Qv * (1.0 / 2e-6)).norm() * h * h;
            Vec3 rv = Q - P;
            double r = rv.norm();
            if (r < 1e-9) continue;
            src_int += dS / r;
            sa += rv.dot(n) / (r * r * r) * dS;
        }
}

static void check(const char* name, const Quad& q, const Vec3& P) {
    double src_bf, sa_bf;
    brute(q, P, src_bf, sa_bf);
    double src_k = source_potential(q, P);     // = -(1/4pi) INT 1/r dS
    double src_k_int = -PI4 * src_k;           // recovered INT 1/r dS
    double dbl_k = doublet_potential(q, P);    // = -solid_angle/4pi
    double sa_k = -PI4 * dbl_k;                 // recovered solid angle
    std::printf("%-18s P=(%.2f,%.2f,%.2f)\n", name, P.x, P.y, P.z);
    std::printf("   source INT: kernel=%+.5f  brute=%+.5f  d=%+.1e\n",
                src_k_int, src_bf, src_k_int - src_bf);
    std::printf("   solidangle: kernel=%+.5f  brute=%+.5f  d=%+.1e   doublet=%+.5f\n",
                sa_k, sa_bf, sa_k - sa_bf, dbl_k);
}

// expose solid_angle for the probe
namespace aero { namespace panel { double solid_angle_pub(const Quad& q, const Vec3& P); } }

int main() {
    // Unit square in z=0, CCW about +z, centered at origin.
    Quad q{{ Vec3(-1,-1,0), Vec3(1,-1,0), Vec3(1,1,0), Vec3(-1,1,0) }};
    std::printf("area=%.4f normal=(%.2f,%.2f,%.2f)\n",
                quad_area(q), quad_normal(q).x, quad_normal(q).y, quad_normal(q).z);

    check("on-axis +z near", q, Vec3(0, 0, 0.5));
    check("on-axis -z near", q, Vec3(0, 0, -0.5));
    check("on-axis far",     q, Vec3(0, 0, 5.0));
    check("off-axis",        q, Vec3(0.3, -0.4, 0.7));
    check("far off-axis",    q, Vec3(2.0, 1.5, 3.0));

    // Far-field sanity: source INT -> A/R, solid angle -> A*cos/R^2.
    double R = 20.0, A = 4.0;
    check("very far +z", q, Vec3(0, 0, R));
    std::printf("   expect src INT ~ A/R = %.5f, solid angle ~ A/R^2 = %.5f\n",
                A / R, A / (R * R));
    return 0;
}

// Provide external linkage shim (solid_angle is file-static in the .cpp; we
// re-derive via doublet_potential instead). Defined to satisfy the declaration.
namespace aero { namespace panel {
double solid_angle_pub(const Quad& q, const Vec3& P) { return -PI4 * doublet_potential(q, P); }
}}
