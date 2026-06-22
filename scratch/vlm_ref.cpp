// vlm_ref.cpp — self-contained, textbook horseshoe VLM for a flat rectangular
// wing, to establish the correct Prandtl-ish lift slope independently of the
// project code. Uniform spanwise panels, single chordwise panel (lifting-line
// style), mirror handled by explicitly placing both half-wings' horseshoes.
#include <cstdio>
#include <cmath>
#include <vector>

static const double PI = 3.14159265358979323846;
static const double FOURPI = 4.0 * PI;

// z-velocity at (xc,yc,0) from unit horseshoe: bound x=xb from yL..yR,
// right trailing leg at yR (+1), left trailing leg at yL (-1), legs -> +inf x.
double w_horseshoe(double xc, double yc, double xb, double yL, double yR) {
    double dx = xc - xb, dyL = yc - yL, dyR = yc - yR;
    double r1 = std::sqrt(dx*dx + dyL*dyL), r2 = std::sqrt(dx*dx + dyR*dyR);
    double eps = 1e-12;
    double wb = 0.0;
    if (std::fabs(dx) > eps) wb = -1.0/(FOURPI*dx) * (dyL/r1 - dyR/r2);
    double wt = 0.0;
    if (std::fabs(dyR) > eps) wt += 1.0/(FOURPI*dyR) * (1.0 + dx/r2);
    if (std::fabs(dyL) > eps) wt -= 1.0/(FOURPI*dyL) * (1.0 + dx/r1);
    return wb + wt;
}

void run(double AR, int Nhalf) {
    double c = 0.2;
    double b = AR * c;                 // full span
    double s = 0.5 * b;                // semi span
    double S = b * c;
    // Uniform half-wing panels [0, s]; mirror via explicit left horseshoes.
    int N = Nhalf;
    std::vector<double> yL(N), yR(N), yc(N), xb(N), xc(N);
    double dy = s / N;
    for (int i = 0; i < N; ++i) {
        yL[i] = i*dy; yR[i] = (i+1)*dy; yc[i] = 0.5*(yL[i]+yR[i]);
        xb[i] = 0.25*c; xc[i] = 0.75*c;
    }
    // AIC: downwash at i from right horseshoe j + mirrored left horseshoe j.
    std::vector<std::vector<double>> A(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double wr = w_horseshoe(xc[i], yc[i], xb[j], yL[j], yR[j]);
            double wm = w_horseshoe(xc[i], yc[i], xb[j], -yR[j], -yL[j]);
            A[i][j] = wr + wm;
        }
    // Solve A G = rhs, rhs = -(alpha) (flat plate, twist=0, aL0=0). Use alpha=1 rad
    // then CL is linear so slope = CL/alpha; small-angle kernel is linear anyway.
    double alpha = 1.0 * PI/180.0;
    std::vector<double> rhs(N, -alpha), G(N);
    // Gaussian elimination
    std::vector<std::vector<double>> M = A;
    for (int k = 0; k < N; ++k) {
        int p = k; double mx = std::fabs(M[k][k]);
        for (int r = k+1; r < N; ++r) if (std::fabs(M[r][k])>mx){mx=std::fabs(M[r][k]);p=r;}
        std::swap(M[k], M[p]); std::swap(rhs[k], rhs[p]);
        for (int r = k+1; r < N; ++r) {
            double f = M[r][k]/M[k][k];
            for (int cc = k; cc < N; ++cc) M[r][cc] -= f*M[k][cc];
            rhs[r] -= f*rhs[k];
        }
    }
    for (int k = N-1; k >= 0; --k) {
        double sum = rhs[k];
        for (int cc = k+1; cc < N; ++cc) sum -= M[k][cc]*G[cc];
        G[k] = sum/M[k][k];
    }
    double sumG = 0.0; for (int i=0;i<N;++i) sumG += G[i]*dy;
    double CL = 4.0 * sumG / S;          // both halves, V=1
    double slope = CL/alpha;
    double prandtl = 2.0*PI*AR/(AR+2.0);
    std::printf("AR=%.1f N=%d  CL=%.4f slope=%.4f prandtl=%.4f err=%+.1f%%\n",
                AR, N, CL, slope, prandtl, 100.0*(slope-prandtl)/prandtl);
}

int main() {
    std::printf("--- independent textbook horseshoe VLM ---\n");
    for (double ar : {4.0,6.0,8.0,10.0}) run(ar, 40);
    return 0;
}
