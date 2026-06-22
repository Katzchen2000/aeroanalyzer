// linalg.h — self-contained dense linear algebra for AeroAnalyzer Pro.
//
// Rationale (see README "Numerical strategy"): the panel-method AIC matrix is a
// small, dense, non-symmetric system (~400-1600 unknowns at 20 stations). A
// direct partial-pivot LU is robust and fast there, so we ship our own rather
// than depend on Eigen to compile the scaffold. When you wire Eigen in for the
// Milestone-3 panel solver, swap solve_lu() for Eigen::PartialPivLU behind the
// same call site — nothing else has to change.
#pragma once
#include <vector>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <algorithm>

namespace aero {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm() const { return std::sqrt(dot(*this)); }
};

// Row-major dense matrix.
struct Matrix {
    std::size_t rows = 0, cols = 0;
    std::vector<double> a;
    Matrix() = default;
    Matrix(std::size_t r, std::size_t c) : rows(r), cols(c), a(r * c, 0.0) {}
    double& operator()(std::size_t i, std::size_t j) { return a[i * cols + j]; }
    double operator()(std::size_t i, std::size_t j) const { return a[i * cols + j]; }
};

struct LUFactorization {
    Matrix LU;
    std::vector<std::size_t> piv;
    bool singular = true;
};

// Factorize square matrix A in place, returning the factorization struct.
// piv[k] stores the row index that was swapped INTO row k at step k,
// so solve_lu_factored can replay the swaps in order during forward substitution.
inline LUFactorization factor_lu(Matrix A) {
    LUFactorization fact;
    fact.LU = std::move(A);
    const std::size_t n = fact.LU.rows;
    if (n == 0 || fact.LU.cols != n) return fact;
    fact.piv.resize(n);
    for (std::size_t i = 0; i < n; ++i) fact.piv[i] = i;  // identity — overwritten below

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        double maxv = std::fabs(fact.LU(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            double v = std::fabs(fact.LU(i, k));
            if (v > maxv) { maxv = v; p = i; }
        }
        if (maxv < 1e-14) return fact;  // singular
        fact.piv[k] = p;   // record swap partner (may equal k — no swap)
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(fact.LU(k, j), fact.LU(p, j));
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            double f = fact.LU(i, k) / fact.LU(k, k);
            fact.LU(i, k) = f;
            for (std::size_t j = k + 1; j < n; ++j) fact.LU(i, j) -= f * fact.LU(k, j);
        }
    }
    fact.singular = false;
    return fact;
}

// Solve LU x = b using a cached factorization.
inline bool solve_lu_factored(const LUFactorization& fact, std::vector<double> b, std::vector<double>& x) {
    if (fact.singular) return false;
    const std::size_t n = fact.LU.rows;
    if (b.size() != n) return false;
    
    // Forward substitution with pivoting
    for (std::size_t k = 0; k < n; ++k) {
        if (fact.piv[k] != k) {
            std::swap(b[k], b[fact.piv[k]]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            b[i] -= fact.LU(i, k) * b[k];
        }
    }
    
    // Backward substitution
    x.assign(n, 0.0);
    for (std::size_t ii = n; ii-- > 0;) {
        double s = b[ii];
        for (std::size_t j = ii + 1; j < n; ++j) s -= fact.LU(ii, j) * x[j];
        x[ii] = s / fact.LU(ii, ii);
    }
    return true;
}

// Backwards-compatible solve_lu
inline bool solve_lu(Matrix A, std::vector<double> b, std::vector<double>& x) {
    LUFactorization fact = factor_lu(std::move(A));
    return solve_lu_factored(fact, std::move(b), x);
}

// Least-squares solve of an overdetermined A x = b via normal equations
// (A^T A) x = A^T b. Adequate for the small, well-conditioned CST fit; for the
// panel solver use a pivoted factorization instead.
inline bool lstsq(const Matrix& A, const std::vector<double>& b,
                  std::vector<double>& x) {
    const std::size_t m = A.rows, n = A.cols;
    if (b.size() != m || m < n) return false;
    Matrix N(n, n);
    std::vector<double> rhs(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < m; ++k) s += A(k, i) * A(k, j);
            N(i, j) = s;
        }
        double r = 0.0;
        for (std::size_t k = 0; k < m; ++k) r += A(k, i) * b[k];
        rhs[i] = r;
    }
    return solve_lu(N, rhs, x);
}

}  // namespace aero
