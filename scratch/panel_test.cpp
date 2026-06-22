#include <iostream>
#include <vector>
#include <cmath>

struct Vec3 {
    double x, y, z;
    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
    Vec3 cross(const Vec3& o) const { return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x}; }
    double norm() const { return std::sqrt(dot(*this)); }
};

double solid_angle_triangle(const Vec3& r1, const Vec3& r2, const Vec3& r3) {
    double n1 = r1.norm(), n2 = r2.norm(), n3 = r3.norm();
    double num = r1.dot(r2.cross(r3));
    double den = n1*n2*n3 + r1.dot(r2)*n3 + r2.dot(r3)*n1 + r3.dot(r1)*n2;
    return 2.0 * std::atan2(num, den);
}

double integral_1_over_r(const std::vector<Vec3>& pts, const Vec3& P, const Vec3& n) {
    double I = 0.0;
    int N = pts.size();
    for (int i = 0; i < N; ++i) {
        Vec3 p1 = pts[i];
        Vec3 p2 = pts[(i+1)%N];
        Vec3 r1 = p1 - P;
        Vec3 r2 = p2 - P;
        double d = (p2 - p1).norm();
        if (d < 1e-14) continue;
        double num = r1.norm() + r2.norm() + d;
        double den = r1.norm() + r2.norm() - d;
        if (den < 1e-14) den = 1e-14;
        double term1 = (r1.cross(r2).dot(n)) / d;
        I += term1 * std::log(num / den);
    }
    double z = (pts[0] - P).dot(n);
    double omega = 0.0;
    for (int i = 1; i < N-1; ++i) {
        omega += solid_angle_triangle(pts[0]-P, pts[i]-P, pts[i+1]-P);
    }
    I -= std::abs(z) * std::abs(omega); // wait, omega could be negative, so abs(omega)
    return I;
}

int main() {
    std::vector<Vec3> panel = {
        {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}
    };
    Vec3 n = {0, 0, 1};
    
    Vec3 P1 = {0, 0, 1}; // Point above
    Vec3 P2 = {0, 0, 0}; // Point on the panel
    
    double I1 = integral_1_over_r(panel, P1, n);
    double I2 = integral_1_over_r(panel, P2, n);
    
    std::cout << "I1 (z=1): " << I1 << "\n";
    std::cout << "I2 (z=0): " << I2 << "\n";
    
    // Analytic for square at z=0 center: 
    // integral dx dy / sqrt(x^2+y^2) over [-1,1]x[-1,1]
    // = 8 * log(1 + sqrt(2)) = 7.0514
    std::cout << "Analytic I2: " << 8.0 * std::log(1.0 + std::sqrt(2.0)) << "\n";
    
    return 0;
}
