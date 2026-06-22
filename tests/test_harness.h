// Minimal dependency-free test harness (no GoogleTest needed).
#pragma once
#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <cmath>

namespace test {
struct Case { std::string name; std::function<void()> fn; };
inline std::vector<Case>& registry() { static std::vector<Case> r; return r; }
inline int& fails() { static int f = 0; return f; }
struct Reg {
    Reg(const std::string& n, std::function<void()> f) { registry().push_back({n, f}); }
};
inline void check(bool c, const char* e, const char* f, int l) {
    if (!c) { fails()++; std::cerr << "  FAIL " << f << ":" << l << "  " << e << "\n"; }
}
inline void check_near(double a, double b, double tol, const char* e,
                       const char* f, int l) {
    if (std::fabs(a - b) > tol || std::isnan(a) || std::isnan(b)) {
        fails()++;
        std::cerr << "  FAIL " << f << ":" << l << "  " << e << "  (" << a
                  << " vs " << b << ", tol " << tol << ")\n";
    }
}
}  // namespace test

#define TEST(n) static void n(); static test::Reg reg_##n(#n, n); static void n()
#define CHECK(c) test::check((c), #c, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tol) \
    test::check_near((a), (b), (tol), #a " ~= " #b, __FILE__, __LINE__)
