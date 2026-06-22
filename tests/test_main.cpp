#include "test_harness.h"

int main() {
    int n = 0;
    for (auto& c : test::registry()) {
        std::cout << "[ RUN ] " << c.name << "\n";
        c.fn();
        ++n;
    }
    if (test::fails() > 0) {
        std::cerr << "\n" << test::fails() << " check(s) FAILED across " << n
                  << " tests.\n";
        return 1;
    }
    std::cout << "\nAll " << n << " tests passed.\n";
    return 0;
}
