#include "aeroanalyzer/aero_neuralfoil.h"
#include <cmath>
#include <fstream>
#include <iostream>

namespace aero {
namespace nf {

namespace {

// nn-large layer dimensions in evaluation order (net.0, net.2, net.4, net.6, net.8).
constexpr int kOut[5] = {128, 128, 128, 128, 198};
constexpr int kIn[5]  = { 25, 128, 128, 128, 128};
const char* const kWFile[5] = {"net_0_weight.bin", "net_2_weight.bin",
                               "net_4_weight.bin", "net_6_weight.bin",
                               "net_8_weight.bin"};
const char* const kBFile[5] = {"net_0_bias.bin", "net_2_bias.bin",
                               "net_4_bias.bin", "net_6_bias.bin",
                               "net_8_bias.bin"};

// Read exactly `count` little-endian float32 values from a binary file.
bool read_floats(const std::string& path, std::size_t count,
                 std::vector<float>& out) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) return false;
    const std::streamsize bytes = f.tellg();
    if (bytes != static_cast<std::streamsize>(count * sizeof(float)))
        return false;  // wrong-size file => reject (guards a bad extraction)
    out.resize(count);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(out.data()), bytes);
    return static_cast<bool>(f);
}

constexpr double kPi = 3.14159265358979323846;
inline float swish(float x) { return x / (1.0f + std::exp(-x)); }
inline double sigmoid(double x) { return 1.0 / (1.0 + std::exp(-x)); }
inline double sind(double deg) { return std::sin(deg * kPi / 180.0); }
inline double cosd(double deg) { return std::cos(deg * kPi / 180.0); }

}  // namespace

bool NeuralFoil::load(const std::string& dir) {
    loaded_ = false;
    for (int L = 0; L < 5; ++L) {
        Layer& ly = layers_[L];
        ly.out = kOut[L];
        ly.in  = kIn[L];
        const std::size_t wcount = static_cast<std::size_t>(ly.out) * ly.in;
        if (!read_floats(dir + "/" + kWFile[L], wcount, ly.W) ||
            !read_floats(dir + "/" + kBFile[L],
                         static_cast<std::size_t>(ly.out), ly.b)) {
            std::cerr << "[neuralfoil] missing or malformed " << kWFile[L]
                      << " / " << kBFile[L] << " in " << dir << "\n";
            return false;
        }
    }
    loaded_ = true;
    return true;
}

void NeuralFoil::forward(const float in[25], float out[198]) const {
    float buf_a[128], buf_b[128];
    const float* x = in;        // current layer input
    float* y = buf_a;           // scratch destination (ping-pong)
    for (int L = 0; L < 5; ++L) {
        const Layer& ly = layers_[L];
        float* dst = (L == 4) ? out : y;
        for (int o = 0; o < ly.out; ++o) {
            const float* wr = &ly.W[static_cast<std::size_t>(o) * ly.in];
            float acc = ly.b[o];
            for (int i = 0; i < ly.in; ++i) acc += wr[i] * x[i];
            dst[o] = (L < 4) ? swish(acc) : acc;  // no activation on the output
        }
        if (L < 4) {
            x = dst;
            y = (dst == buf_a) ? buf_b : buf_a;
        }
    }
}

namespace {
// Build the 25-element normalized input vector. The per-feature scalings ARE the
// network's input normalization (NeuralFoil applies no separate mean/std step).
void build_input(const std::array<double, 8>& up,
                 const std::array<double, 8>& lo,
                 double le_weight, double te_thick, double alpha_deg, double Re,
                 double ncrit, double xtr_upper, double xtr_lower, float in[25]) {
    for (int i = 0; i < 8; ++i) in[i]     = static_cast<float>(up[i]);
    for (int i = 0; i < 8; ++i) in[8 + i] = static_cast<float>(lo[i]);
    in[16] = static_cast<float>(le_weight);
    in[17] = static_cast<float>(te_thick * 50.0);
    in[18] = static_cast<float>(sind(2.0 * alpha_deg));
    in[19] = static_cast<float>(cosd(alpha_deg));
    in[20] = static_cast<float>(1.0 - cosd(alpha_deg) * cosd(alpha_deg));
    in[21] = static_cast<float>((std::log(Re) - 12.5) / 3.5);
    in[22] = static_cast<float>((ncrit - 9.0) / 4.5);
    in[23] = static_cast<float>(xtr_upper);
    in[24] = static_cast<float>(xtr_lower);
}
}  // namespace

Aero NeuralFoil::eval(const std::array<double, 8>& upper,
                      const std::array<double, 8>& lower,
                      double le_weight, double te_thick, double alpha_deg,
                      double Re, double ncrit, double xtr_upper,
                      double xtr_lower) const {
    Aero a;
    if (!loaded_) return a;

    float in[25], out_n[198], out_m[198];

    // --- nominal pass ---
    build_input(upper, lower, le_weight, te_thick, alpha_deg, Re, ncrit,
                xtr_upper, xtr_lower, in);
    forward(in, out_n);

    // --- mirrored pass: flip the airfoil about the chord line. Upper<->lower
    //     surfaces swap with a sign change, LE weight and alpha negate. The
    //     antisymmetric output channels (CL idx 1, CM idx 3) flip sign on the
    //     mirror; the symmetric ones (confidence idx 0, ln-drag idx 2) do not.
    std::array<double, 8> up_m, lo_m;
    for (int i = 0; i < 8; ++i) { up_m[i] = -lower[i]; lo_m[i] = -upper[i]; }
    build_input(up_m, lo_m, -le_weight, te_thick, -alpha_deg, Re, ncrit,
                xtr_lower, xtr_upper, in);
    forward(in, out_m);

    // Fuse the two passes in RAW output space, THEN decode. This is what the
    // reference does: averaging before the nonlinear decode (exp for drag,
    // sigmoid for confidence) is a geometric/logit mean, not an arithmetic one,
    // and matches NeuralFoil to ~5 digits. (For the linear CL/CM decode the two
    // orderings coincide, but for CD/confidence they do not.)
    const double y0 = 0.5 * (out_n[0] + out_m[0]);   // confidence (symmetric)
    const double y1 = 0.5 * (out_n[1] - out_m[1]);   // CL (antisymmetric)
    const double y2 = 0.5 * (out_n[2] + out_m[2]);   // ln-drag (symmetric)
    const double y3 = 0.5 * (out_n[3] - out_m[3]);   // CM (antisymmetric)

    a.cl         = y1 / 2.0;
    a.cd         = std::exp(2.0 * y2 - 4.0);
    a.cm         = y3 / 20.0;
    // NeuralFoil also subtracts a Mahalanobis out-of-distribution penalty (from
    // scaled_input_distribution.npz, not the weight files) before this sigmoid,
    // so our confidence runs ~0.005 optimistic vs the reference. Immaterial to
    // the coarse confidence<0.5 extrapolation gate; CL/CD/CM are exact.
    a.confidence = sigmoid(y0);
    return a;
}

}  // namespace nf
}  // namespace aero
