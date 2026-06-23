// aero_neuralfoil.h - native C++ forward pass of the NeuralFoil "nn-large" model.
//
// NeuralFoil (Sharpe, arXiv:2503.16323; github.com/peterdsharpe/NeuralFoil) is a
// physics-informed neural surrogate of XFOIL. nn-large is a plain MLP:
//   25 -> 128 -> 128 -> 128 -> 128 -> 198  with Swish activations between the
// hidden layers and none on the output. We load the pretrained float32 weights
// (extracted from nn-large.npz to row-major (out,in) .bin files, see
// data/Neurafoilbin/) and run the forward pass + I/O glue in-process. This
// replaces the offline XFOIL polar table as the viscous surrogate's default
// engine: continuous, C-infinity smooth, no GPL dependency, valid across the
// whole Kulfan shape space.
//
// The glue (input construction, output decode, vertical-symmetry averaging) is
// transcribed from NeuralFoil's Python inference path and is the only part not
// contained in the weight files; it is pinned by the golden-vector unit test.
#pragma once
#include <array>
#include <string>
#include <vector>

namespace aero {
namespace nf {

// Decoded aerodynamics for one (airfoil, alpha, Re) operating point.
struct Aero {
    double cl = 0.0;
    double cd = 0.0;
    double cm = 0.0;
    double confidence = 0.0;  // sigmoid(output 0): NeuralFoil's own in-distribution
                              // trust signal in [0,1] (low => extrapolating)
};

class NeuralFoil {
public:
    // Loads the 10 weight/bias .bin files (net_{0,2,4,6,8}_{weight,bias}.bin)
    // from dir. Returns false and leaves the model unloaded if any file is
    // missing or the wrong size for nn-large.
    bool load(const std::string& dir);
    bool loaded() const { return loaded_; }

    // Raw 25 -> 198 forward pass (normalized inputs in, raw outputs out). Public
    // so the golden-vector test can exercise the matmul/activation in isolation.
    void forward(const float in[25], float out[198]) const;

    // Evaluate at a Kulfan airfoil (8 upper + 8 lower weights LE->TE, leading-edge
    // modification weight, trailing-edge thickness fraction) and operating point.
    // alpha in degrees. Runs the nominal pass plus the vertically-mirrored pass
    // and averages them (NeuralFoil's symmetry trick) for a sign-correct result.
    Aero eval(const std::array<double, 8>& upper,
              const std::array<double, 8>& lower,
              double le_weight, double te_thick,
              double alpha_deg, double Re,
              double ncrit, double xtr_upper, double xtr_lower) const;

private:
    struct Layer {
        int out = 0, in = 0;
        std::vector<float> W;  // row-major (out x in)
        std::vector<float> b;  // (out)
    };

    bool loaded_ = false;
    std::array<Layer, 5> layers_;  // net.0, net.2, net.4, net.6, net.8
};

}  // namespace nf
}  // namespace aero
