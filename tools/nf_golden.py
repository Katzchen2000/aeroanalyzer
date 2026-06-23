"""Generate golden CL/CD/CM values from the reference NeuralFoil (Python).

Run this on the machine where you extracted the weights (it needs the `neuralfoil`
and `aerosandbox` packages, same env you used for the .bin extraction):

    python tools/nf_golden.py

Paste the printed rows into tests/test_neuralfoil.cpp (the `ref[]` table in the
neuralfoil_golden test). The C++ port must match these to <1%.

The inputs below are IDENTICAL to the C++ golden test: an explicit 8+8 Kulfan
airfoil, leading_edge_weight=0, TE_thickness=0.004, Re=3e5, n_crit=9, free
transition (xtr_upper=xtr_lower=1), model_size="large".
"""
import numpy as np
import neuralfoil as nf
import aerosandbox as asb

UP = np.array([0.18, 0.20, 0.18, 0.16, 0.15, 0.16, 0.14, 0.12])
LO = np.array([-0.10, -0.12, -0.10, -0.08, -0.05, -0.02, 0.02, 0.05])
TE = 0.004
RE = 3.0e5
NCRIT = 9.0
ALPHAS = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0]

af = asb.KulfanAirfoil(
    upper_weights=UP,
    lower_weights=LO,
    leading_edge_weight=0.0,
    TE_thickness=TE,
)

print("# alpha      CL          CD          CM        confidence")
for a in ALPHAS:
    r = nf.get_aero_from_airfoil(
        af, alpha=a, Re=RE,
        model_size="large",
        n_crit=NCRIT,
        xtr_upper=1.0, xtr_lower=1.0,
    )
    # NeuralFoil returns 0-d numpy arrays; pull scalars.
    cl = float(np.ravel(r["CL"])[0])
    cd = float(np.ravel(r["CD"])[0])
    cm = float(np.ravel(r["CM"])[0])
    conf = float(np.ravel(r["analysis_confidence"])[0])
    print(f"{{{a:5.1f}, {cl:9.5f}, {cd:9.6f}, {cm:9.5f}}},  # conf={conf:.3f}")

# Also dump the raw 198-vector for one operating point so the bare forward pass
# (no symmetry averaging) can be checked in isolation if a mismatch appears.
# This calls the underlying net the same way main.py does for the nominal pass.
try:
    x = np.array([
        *UP, *LO, 0.0, TE * 50.0,
        np.sin(np.deg2rad(2 * 2.0)), np.cos(np.deg2rad(2.0)),
        1 - np.cos(np.deg2rad(2.0)) ** 2,
        (np.log(RE) - 12.5) / 3.5, (NCRIT - 9) / 4.5, 1.0, 1.0,
    ], dtype=np.float32)
    weights = np.load(
        __import__("os").path.join(__import__("os").path.dirname(nf.__file__),
                                   "nn_weights_and_biases", "nn-large.npz"))
    h = x
    for layer in [0, 2, 4, 6, 8]:
        W = weights[f"net.{layer}.weight"]
        b = weights[f"net.{layer}.bias"]
        h = W @ h + b
        if layer != 8:
            h = h / (1.0 + np.exp(-h))  # swish
    print("# raw forward out[0..5] =", " ".join(f"{v:.5f}" for v in h[:6]))
except Exception as e:  # noqa: BLE001
    print("# (raw forward dump skipped:", e, ")")
