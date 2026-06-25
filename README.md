# AeroAnalyzer Pro

A multi-objective optimizer for a 3D-printed **flying wing**. It evolves a population
of designs with **NSGA-II** and returns a **Pareto front** trading three objectives:

- **Drag** (N, trimmed level cruise at 15 m/s)
- **Mass** (kg, structure + payload)
- **Static-margin deviation** from a target 6-8% MAC band

It uses constraint-domination gates such as TE thickness, tip Reynolds floor,
spar/OML clearance, hinge moment, high-alpha neutral-point migration, tip-stall,
roll authority, hardware keep-out, section validity, adverse yaw, and static-margin
floors.

---

## Contents

1. [Milestone status](#milestone-status-read-this-first)
2. [Build it yourself](#build-it-yourself)
3. [Typical workflow & outputs](#typical-workflow)
4. [Architecture & data flow](#architecture--data-flow)
5. [Implementation details (per phase)](#implementation-details)
6. [Objectives & constraints](#objectives--constraints)
7. [Airfoil seeds](#airfoil-seeds-ancestry-injection)
8. [Viscous engine: NeuralFoil](#viscous-engine-neuralfoil-runtime-default)
9. [Numerical strategy](#numerical-strategy-deliberate-choices)
10. [Performance & cost](#performance--cost)
11. [Configuration reference](#configuration-reference)
12. [Validation & tests](#validation--tests)
13. [Known limitations / approximations](#known-limitations--approximations)
14. [Roadmap](#roadmap)
15. [Equations reference](#equations-reference)

---

## Milestone status (read this first)

| Area | Status |
|------|--------|
| Geometry | Implemented and tested: CST sections, root/tip lofting, seeding, tip TE closure |
| Mass properties | Implemented and tested: sectional integration, CG, spar and hardware clearance outputs |
| Stability and trim | Implemented and tested: Newton trim, SM objective, high-alpha watchdogs |
| NSGA-II | Implemented and tested: fronts, crowding, SBX/PM, constraint-domination |
| Viscous engine | NeuralFoil is the default; table and analytic backends remain fallbacks |
| 3D aerodynamics | Morino panel solver is the default model; VLM remains the fallback |
| Control | Mode-aware roll authority, hinge moments, and adverse-yaw derivative are computed |
| Dynamic stability | Dutch-roll damping ratio and phugoid damping ratio computed per design via strip-theory Cn_β/Cn_r + Lanchester phugoid; reported in pareto.csv; opt-in cv gate |
| Release gates | Safety gates and surrogate-confidence gate are on by default |

## Current state

AeroAnalyzer is past the original prototype milestones. The live architecture is a 26-gene root/tip CST flying-wing optimizer using the Morino panel solver and NeuralFoil by default. Section shape lofts continuously from a root CST profile to a tip CST profile — the planform outline (LE/TE) is a straight swept taper, but the cross-sections are organically varying. The remaining open work is larger features: richer planform geometry (curved LE/TE, 3+ spanwise sections) and UX (visualization, checkpoint/restart, sensitivity sweeps).

3D aerodynamics defaults to the Morino Dirichlet panel solver (aero_model = panel, src/aero_panel.cpp): a constant-strength doublet/source method with strip viscous coupling. The earlier finite-wing analytic model remains available as aero_model = vlm behind the same solve signature.

## Build it yourself

Two no-CMake build scripts. Both need **Eigen** (header-only) — the default Morino panel
solver uses it for its dense AIC factorization, so it is a hard dependency of the core lib.
**`build_mingw.ps1` is the primary path** (MSYS2 `g++`):

```powershell
# from the project root (needs MSYS2 mingw-w64 g++ + Eigen; paths set at the top of the script)
powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test   # build + run the unit-test gates
powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Run    # build + run the optimizer
```

`build.ps1` (MSVC `cl.exe`, no CMake) is the alternate path and builds the same full app
(panel solver included), given the Eigen path set at the top of the script:

```powershell
powershell -ExecutionPolicy Bypass -File build.ps1 -Test -Run -OpenMP   # flags combine
```

Either script builds two executables into `build\`: `aeroanalyzer.exe` (optimizer) and
`unit_tests.exe`.

> **Editing `build*.ps1`:** keep them **pure ASCII**. PowerShell 5.1 reads non-BOM files as
> ANSI, and a UTF-8 em-dash decodes to a curly quote that PS treats as a string
> delimiter — which silently breaks the parse. (This bit us once already.)
>
> **Eigen + optimization level:** Eigen 5.0.0 fails to *link* at `-O0` with this MinGW
> toolchain; build at `-O2` (what `build_mingw.ps1` uses) for anything that links Eigen.

**CMake alternative** (after `winget install Kitware.CMake`): Eigen is wired via
`find_package` + `AERO_HAVE_EIGEN`.

```powershell
cmake -B build -S .
cmake --build build --config Release
ctest --test-dir build -C Release
```

Apart from Eigen (panel solver) the runtime is dependency-light — self-contained linear
algebra in `include/aeroanalyzer/linalg.h` backs the VLM path. OpenMP is auto-detected and
used only to parallelize candidate evaluation (the panel cache is `thread_local`, so the
parallel run stays byte-identical to serial).

---

## Typical workflow

```powershell
# 1. (optional) drop reflexed airfoil .dat files into data\airfoils\
# 2. build everything (the NeuralFoil viscous engine loads at runtime, no offline step)
powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test
# 3. run the optimizer
build\aeroanalyzer.exe config\baseline.cfg
# 4. inspect out\pareto.csv and the AVL decks (out\*.avl)
```

Outputs land in `out\`:

- **`pareto.csv`** — every feasible front-0 design, one row each. Columns:
  `idx, drag_N, mass_kg, sm_dev, static_margin, span_m, AR, root_c, tip_c, sweep_deg,
  washout_deg, CL, CD, hinge_kgcm, roll_helix, mode, dutch_roll_zeta, phugoid_zeta`.
- **`min_drag.avl`, `min_mass.avl`, `knee.avl`** (+ matching `.dat` section files) —
  AVL decks for the three incumbents (the two single-objective extremes and the
  normalized-L2 "knee" of the front), ready to open in **AVL** for an independent
  cross-check. The knee is the most balanced compromise design.

---

## Architecture & data flow

Modular by design so each phase is independently testable (a "monolithic main" would
make the panel solver impossible to unit-test). `main.cpp` is a thin driver.

```
                 config/baseline.cfg
                         │
   data/airfoils/*.dat ──┤      data/Neurafoilbin/*.bin
   + NACA seeds          │              │  (NeuralFoil nn-large weights)
                         ▼              ▼
                    ┌──────────┐   ┌───────────┐
                    │  seeds   │   │ aero_     │
                    │  (CST    │   │ viscous   │  NeuralFoil per (shape, Re)
                    │  fit +   │   │ Surrogate │  (Table/Analytic fallbacks)
                    │  bounds) │   └─────┬─────┘
                    └────┬─────┘         │
                         │  seed genomes │
                         ▼               │
   ga::NSGA2  ──────►  Evaluator ────────┤
   (fronts, crowding,    │ per candidate:│
    SBX, mutation,       ▼               ▼
    constraint-dom)   geom.decode → loft → massprops → stability.trim
                                                           │ (calls aero_potential.solve,
                                                           │  which queries the surrogate)
                                                           ▼
                                          objectives {drag, mass, SMdev} + cv
                         │
                         ▼
                 out/pareto.csv + AVL decks (avl_export)
```

Per-candidate evaluation lives in `Evaluator::run()` (`src/evaluate.cpp`) and returns an
`EvalResult` (geometry + mass + aero + objectives + constraint violation). The GA only
ever touches `Candidate` (genes + objectives + `cv`); everything physical is hidden
behind the evaluator. The surrogate is `const`/read-only after load, which keeps the
OpenMP parallel-evaluation loop race-free.

### Project layout

```
include/aeroanalyzer/   public headers (linalg, engine_core, geom, massprops,
                        aero_potential, aero_viscous, stability, ga, evaluate,
                        airfoil_io, seeds, avl_export, config)
src/                    implementations (one .cpp per header that needs one)
app/main.cpp            optimizer driver (thin)
tools/nf_golden.py      Python oracle to regenerate the NeuralFoil golden test vectors
tests/                  dependency-free unit tests + test_harness.h
config/baseline.cfg     all tunables (key = value)
data/airfoils/          drop seed .dat files here
data/Neurafoilbin/      NeuralFoil nn-large weight .bin files
data/surrogates/        pre-computed polar_coeffs.csv (Table-backend fallback)
build/                  objects + the two .exe (created by build.ps1)
out/                    pareto.csv + AVL decks (created at run)
```

---

## Implementation details

### Design genome (26 genes)

The GA optimizes a fixed-length real vector decoded by geom::decode(). Bounds come from geom::default_genome(); CST and TE bounds are widened at startup to contain seed airfoils. The genome is now root/tip aware, so section shape can vary along the span.

| Idx | Gene | Meaning | Default bounds |
|----:|------|---------|----------------|
| 0 | root_chord_m | root chord (m) | 0.18 - 0.35 |
| 1 | taper_ratio | tip/root chord | 0.30 - 0.90 |
| 2 | semi_span_m | half span (m) | 0.45 - 0.80 |
| 3 | le_sweep_deg | leading-edge sweep (deg) | 8 - 30 |
| 4 | washout_deg | tip twist (deg; negative = washout) | -6 - 0 |
| 5 | battery_x_m | battery-box CG x (m), the trim handle | 0.00 - 0.22 |
| 6-9 | wu0..wu3 | root upper CST weights | ~0.06 - 0.34 |
| 10-13 | wl0..wl3 | root lower CST weights; aft weight drives reflex | ~-0.22 - 0.20 |
| 14 | te_frac | root trailing-edge thickness fraction | 0.002 - 0.010 |
| 15 | mode | <0.5 = elevon, >=0.5 = split control | 0 - 1 |
| 16 | cs_chord_frac | control-surface chord fraction | 0.15 - 0.35 |
| 17 | ail_span_frac | aileron inboard edge, fraction of semi-span | 0.40 - 0.80 |
| 18-21 | tip_wu0..tip_wu3 | tip upper CST weights | root bounds mirrored |
| 22-25 | tip_wl0..tip_wl3 | tip lower CST weights | root bounds mirrored |

### Phase 1 - Geometry (geom.cpp)

- CST airfoils use the Kulfan form y(x) = C(x) * S(x) +/- x * te/2 with N1 = 0.5 and N2 = 1.0.
- fit_cst solves the least-squares fit used by ancestry injection and the CST round-trip tests.
- loft maps the wing across cosine-spaced span stations, interpolating root to tip CST weights while applying taper, sweep, and washout.
- evaluate.cpp applies te_thick_tip_frac after decode so the tip can close without a physical endplate.
- Planform integrals produce S_ref, MAC, MAC-LE x, span, and AR.

### Phase 2 - Mass properties (massprops.cpp)

- Section area and perimeter come directly from the local CST shape.
- Structural mass uses shell thickness plus a root-to-tip infill gradient.
- Motor, battery, avionics, and servos are point masses for total mass and CG.
- Spar clearance and hardware clearance are computed in massprops.cpp and enforced in evaluate.cpp.

### Phase 3 - Aerodynamics (aero_potential.cpp and aero_panel.cpp)

- aero_model = panel dispatches to the Morino Dirichlet panel solver.
- aero_model = vlm remains the fast analytic fallback for exploratory sweeps.
- The panel solver builds a dense AIC matrix once per geometry and reuses it through trim.
- The relaxed-wake pass is reserved for detailed incumbent re-evaluation; the GA hot path keeps it off.

### Phase 4 - Viscous engine (aero_viscous.cpp and aero_neuralfoil.cpp)

- NeuralFoil nn-large is the default section polar engine.
- The native C++ forward pass loads extracted .bin weights from data/Neurafoilbin.
- The legacy table backend and analytic backend remain fallbacks.
- The backend still presents the existing compact polar query contract to trim and drag code.

### Phase 5 - Stability and trim (stability.cpp)

- Trim is Newton-Raphson on alpha and elevator deflection for CL = W/(qS) and Cm = 0.
- Static-margin objective is zero inside the configured band and distance-to-band outside it.
- High-alpha panel re-evaluation supplies the pitch-up and tip-stall watchdogs.

### NSGA-II (ga.cpp)

- Constraint-domination follows Deb 2002.
- Fast nondominated sorting and crowding distance preserve front diversity.
- SBX crossover and bounded polynomial mutation are the variation operators.
- The initial population mixes seed elites, seeded hybrids, and random explorers.

## Objectives & constraints

Objectives are all minimized in src/evaluate.cpp.

| # | Objective | Definition |
|---|-----------|------------|
| OBJ_DRAG | Drag force | CD * q * S_ref at trimmed cruise |
| OBJ_MASS | Mass | Total structure plus payload |
| OBJ_SM | Stability | Static-margin deviation from sm_band_lo..sm_band_hi |

Currently enforced constraints are summed into cv, where cv = 0 means feasible.

| Constraint | Trigger | Weight |
|------------|---------|-------:|
| TE thickness clamp | te_frac * root_chord < te_thick_min_mm | 50 |
| Tip Reynolds floor | Re_tip < re_tip_min | 5 |
| High-alpha NP migration | x_np_high < x_cg | 60 |
| Tip-stall watchdog | tips cap before root | 20 |
| Trim convergence | Newton failed to trim | 100 |
| Static stability floor | SM < 0 | 25 |
| SM band floor | SM < sm_band_lo | sm_floor_penalty |
| Roll authority | roll_helix < roll_helix_min | 30 |
| Section validity | interior thickness < min_thickness_frac | 40 |
| Adverse yaw | adverse_yaw_penalty > 0 and cn_da > 0 | adverse_yaw_penalty |

Safety quantities now enforced by default:

- Spar/OML clearance is computed in massprops.cpp and penalized below the 1 mm target.
- Hinge moment is computed by control.cpp and penalized above hinge_moment_max.
- Hardware keep-out is computed in massprops.cpp and penalized when motor or avionics clearance goes negative.

## Airfoil seeds (ancestry injection)

The GA is seeded from real airfoils, then morphs them via the CST genes:

- **Uploaded .dat** - drop Selig or Lednicer files in data/airfoils (auto-detected, any chord scale; normalized to [0,1] on load). Reflexed sections are good flying-wing seeds.
- **Generated NACA 4-digit** - listed in seed_naca; the baseline uses 0012,2412,4412,4421,4418. Five-digit sections are not generated yet, so add those as .dat files instead.

Each seed is fit to CST as 4 upper + 4 lower weights. The CST gene bounds are widened by cst_bound_margin to contain every seed for both root and tip section genes. The initial population mixes pure seed elites, seeded root/tip hybrids with randomized planforms, and random explorers.

The hybrids jitter their root and tip airfoil genes around the seed by cst_seed_jitter, a fraction of each CST gene range. Without that jitter, the seeded half collapses onto only as many distinct airfoil shapes as there are seeds, which starves the CST search. The elites stay pure as known-good anchors; set cst_seed_jitter = 0 to recover exact-seed behavior.

---

## Viscous engine: NeuralFoil (runtime, default)

The default viscous engine is a native C++ forward pass of **NeuralFoil `nn-large`**
(Sharpe, [arXiv:2503.16323](https://arxiv.org/abs/2503.16323);
[github.com/peterdsharpe/NeuralFoil](https://github.com/peterdsharpe/NeuralFoil)) — a
physics-informed neural surrogate of XFOIL. It runs in-process (`src/aero_neuralfoil.cpp`,
namespace `aero::nf`) and **replaces the offline XFOIL polar table for runtime queries**:
continuous, C-infinity smooth, GPL-free, and valid across the whole Kulfan shape space (no
sparse-table plateaus or convex-hull guard needed — NeuralFoil's own `analysis_confidence`
output flags extrapolation instead).

- **Model:** a plain MLP `25 → 128 → 128 → 128 → 128 → 198` with **Swish** activations
  between hidden layers (none on the output). The pretrained float32 weights are extracted
  from `nn-large.npz` into row-major `(out, in)` `.bin` files under `data/Neurafoilbin/`
  (`net_{0,2,4,6,8}_{weight,bias}.bin`), loaded at startup.
- **Inputs (25):** 8 upper + 8 lower Kulfan weights, leading-edge weight, `TE·50`,
  `sin2α`, `cosα`, `sin²α`, `(lnRe−12.5)/3.5`, `(ncrit−9)/4.5`, `xtr_upper`, `xtr_lower`.
  These per-feature scalings *are* the input normalization (no separate mean/std step).
- **Outputs used:** `confidence = σ(y₀)`, `CL = y₁/2`, `CD = exp(2y₂−4)`, `CM = y₃/20`
  (the remaining 192 channels are boundary-layer arrays, unused here). Each call averages a
  nominal pass and a vertically-mirrored pass (NeuralFoil's symmetry trick).
- **Geometry bridge:** the GA's sections are 4+4 order-3 CST; NeuralFoil wants 8+8 order-7
  Kulfan + an LE weight. The 8-weight set is obtained by re-fitting via `geom::fit_cst(…,
  order=7,…)` (LE weight = 0), memoized per candidate.
- **Polar bridge:** to keep the existing `query(shape, cl, Re)` contract, the backend runs a
  small alpha sweep through the net and fits the same compact polar `cd = cd0 + k·cl²`
  (`cl_max/cl_min` from the swept envelope), so the trim/aero call sites are unchanged.

Select the engine with `viscous_backend` (`neuralfoil` | `table` | `analytic`) and point
`neuralfoil_dir` at the weights. If the weights are missing, it falls back to the
pre-computed polar table below, then to the analytic polar — the optimizer always runs.

**Verification:** `tests/test_neuralfoil.cpp` checks physical sanity (≈2π lift slope, a
positive drag bucket, up/down symmetry) and a golden cross-check against the reference Python
NeuralFoil (regenerate the reference with `python tools/nf_golden.py`).

## Table backend (read-only polar CSV) — *fallback*

`viscous_backend = table` reads a pre-computed polar table,
`data\surrogates\polar_coeffs.csv`, and interpolates the 5 coefficients across (shape, Re)
with an **inverse-distance scheme** (distances normalized per dimension by the trained hull)
and a **convex-hull guard** — queries outside the trained region are flagged (`clamped`),
never silently extrapolated.

The offline XFOIL generator that originally produced this CSV (`build_surrogate.exe` and the
in-process `aero::xfoil` solver) has been **retired** in favor of the NeuralFoil engine
above; the committed `polar_coeffs.csv` remains as a reproducible fallback and cross-check.

**`polar_coeffs.csv` format:** comment lines start with `#`; data rows are
`wu0,wu1,wu2,wu3,wl0,wl1,wl2,wl3,Re,cd0,k,cl_max,cl_min,cm0` (8 shape + Re + 5 coeffs).

If neither the NeuralFoil weights nor the table CSV are present, the runtime falls back to
the analytic polar with a warning, so the optimizer always runs.

---

## Numerical strategy (deliberate choices)

- **No `-ffast-math` / `/fp:fast`.** It implies finite-math-only, which lets the compiler
  delete the NaN/Inf watchdog checks the trim and solver rely on, and makes tolerances
  non-deterministic (and the GA non-reproducible). We use `/fp:precise` (`-O3
  -fno-math-errno` on GCC/Clang).
- **Direct LU, not BiCGSTAB.** The panel AIC matrix is small, dense, and not diagonally
  dominant; a partial-pivot LU (`linalg.h`, later Eigen `PartialPivLU`) is robust and
  reuses one factorization across the trim Newton's right-hand sides. A Jacobi-
  preconditioned iterative solver under fast-math is the opposite of that.
- **Reproducibility.** Evaluation is deterministic; RNG is seeded (`ga_seed`). The same
  config + seed yields a byte-identical `pareto.csv`, **even with OpenMP** (verified).

---

## Performance & cost

The per-candidate aerodynamic cost is dominated by **building the dense panel AIC matrix**
(transcendental influence kernels), rebuilt once per geometry. Everything else — LU
back-substitution, strip viscous coupling, mass/geometry — is comparatively free.

**Per-operation, single-thread** (Morino panel, measured via `scratch/panel_timing.cpp`):

| Operation | `panel_chordwise = 10` (default) | `= 6` |
|-----------|---------------------------------:|------:|
| Cold panel `solve()` (1 AIC build) | ~256 ms | ~92 ms |
| Trimmed `trim()` — **wake-frozen (default)** | ~266 ms | ~104 ms |
| Trimmed `trim()` — legacy (rebuild per Newton step) | ~3.4 s | ~1.5 s |
| VLM `trim()` (fallback model) | ~5 ms | ~5 ms |

**Full optimizer run** (`ga_pop 120 × ga_generations 150` = 18 000 trims, 12 threads):

| Model | Projected wall-clock |
|-------|---------------------:|
| **Panel, nc = 10 (default)** | **~13 min** |
| Panel, nc = 6 | ~5 min |
| VLM fallback | ~10 s |

Build/test turnaround: a full `build_mingw.ps1 -Test` (clean compile of all `src/` + Eigen
at `-O2`, then the gates) is ~2–3 min, dominated by the Eigen template compile.

### The wake-freeze optimization

`stability::trim()` solves `(α, δ_e)` with a Newton iteration whose per-iteration residual
**and** finite-difference Jacobian each call `solve()`. The panel AIC cache is keyed on the
wake angle, so every α change used to force a full rebuild — making a trim ~13× the cost of
one solve. But the wake **inclination** has negligible effect on attached loading (verified:
`dCL ~1e-11`, `dx_np ~15 µm` between a wake frozen at the initial α and one tracking α), so
`trim()` now **freezes the wake for the whole Newton solve** (`panel_trim_freeze_wake = 1`,
default): exactly **one AIC build per candidate**. The RHS still uses the live α, so the FD
derivative stays exact. This took the full GA from ~85 min to ~13 min (measured: 12T, nc=10, pop=120, 150 gens). Set
`panel_trim_freeze_wake = 0` (and `panel_trim_freeze_jac = 0`) to restore the legacy
rebuild-every-step behavior for comparison.

### Tuning levers (biggest first)

- **`panel_chordwise`** — the dominant cost knob; AIC assembly grows steeply with panel
  count. `10` is the converged default; `6` is ~2.6× faster and still passes the chordwise-
  convergence gate.
- **`aero_model = vlm`** — the analytic fallback is ~40× faster than the panel, but does not
  track washout (it reports a near-flat span efficiency); use it for quick exploratory sweeps,
  the panel for the final front.
- **OpenMP** — candidate evaluation parallelizes near-linearly; the `thread_local` panel
  cache keeps the parallel run byte-identical to serial.
- **`ga_pop` / `ga_generations`** — wall-clock scales linearly with their product.

---

## Configuration reference

All tunables live in `config\baseline.cfg` (`key = value`, `#` comments) — no recompile to
retune. A second path can be passed: `aeroanalyzer.exe path\to\my.cfg`.

| Key | Default | Meaning |
|-----|---------|---------|
| `v_cruise` | 15.0 | cruise speed (m/s) |
| `target_mass_aux` | 0.045 | avionics block mass (kg) |
| `seed_airfoils_dir` | `data/airfoils` | folder scanned for `*.dat` seeds |
| `seed_naca` | `0012,2412,4412` | generated NACA 4-digit seeds |
| `cst_bound_margin` | 0.04 | CST gene morph room widened around seeds |
| `cst_seed_jitter` | 0.25 | spread hybrids' CST genes around the seed (fraction of each gene's range); `0` = every seeded genome shares its seed's exact airfoil |
| `mass_motor` | 0.060 | motor mass (kg), at Y=0 TE keep-out |
| `mass_battery` | 0.210 | battery mass (kg), X-shiftable |
| `mass_servo_each` | 0.012 | per-servo mass (kg) |
| `material_density` | 1070 | ASA density (kg/m³) |
| `shell_thickness` | 0.0012 | wall thickness (m) |
| `infill_root` / `infill_tip` | 0.10 / 0.03 | infill gradient (fraction) |
| `spar_diameter` | 0.006 | spar cylinder diameter (m) |
| `spar_root_frac` | 0.15 | spar chordwise anchor (fraction) |
| `spar_straight` | 0 | 0 = follows swept chord; 1 = literal straight tube |
| `te_thick_min_mm` | 0.8 | trailing-edge thickness floor (mm) |
| `re_tip_min` | 100000 | tip Reynolds floor |
| `sm_band_lo` / `sm_band_hi` | 0.06 / 0.08 | target static-margin band |
| `sm_floor_penalty` | 30 | hard cv-gate weight forcing `SM ≥ sm_band_lo`; `0` = soft objective only |
| `thrust_z_offset` | 0.020 | thrust line above CG (m) |
| `hinge_moment_max` | 10 | fatal servo torque (kg-cm) |
| `roll_helix_min` | 0.05 | steady helix pb/2V floor (0 = disabled) |
| `aileron_deflect_max_deg` | 20 | max aileron throw (°) |
| `aileron_diff_ratio` | 4 | up/down differential ratio |
| `motor_diameter` | 0.028 | outrunner can diameter for keep-out (m) |
| `avionics_half_h` | 0.012 | avionics block half-height for keep-out (m) |
| `sweep_crossflow_deg` | 25.0 | LE-sweep crossflow threshold (°) |
| `crossflow_factor` | 1.15 | outer-span profile-drag multiplier |
| `ncrit` | 4.0 | viscous transition Ncrit (printed roughness) |
| `confidence_threshold` | 0.5 | soft cv gate for low-confidence NeuralFoil/table polars; `0` disables |
| `viscous_backend` | `neuralfoil` | viscous engine: `neuralfoil` (default), `table`, or `analytic` |
| `neuralfoil_dir` | `data/Neurafoilbin` | NeuralFoil nn-large weight `.bin` directory |
| `aero_model` | `panel` | 3D model: `panel` (Morino, default) or `vlm` (analytic fallback) |
| `high_alpha_deg` | 12 | angle of attack (°) for the high-α capped solve that evaluates NP migration and tip-stall; panel model only |
| `panel_chordwise` | 10 | chordwise panels per surface (dominant cost knob; see [Performance](#performance--cost)) |
| `panel_chord_spacing` | `halfcosine` | chordwise node spacing: `halfcosine` (convergent) or `cosine` (legacy) |
| `panel_wake_chords` | 20 | trailing wake length, in root chords |
| `panel_trim_freeze_wake` | 1 | freeze the wake across the trim Newton solve (~13× faster; 0 = legacy) |
| `panel_trim_freeze_jac` | 1 | freeze the wake across the trim FD Jacobian probes only |
| `avl_exe` | `tools/bin/avl352.exe` | standalone AVL path (cross-check; manual for now) |
| `ga_pop` | 120 | population size |
| `ga_generations` | 150 | generations |
| `ga_seed` | 1 | RNG seed (reproducibility) |
| `ga_eta_cx` / `ga_eta_mut` | 15 / 20 | SBX / polynomial-mutation distribution indices |
| `ga_pcx` | 0.9 | crossover probability |
| `ga_pmut` | −1 | mutation prob (<0 ⇒ 1/n_genes) |
| `n_stations` | 20 | spanwise stations |

---

## Validation & tests

`build_mingw.ps1 -Test` builds and runs the **dependency-free gates** (`tests/`, harness in
`test_harness.h`). Grouped:

- **Geometry** — CST round-trip, symmetric→zero camber, cambered→negative `α_L0`,
  reflex→`cm_ac` shift, rectangular planform area/AR, cosine lofting.
- **Mass** — solid-infill volume vs closed form, structure-only CG at ~0.42c, battery
  shift moves CG, spar clearance bounded.
- **Aero (Milestone-3 gates)** — run for **both** models: Prandtl lift slope `2π·AR/(AR+2)`,
  Oswald range, induced-drag consistency, lift increases with α, **trim converges** to
  `CL_req`/`Cm=0`, SM-objective band, surrogate hull clamp. Panel-specific: panel lift slope
  vs Prandtl, elliptic `e≈1`, chordwise convergence, washout loading sign/slope survival,
  quarter-chord neutral point, geometry-change cache invalidation.
- **Surrogate** — table load + IDW query + hull/cl clamps, analytic fallback when absent.
- **NeuralFoil** — physics sanity (monotonic lift, positive drag bucket, confidence in [0,1]), up/down symmetry for symmetric sections, golden cross-check against Python NeuralFoil (`tools/nf_golden.py`).
- **Airfoil I/O** — NACA 0012 max-thickness ≈ 12%, NACA 2412 cambered, NACA→CST round-trip,
  Selig + Lednicer parse.
- **Seeds** — NACA seeds load, bounds widen to contain seeds, seed genome decodes back to
  the seed shape, empty-seed graceful path.
- **NSGA-II** — constraint-domination rules, known-front non-dominated sort, crowding
  boundaries = ∞.

**AVL cross-check:** automated via `validate_avl.ps1`. Put `avl352.exe` in `tools\bin\`
(or pass `-AvlExe <path>`), then:

```powershell
powershell -ExecutionPolicy Bypass -File validate_avl.ps1
```

Runs AVL over the three incumbents at their panel-trimmed alpha, compares span efficiency
and neutral-point location, prints a PASS/FAIL table, and exits non-zero on any failure.
CL delta (~30%) is printed as INFO — the expected gap between viscous+trimmed panel and
inviscid clean AVL. Add `-Run` to regenerate the incumbents first.

**Reproducibility / OpenMP determinism** are both verified (same seed ⇒ identical
`pareto.csv`, serial and parallel).

---

## Known limitations / approximations

These are the remaining sharp edges rather than the model-of-record path:

- Control derivatives are strip-theory approximations intended for ranking, not final servo certification.
- Spar/OML clearance, hinge moment, and hardware keep-out are computed but not currently fatal gates.
- Dynamic stability (Dutch-roll + phugoid) is computed via single-DOF strip-theory / Lanchester approximations. Dutch-roll ignores roll coupling (Cl_β / dihedral effect); upgrade to a 2-DOF lateral state matrix if roll coupling matters at validation.
- Planform outline is a single straight swept taper (linear LE/TE). Section shape lofts organically from root to tip CST profiles. Curved planforms (non-straight LE/TE) and 3+ spanwise airfoil breakpoints are future geometry work.
- Battery dimensions are not explicit design variables yet.
- Battery dimensions are not explicit design variables yet.

## Roadmap

The ordered backlog below is the current specification for what remains.

### 1. Automate AVL validation ✓ DONE
- `validate_avl.ps1` drives AVL over the three incumbents; exits 0 if span efficiency (10%) and neutral-point location (3% MAC) are within tolerance. CL delta (~30%) is logged as INFO.
- Confirmed PASS: e within 3.5%, Xnp within 3.1 mm (gate = 5.6 mm) across all three incumbents.

### 2. Add dynamic stability ✓ DONE
- Dutch-roll damping ratio and phugoid damping ratio are computed per design.
- Strip-theory Cn_β (weathercock from LE sweep) and Cn_r (yaw rate damping); phugoid via Lanchester.
- Yaw inertia Izz computed from structural strips + point masses in massprops.
- Both metrics reported in pareto.csv. Opt-in cv gate via dynamic_stab_penalty, dutch_roll_zeta_min, phugoid_zeta_min (all default 0 = report-only, no feasibility impact).
- Calibration knobs: cn_beta_scale, cn_r_scale.

### 3. Expand geometry
- Add curved planform controls and/or more than two spanwise airfoil sections.
- Acceptance: genome, seed generation, decode, lofting, and tests cover the extra freedom.

### 4. Improve UX and analysis
- Add Pareto visualization.
- Add checkpoint/restart.
- Add sensitivity sweeps and export helpers.
- Acceptance: runs can be resumed and inspected without hand-editing CSVs.

---

## Equations reference

Quick index of the core relations:

| Quantity | Expression | Source |
|----------|------------|--------|
| CST surface | y = x^N1(1-x)^N2 * sum(w_i B_i(x)) +/- x*te/2 | geom.cpp |
| Induced drag | CDi = CL^2 / (pi e AR) | aero_potential.cpp |
| Section drag polar | cd = cd0 + k*cl^2, clamped to cl_min..cl_max | aero_viscous.cpp |
| Static margin | SM = (x_np - x_cg) / MAC | aero_potential.cpp |
| Trim | solve CL = W/(qS), Cm = 0 for alpha and delta_e | stability.cpp |
| Required CL | CL_req = W / (0.5 rho V^2 S) | stability.cpp |

---

AeroAnalyzer Pro - panel solver as model of record, NeuralFoil as the default viscous backend. Keep the tests green and treat the ordered backlog above as the source of truth for remaining work.
