# AeroAnalyzer Pro

A multi-objective optimizer for a 3D-printed **flying wing**. It evolves a population
of designs with **NSGA-II** and returns a **Pareto front** trading three objectives:

- **Drag** (N, trimmed level cruise at 15 m/s)
- **Mass** (kg, structure + payload)
- **Static-margin deviation** from a target 6–8% MAC band

…subject to hard constraints (TE thickness, tip Reynolds floor, spar/OML clearance,
servo hinge-moment limit, high-α neutral-point migration, tip-stall) handled with
Deb's **constraint-domination**.

---

## Contents

1. [Milestone status](#milestone-status-read-this-first)
2. [Build it yourself](#build-it-yourself)
3. [Typical workflow & outputs](#typical-workflow)
4. [Architecture & data flow](#architecture--data-flow)
5. [Implementation details (per phase)](#implementation-details)
6. [Objectives & constraints](#objectives--constraints)
7. [Airfoil seeds](#airfoil-seeds-ancestry-injection)
8. [Viscous surrogate](#viscous-surrogate-xfoil-offline)
9. [Numerical strategy](#numerical-strategy-deliberate-choices)
10. [Configuration reference](#configuration-reference)
11. [Validation & tests](#validation--tests)
12. [Known limitations / approximations](#known-limitations--approximations)
13. [Roadmap](#roadmap)
14. [Equations reference](#equations-reference)

---

## Milestone status (read this first)

| Phase | Status |
|------|--------|
| Geometry (CST, lofting, seeding) | ✅ implemented + tested |
| Mass properties (sectional integration, spar clearance) | ✅ implemented + tested |
| Stability & trim (Newton, SM objective, watchdogs) | ✅ implemented + tested |
| NSGA-II (fronts, crowding, SBX/PM, constraint-domination) | ✅ implemented + tested |
| Airfoil seeds (.dat + NACA) + ancestry injection | ✅ implemented + tested |
| Viscous surrogate (XFOIL-calibrated, shape-parameterized) | ✅ runtime + offline tool |
| **3D aerodynamics** | ✅ **Morino panel solver (default)** — AVL-validated; VLM is the fallback |

3D aerodynamics defaults to the **Morino Dirichlet panel solver** (`aero_model = panel`,
`src/aero_panel.cpp`): a constant-strength doublet/source method with strip viscous
coupling. It passes every analytic gate in `tests/test_aero.cpp` (lift slope
`2π·AR/(AR+2)`, elliptic induced drag, trim convergence, quarter-chord neutral point) and
is cross-checked against **AVL** on the `knee` / `min_drag` / `min_mass` decks — CLα within
~0.1–2.5 %, neutral point within a few % MAC. The earlier finite-wing analytic model
(Prandtl lift slope + strip coupling) remains available as the documented fallback via
`aero_model = vlm` in `aero_potential::solve()`, behind the same signature.

---

## Build it yourself

You have the VS2022 **Build Tools** (`cl.exe`) but no CMake — so the primary path needs
no CMake:

```powershell
# from the project root
powershell -ExecutionPolicy Bypass -File build.ps1 -Test    # build + run unit tests
powershell -ExecutionPolicy Bypass -File build.ps1 -Gen     # build + generate viscous surrogate
powershell -ExecutionPolicy Bypass -File build.ps1 -Run     # build + run the optimizer
powershell -ExecutionPolicy Bypass -File build.ps1 -OpenMP  # multicore GA evaluation
```

Flags combine, e.g. `build.ps1 -Test -Gen -Run`. `build.ps1` locates `vcvars64.bat` and
drives `cl.exe` directly. It builds three executables into `build\`:
`aeroanalyzer.exe` (optimizer), `build_surrogate.exe` (offline surrogate generator),
`unit_tests.exe`.

> **Editing `build.ps1`:** keep it **pure ASCII**. PowerShell 5.1 reads non-BOM files as
> ANSI, and a UTF-8 em-dash decodes to a curly quote that PS treats as a string
> delimiter — which silently breaks the parse. (This bit us once already.)

**CMake alternative** (after `winget install Kitware.CMake`):

```powershell
cmake -B build -S .
cmake --build build --config Release
ctest --test-dir build -C Release
```

The runtime binary has **zero external dependencies** (self-contained linear algebra in
`include/aeroanalyzer/linalg.h`). Eigen is wired as *optional* in CMake (`find_package`
+ `AERO_HAVE_EIGEN`) for the future panel solver but is not required to build today.
OpenMP is auto-detected and used only to parallelize candidate evaluation.

---

## Typical workflow

```powershell
# 1. (optional) drop reflexed airfoil .dat files into data\airfoils\
# 2. build everything and generate the viscous surrogate
powershell -ExecutionPolicy Bypass -File build.ps1 -Test -Gen
# 3. run the optimizer
build\aeroanalyzer.exe config\baseline.cfg
# 4. inspect out\pareto.csv and the AVL decks (out\*.avl)
```

Outputs land in `out\`:

- **`pareto.csv`** — every feasible front-0 design, one row each. Columns:
  `idx, drag_N, mass_kg, sm_dev, static_margin, span_m, AR, root_c, tip_c, sweep_deg,
  washout_deg, CL, CD, hinge_kgcm, mode`.
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
   data/airfoils/*.dat ──┤      data/surrogates/polar_coeffs.csv
   + NACA seeds          │              │  (built offline by build_surrogate.exe)
                         ▼              ▼
                    ┌──────────┐   ┌───────────┐
                    │  seeds   │   │ aero_     │
                    │  (CST    │   │ viscous   │  RBF/IDW over (shape, Re)
                    │  fit +   │   │ Surrogate │
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
tools/build_surrogate.cpp   offline XFOIL surrogate generator (separate exe)
tests/                  dependency-free unit tests (31 gates) + test_harness.h
config/baseline.cfg     all tunables (key = value)
data/airfoils/          drop seed .dat files here
data/surrogates/        generated polar_coeffs.csv
build/                  objects + the three .exe (created by build.ps1)
out/                    pareto.csv + AVL decks (created at run)
```

---

## Implementation details

### Design genome (16 genes)

The GA optimizes a fixed-length real vector decoded by `geom::decode()`. Bounds come from
`geom::default_genome()`; the CST and TE bounds are **widened at startup** to contain the
seed airfoils (see [Airfoil seeds](#airfoil-seeds-ancestry-injection)).

| Idx | Gene | Meaning | Default bounds |
|----:|------|---------|----------------|
| 0 | `root_chord_m` | root chord (m) | 0.18 – 0.35 |
| 1 | `taper_ratio` | tip/root chord | 0.30 – 0.90 |
| 2 | `semi_span_m` | half span (m) | 0.45 – 0.80 |
| 3 | `le_sweep_deg` | leading-edge sweep (°) | 8 – 30 |
| 4 | `washout_deg` | tip twist (°, −=washout) | −6 – 0 |
| 5 | `battery_x_m` | battery-box CG x (m) — the SM trim handle | 0.00 – 0.22 |
| 6–9 | `wu0..wu3` | upper CST (Bernstein) weights | ~0.06 – 0.34 |
| 10–13 | `wl0..wl3` | lower CST weights (aft drives **reflex**) | ~−0.22 – 0.20 |
| 14 | `te_frac` | trailing-edge thickness (fraction of chord) | 0.002 – 0.010 |
| 15 | `mode` | <0.5 = elevon, ≥0.5 = split control | 0 – 1 |

### Phase 1 — Geometry (`geom.cpp`)

- **CST (Kulfan) airfoil:** `y(x) = C(x)·S(x) + x·(te/2)`, class function `C = x^N1 (1−x)^N2`
  with `N1=0.5, N2=1.0` (round LE, sharp TE), shape `S` a Bernstein-polynomial sum over
  4 weights per surface. Independent upper/lower weights make **reflex** representable —
  essential for a tailless wing's pitch trim.
- **CST fit** (`fit_cst`) solves a small linear least-squares (normal equations via
  `linalg::lstsq`) for the weights given sampled `(x, z)` points. Used for ancestry
  injection and validated by a round-trip test.
- **Thin-airfoil coefficients** (`thin_airfoil`) integrate the camber-line slope to get
  zero-lift angle `α_L0`, aerodynamic-center moment `cm_ac` (positive for reflex → enables
  tailless trim), and `cl_α = 2π`.
- **Lofting** (`loft`) maps the section across **20 cosine-spaced** spanwise stations
  (clustered at root and tip), applying linear taper, **sheared sweep** (`x_le = y·tanΛ`,
  which preserves streamwise sections), and linear washout. Each station carries a strip
  width for trapezoidal integration.
- **Planform** integrals give `S_ref`, MAC, MAC-LE x, span, and AR from the stations.

### Phase 2 — Mass properties (`massprops.cpp`)

Cross-sectional lofting (**not** B-Rep — no boolean kernel, no voxels):

- Normalized section area `Â` and perimeter `P̂` integrated once from the CST shape.
- Per station: `enclosed = Â·c²`, `shell = P̂·c·t_shell`, and a **linear infill gradient**
  (root→tip) gives an equivalent solid area `shell + infill·(enclosed − shell)`. Spanwise
  **trapezoidal** accumulation yields structural volume → mass (×2 for both halves).
- **Point masses** (motor at root TE, battery at `battery_x`, avionics block over 20–50%
  span, 2 servos) are added with chordwise centroids to compute total mass and **CG**.
- **Analytic spar clearance:** at each station, clearance = local half-thickness at the
  spar's chordwise position − spar radius (negative = OML breach). The spar follows the
  swept 15%-local-chord line by default (`spar_straight = 0`); a literal straight tube
  (`spar_straight = 1`) breaches the OML on swept wings and is what made every early
  candidate infeasible.

### Phase 3 — Aerodynamics (`aero_potential.cpp`) — *VLM fallback*

> The default model is now the Morino panel solver (`src/aero_panel.cpp`, `aero_model =
> panel`); the analytic model below is the `aero_model = vlm` fallback, dispatched from the
> same `solve(w, mp, surr, cfg, alpha, delta_e)` signature.

- 3D lift slope `a = a0 / (1 + a0/(π·e·AR))` (Prandtl), `a0 = 2π`, span efficiency `e`
  from a taper-based Oswald estimate.
- `CL = a·(α + ½·washout − α_L0) + CL_δ·δ_e`; uniform induced angle `α_i = CL/(π·e·AR)`.
- **Strip viscous coupling** over the 20 stations: local `cl_i = cl_α·(α + twist_i − α_L0
  − α_i)`, queried against the surrogate at **sweep-normal** conditions (`Re = ρ·V·cosΛ·c/μ`,
  `cl_normal = cl_i/cos²Λ`). A **smooth crossflow penalty** multiplies outer-span profile
  drag once LE sweep exceeds ~25°.
- `CD = CDi + CDp`, `CDi = CL²/(π·e·AR)`, `CDp` = strip integral of section `cd`.
- **Pitching moment about CG:** `Cm = cm_ac + CL·(x_cg − x_np)/MAC + Cm_δ·δ_e − Cm_thrust`,
  with `x_np` = area-weighted quarter-chord (captures sweep) and a thrust-line Z-offset
  term. **Hinge moment** (servo torque, kg·cm) and a heuristic **high-α NP migration**
  estimate (swept tips unload, washout relieves) are produced for the watchdogs.

### Phase 4 — Viscous surrogate (`aero_viscous.cpp` + `tools/build_surrogate.cpp`)

See [Viscous surrogate](#viscous-surrogate-xfoil-offline). Compact polar
`cd = cd0 + k·cl²` with `cl_max/cl_min` clamp; runtime interpolates the 5 coefficients
across (shape, Re).

### Phase 5 — Stability & trim (`stability.cpp`)

- **Trim:** Newton–Raphson on `(α, δ_e)` solving `CL = W/(q·S)` and `Cm = 0`
  simultaneously. The 2×2 Jacobian is built by **finite differences** of `solve()` so it
  stays valid when the nonlinear panel solver swaps in. Step-clamped (≤5°) for robustness;
  `state.trimmed = false` if it fails to converge (which the constraints penalize hard).
- **SM objective:** `sm_objective(SM, lo, hi)` = 0 inside the band, else distance to the
  nearest edge (this is the minimized third objective).

### NSGA-II (`ga.cpp`)

- **Constraint-domination** (Deb 2002): feasible beats infeasible; among infeasible the
  smaller total violation wins; among feasible, ordinary Pareto dominance on the 3
  objectives.
- `fast_nondominated_sort` ranks into fronts; `crowding_distance` preserves diversity
  (boundary points get ∞). Selection is a binary tournament on (rank, crowding).
- **Operators:** SBX crossover (`eta_cx`) + bounded polynomial mutation (`eta_mut`,
  per-gene prob `1/n` by default). Elitist truncation: combine parents+offspring, fill the
  next generation by fronts then crowding.
- **Seeding:** the first ~50% of the initial population is taken from seed genomes (elites
  + hybrids), the rest random explorers.

---

## Objectives & constraints

**Objectives (all minimized)** — `src/evaluate.cpp`:

| # | Objective | Definition |
|---|-----------|------------|
| `OBJ_DRAG` | drag force | `CD · q · S_ref` (N) at trimmed cruise |
| `OBJ_MASS` | mass | total structure + payload (kg) |
| `OBJ_SM` | stability | deviation of static margin from the `[sm_band_lo, sm_band_hi]` band |

**Constraints** → summed into a single violation `cv` (0 = feasible), weighted so the
fatal ones dominate the merely marginal:

| Constraint (plan ref) | Trigger | Weight |
|-----------------------|---------|-------:|
| TE thickness clamp (§3) | `te_frac·root_chord < te_thick_min_mm` | 50 |
| Tip Reynolds floor (§3) | `Re_tip < re_tip_min` | 5 |
| Spar/OML clearance (§4) | clearance < 1 mm (continuous); ×30 if breach | 4 / 30 |
| Hinge-moment gate (§7) | `hinge_moment > hinge_moment_max` (fatal) | 40 |
| High-α NP migration (§7) | `x_np_high < x_cg` (pitch-up) | 60 |
| Tip-stall watchdog (§6) | outboard station stalled at cruise | 20 |
| Trim convergence | Newton failed to trim | 100 |
| Static stability floor | `SM < 0` (NP forward of CG) | 25 |

---

## Airfoil seeds (ancestry injection)

The GA is seeded from real airfoils, then morphs them via the CST genes:

- **Uploaded `.dat`** — drop Selig or Lednicer files in `data\airfoils\` (auto-detected,
  any chord scale; normalized to [0,1] on load). Reflexed sections (MH/EH families) are
  good flying-wing seeds.
- **Generated NACA 4-digit** — listed in `seed_naca` (e.g. `0012,2412,4412`); no file
  needed. *(5-digit like 23012 isn't generated yet — drop it as a `.dat` instead.)*

Each seed is fit to CST (4 upper + 4 lower weights). The CST gene bounds are **widened** by
`cst_bound_margin` to contain every seed (so seeds aren't clipped **and** the GA's morph
region matches the surrogate's trained hull). ~50% of the initial population is built from
seeds — 10% elites (seed shape + mid-box planform), 40% hybrids (seed shape + random
planform) — with the rest random explorers.

---

## Viscous surrogate (XFOIL, offline)

The GA reshapes airfoils, so a single fixed polar won't do. `build_surrogate.exe` samples
many CST shapes (seeds + a Latin-Hypercube DoE of `surrogate_samples` around them) × a
Reynolds grid (`surrogate_re`), calibrates a compact polar (`cd0, k, cl_max, cl_min, cm0`,
with `cd = cd0 + k·cl²`) for each, and writes `data\surrogates\polar_coeffs.csv`.

At runtime the surrogate interpolates those coefficients across (shape, Re) with an
**inverse-distance scheme** (distances normalized per dimension by the trained hull) and a
**convex-hull guard** — queries outside the trained region are flagged (`clamped`), never
silently extrapolated. The runtime binary never calls XFOIL.

Two modes (config `surrogate_mode`):

- **`synthetic`** (default) — analytic calibration, **no XFOIL needed**; lets you exercise
  the full table pipeline now. Coefficients still vary with shape (thickness/camber) and
  Re so interpolation is meaningful.
- **`xfoil`** — drives standalone `xfoil.exe` (set `xfoil_exe`, **no Python required**) per
  (shape, Re): writes a `.dat`, scripts `OPER/VPAR N/VISC Re/PACC/ASEQ`, parses the polar
  dump, and least-squares-fits the compact polar. Rows that fail to converge fall back to
  synthetic.

**`polar_coeffs.csv` format:** comment lines start with `#`; data rows are
`wu0,wu1,wu2,wu3,wl0,wl1,wl2,wl3,Re,cd0,k,cl_max,cl_min,cm0` (8 shape + Re + 5 coeffs).

If `polar_coeffs.csv` is absent, the runtime falls back to the analytic polar with a
warning, so the optimizer always runs.

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
| `thrust_z_offset` | 0.020 | thrust line above CG (m) |
| `hinge_moment_max` | 1.2 | fatal servo torque (kg·cm) |
| `sweep_crossflow_deg` | 25.0 | LE-sweep crossflow threshold (°) |
| `crossflow_factor` | 1.15 | outer-span profile-drag multiplier |
| `ncrit` | 4.0 | surrogate transition Ncrit (printed roughness) |
| `surrogate_mode` | `synthetic` | `synthetic` or `xfoil` |
| `xfoil_exe` | `tools/bin/xfoil.exe` | standalone XFOIL path (xfoil mode) |
| `avl_exe` | `tools/bin/avl.exe` | standalone AVL path (cross-check; manual for now) |
| `surrogate_samples` | 300 | DoE shapes sampled around the seeds |
| `surrogate_re` | `100000,150000,200000,300000` | Reynolds grid |
| `ga_pop` | 120 | population size |
| `ga_generations` | 150 | generations |
| `ga_seed` | 1 | RNG seed (reproducibility) |
| `ga_eta_cx` / `ga_eta_mut` | 15 / 20 | SBX / polynomial-mutation distribution indices |
| `ga_pcx` | 0.9 | crossover probability |
| `ga_pmut` | −1 | mutation prob (<0 ⇒ 1/n_genes) |
| `n_stations` | 20 | spanwise stations |

---

## Validation & tests

`build.ps1 -Test` builds and runs **31 dependency-free gates** (`tests/`, harness in
`test_harness.h`). Grouped:

- **Geometry** — CST round-trip, symmetric→zero camber, cambered→negative `α_L0`,
  reflex→`cm_ac` shift, rectangular planform area/AR, cosine lofting.
- **Mass** — solid-infill volume vs closed form, structure-only CG at ~0.42c, battery
  shift moves CG, spar clearance bounded.
- **Aero (Milestone-3 gates)** — Prandtl lift slope `2π·AR/(AR+2)`, Oswald range,
  induced-drag consistency, lift increases with α, **trim converges** to `CL_req`/`Cm=0`,
  SM-objective band, surrogate hull clamp.
- **Surrogate** — table load + IDW query + hull/cl clamps, analytic fallback when absent.
- **Airfoil I/O** — NACA 0012 max-thickness ≈ 12%, NACA 2412 cambered, NACA→CST round-trip,
  Selig + Lednicer parse.
- **Seeds** — NACA seeds load, bounds widen to contain seeds, seed genome decodes back to
  the seed shape, empty-seed graceful path.
- **NSGA-II** — constraint-domination rules, known-front non-dominated sort, crowding
  boundaries = ∞.

**AVL cross-check:** put `avl.exe` in `tools\bin\` (see `tools/bin/README.txt`), open
`out\knee.avl` in it, run `OPER` → `x`, and compare `CLa`, span efficiency, and `Cm`
against the panel solver. This is the standing oracle for the Milestone-3 panel solver
(`scratch/knee_xcheck.cpp` reconstructs all three decks and prints the panel/VLM rows).

**Reproducibility / OpenMP determinism** are both verified (same seed ⇒ identical
`pareto.csv`, serial and parallel).

---

## Known limitations / approximations

These are deliberate placeholders behind clean interfaces, slated for replacement:

- **3D aero is the Morino panel method** (default); the analytic VLM reference model is the
  `aero_model = vlm` fallback.
- **High-α NP migration** is a sweep/α heuristic, not a real high-α solve.
- **Control derivatives** (elevon/elevator effectiveness, hinge moment) use simple flap
  theory — adequate for ranking, not absolute servo sizing.
- **Split-mode roll authority** and the 4:1 differential mixing aren't modeled (only pitch
  trim is solved); the control mode is a gene but only affects the pitch surface today.
- **Hardware keep-outs** (motor zone, avionics 20–50% span placement) are assumed, not
  enforced as hard constraints.
- **Synthetic surrogate** is the default; absolute drag numbers firm up once real XFOIL
  tables are generated.

---

## Roadmap

### Milestone 3 — Morino Dirichlet panel solver *(done — default model)*
`aero_model = panel` (`src/aero_panel.cpp`) is the default 3D constant-strength
doublet/source panel method, dispatched from `aero_potential::solve()` behind the
unchanged signature.
- Dense AIC built per geometry and **direct-LU** factored once (Eigen), reusing the
  factorization across the trim Newton's RHS. The wake is **frozen across the Newton solve**
  (`panel_trim_freeze_wake`, default on): the wake inclination is irrelevant to attached
  loading (verified: dCL ~1e-11), so rebuilding it per alpha step was pure cost — freezing
  it gives ~13× faster trims (one AIC build per candidate).
- Trailing-edge wake via doublet-strength difference, freestream-aligned then frozen.
- Strip viscous coupling and the `solve()` signature unchanged.
- **Validated:** passes every `tests/test_aero.cpp` analytic gate (panel-mode lift slope,
  elliptic `e≈1`, trim, quarter-chord NP) **and** matches AVL within a few percent on the
  knee/min_drag/min_mass decks. The VLM remains the documented fallback (`aero_model = vlm`).

### Milestone 4 — Real viscous tables
Run `build_surrogate.exe` with `surrogate_mode = xfoil` over the seed + DoE shapes; widen
the Reynolds grid; verify NACA 0012 @ Re 2e5 gives `cd0 ≈ 0.012–0.02`. Add a fallback
log/report of XFOIL convergence rate.

### Milestone 5 — High-α & nonlinear stability
Replace the NP-migration heuristic with a high-α evaluation from the panel solver; add a
proper post-stall `cl` cap per station; re-validate the tip-stall watchdog.

### Milestone 6 — Full control & hardware constraints
Model split-mode roll authority and the 4:1 differential mixing; size servos from real
hinge moments; add motor keep-out and avionics-placement hard constraints; optionally make
control-surface chord/span design variables.

### Milestone 7 — UX & analysis
- Pareto-front visualization (export to a plotting-friendly format / small viewer).
- Resume/restart from a saved population; checkpointing.
- Sensitivity / trade-study mode (sweep one variable, hold the rest).
- Optional STL/STEP export of the knee design for printing.

### Backlog / ideas
- Eigen-accelerated AIC solve once panel counts grow.
- Multi-point optimization (climb + cruise) instead of single cruise point.
- Manufacturing constraints (min wall, overhang) folded into the mass model.
- Battery/CG as an inner solve targeting mid-band SM instead of a free gene.

---

## Equations reference

Quick index of the core relations (see the cited source for context):

| Quantity | Expression | Source |
|----------|------------|--------|
| CST surface | `y = x^N1(1−x)^N2 · Σ wᵢ Bᵢ(x) ± x·te/2` | `geom.cpp` |
| 3D lift slope | `a = a0 / (1 + a0/(π e AR))` | `aero_potential.cpp` |
| Induced drag | `CDi = CL² / (π e AR)` | `aero_potential.cpp` |
| Section drag polar | `cd = cd0 + k·cl²` (clamped to `cl_min..cl_max`) | `aero_viscous.cpp` |
| Sweep-normal query | `Re = ρ V cosΛ c / μ`, `cl_n = cl/cos²Λ` | `aero_potential.cpp` |
| Static margin | `SM = (x_np − x_cg) / MAC` | `aero_potential.cpp` |
| Pitch moment | `Cm = cm_ac + CL(x_cg−x_np)/MAC + Cm_δ δ − Cm_thrust` | `aero_potential.cpp` |
| Trim | solve `CL = W/(qS)`, `Cm = 0` for `(α, δ_e)` (Newton) | `stability.cpp` |
| Required CL | `CL_req = W / (½ ρ V² S)` | `stability.cpp` |

---

*AeroAnalyzer Pro — Milestone 3 (Morino panel solver is the default model of record).
Contributions/edits welcome; keep the `tests/test_aero.cpp` gates green and `build.ps1`
ASCII-only.*
