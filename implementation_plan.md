# P5A + P5B + P5C + P6 — Final Implementation Plan

**Status: READY FOR EXECUTION — all design questions answered.**

---

## Summary of Design Decisions

| Question | Answer |
|----------|--------|
| Number of sections K | **5** |
| η breakpoints | **{0.0, 0.5, 0.75, 0.875, 1.0}** — 3 sections in the outer 25% for tip formation |
| CMA-ES role | **Supplementary offspring generator (helper)** — NSGA-II loop unchanged |
| Chord-collapse cv gate | **Yes** — penalize any station chord below a configurable floor |
| Checkpoint CMA state | **Yes** — full C matrix + σ + evolution paths serialized |

---

## Section Layout Rationale

5 sections, η = `{0.0, 0.5, 0.75, 0.875, 1.0}`:

```
Root ──────────────── 50% ─────── 75% ── 87.5% ── Tip
  s0        s1                    s2      s3       s4
  |<────────── inner 75% ─────────>|<── outer 25% ──>|
         (2 sections, 1 segment)      (3 sections,
                                       2 segments)
```

The inner 75% gets 2 sections (1 blend segment) — wide coverage with few genes.
The outer 25% gets 3 sections (2 blend segments of Δη=0.125 each) — dense sampling
of the tip, where chord is small, Re is low, reflex is critical, and roll control lives.

The tip cluster gives the optimizer independent control over:
- The main tip airfoil shape (s3 at η=0.875, the structural tip chord midpoint)
- The geometric wingtip profile (s4 at η=1.0, can close to a thin wafer)
- The transition from panel airfoil to tip (s2→s3, captures the "carve" into a winglet-like form)

**Genome: 52 genes total.**

| Group | Genes | Count |
|-------|-------|------:|
| Planform | `root_chord`, `taper`, `semi_span`, `le_sweep`, `washout`, `battery_x` | 6 |
| Control surface | `te_frac`, `mode`, `cs_chord_frac`, `ail_span_frac` | 4 |
| Organic planform | `le_bow`, `te_bow` | 2 |
| s0 root CST | `s0_wu0..3`, `s0_wl0..3` | 8 |
| s1 (η=0.50) CST | `s1_wu0..3`, `s1_wl0..3` | 8 |
| s2 (η=0.75) CST | `s2_wu0..3`, `s2_wl0..3` | 8 |
| s3 (η=0.875) CST | `s3_wu0..3`, `s3_wl0..3` | 8 |
| s4 tip CST | `s4_wu0..3`, `s4_wl0..3` | 8 |
| **Total** | | **52** |

---

## P5A — Generalized 5-Section Spanwise Loft

### Key algorithm: piecewise loft with non-uniform η

```cpp
// η breakpoints (config key profile_etas, default shown)
const double ETA[5] = {0.0, 0.5, 0.75, 0.875, 1.0};
const int    K = 5;

// For each lofted station at spanwise fraction t in [0,1]:
// 1. Find segment: largest k such that ETA[k] <= t
int seg = K - 2;
for (int k = 0; k < K-1; ++k) if (ETA[k+1] >= t) { seg = k; break; }
// 2. Local blend parameter within segment
double t_local = (ETA[seg+1] > ETA[seg])
               ? (t - ETA[seg]) / (ETA[seg+1] - ETA[seg])
               : 0.0;
// 3. Linear blend of CST weights
s.af.wu[j] = (1.0-t_local)*sections[seg].wu[j] + t_local*sections[seg+1].wu[j];
```

The cosine clustering of stations (existing `y[i] = semi_span * 0.5*(1-cos(π*i/(n-1)))`)
naturally clusters stations near the tip — where our 3 control sections live. The tip cluster
will have many stations covering short η intervals, so the 3-section shape will be faithfully
reproduced by the lofted mesh.

### Gene enum refactor

The `Gene` enum in `geom.h` is refactored from hardcoded root/tip constants to a
generated layout. Two options:

**Option A — Keep enum, expand it:**
```cpp
enum Gene {
    G_ROOT=0, G_TAPER, G_SEMISPAN, G_SWEEP, G_WASHOUT, G_BATTERY,
    G_TE, G_MODE, G_CS_CHORD, G_AIL_SPAN,
    G_LE_BOW, G_TE_BOW,
    G_S0_WU0, G_S0_WU1, G_S0_WU2, G_S0_WU3,
    G_S0_WL0, G_S0_WL1, G_S0_WL2, G_S0_WL3,
    G_S1_WU0, ..., G_S1_WL3,    // 8 genes
    G_S2_WU0, ..., G_S2_WL3,    // 8 genes
    G_S3_WU0, ..., G_S3_WL3,    // 8 genes
    G_S4_WU0, ..., G_S4_WL3,    // 8 genes
    N_GENES  // = 52
};
```

**Option B — Computed offset function:**
```cpp
constexpr int N_PLANFORM = 12;        // 6 planform + 4 ctrl + 2 bow
constexpr int N_CST_PER_SEC = 8;      // 4 wu + 4 wl
constexpr int N_SECTIONS = 5;
constexpr int N_GENES = N_PLANFORM + N_SECTIONS * N_CST_PER_SEC;  // = 52

inline int G_SEC(int sec, int wu_or_wl, int i) {
    // sec 0..4, wu_or_wl: 0=upper, 1=lower, i 0..3
    return N_PLANFORM + sec * N_CST_PER_SEC + wu_or_wl * 4 + i;
}
```

**Recommendation: Option B.** It's K-agnostic — changing `N_SECTIONS` from 5 to 6 in the
config automatically extends the genome without touching the enum. `seeds.cpp` and
`evaluate.cpp` use `N_PLANFORM + k * N_CST_PER_SEC` offsets instead of `G_TIP_WU0` etc.

### Blast radius

| File | Change |
|------|--------|
| **[MODIFY] [geom.h](file:///c:/Users/kadan/Downloads/aeroanalyzer/include/aeroanalyzer/geom.h)** | Replace `Gene` enum with Option B computed offsets. Add `N_SECTIONS=5`, `N_PLANFORM=12`, `G_LE_BOW`, `G_TE_BOW` constants. `N_GENES = 52`. |
| **[MODIFY] [engine_core.h](file:///c:/Users/kadan/Downloads/aeroanalyzer/include/aeroanalyzer/engine_core.h)** | `WingGeometry`: replace `Airfoil section; Airfoil section_tip;` with `std::vector<Airfoil> sections; // size K`. Add `double le_bow=0, te_bow=0;`. Keep `mid_eta` unused (organic bow takes its place). |
| **[MODIFY] [geom.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/geom.cpp)** | `default_genome()`: loop over K sections, set CST bounds identically for all. Add `le_bow / te_bow` bounds ±0.05 m. `decode()`: populate `w.sections[0..K-1]`, `w.le_bow`, `w.te_bow`. `loft()`: piecewise blend as shown above. |
| **[MODIFY] [seeds.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/seeds.cpp)** | `widen_cst_bounds()`: loop over all K sections. `build_seed_genomes()`: elites copy seed airfoil to ALL K sections; hybrids jitter each independently. |
| **[MODIFY] [evaluate.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/evaluate.cpp)** | Thickness gate: loop over `w.sections[0..K-1]` (not just root+tip). Add chord-collapse gate (see P5B). |
| **[MODIFY] [avl_export.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/avl_export.cpp)** | Write K SECTION entries at η positions, each with its own AFILE. |
| **[MODIFY] [aero_panel.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/aero_panel.cpp)** | `geom_sig_panel()` lines 317-327: already hashes `w.stations` (all lofted sections). No change needed. |
| **[MODIFY] [tests/test_geom.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/tests/test_geom.cpp)** | Add: K=5 decode→loft round-trip; station at η=0.875 gets blend of s3/s4; station at η=0.50 is exactly s1. |
| **[MODIFY] [tests/test_seeds.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/tests/test_seeds.cpp)** | Assert `spec.size() == 52`. |

---

## P5B — Organic Wing (Curved Planform)

### Bow formula

```
x_le(y) = y * tan(Λ) + le_bow * 4 * (y/b) * (1 - y/b)
```
where `b = semi_span`. Applied in `loft()` per station. `te_bow` similarly shifts the TE:
the TE x-position at each station is `x_le + chord`, so `te_bow` is an **additional chord
offset** — it bows the TE relative to the LE:
```
x_te(y) = x_le(y) + chord(y) + te_bow * 4 * (y/b) * (1 - y/b)
```
This is equivalent to a spanwise-varying chord augmentation, which is aerodynamically
meaningful (the optimizer can use it to add effective sweep on the TE without changing the
LE angle).

### Chord-collapse cv gate (Q3 answer)

During `evaluate()`, after lofting:
```cpp
// (14) Chord-collapse gate: no lofted station may have chord below floor.
// Catches le_bow/te_bow combinations that push TE forward of LE.
double chord_min_m = cfg.getd("chord_min_m", 0.03);   // 30 mm floor
for (const auto& s : r.geom.stations)
    if (s.chord < chord_min_m) cv += 50.0 * (chord_min_m - s.chord) / chord_min_m;
```

The `te_bow` augments effective chord but never independently controls x_le — so the LE
never moves backward. The risk is instead `te_bow < 0` making effective chord negative at
the tip. The gate above catches this for any station.

### Impact on downstream modules

The Morino panel `build_mesh()` reads `s.x_le` and `s.chord` from lofted stations — both
now curved. No changes needed to the panel solver. AVL SECTION entries write `x_le` per
section, which also automatically picks up the bow.

---

## P5C — CMA-ES as Supplementary Offspring Helper

### Architecture: "CMA bonus batch"

The NSGA-II loop is **not modified**. CMA-ES is an additional offspring source that runs
alongside SBX+poly-mut, not instead of it.

**Every generation:**
1. NSGA-II generates offspring `Q` via the unchanged SBX+poly-mut operators (size = `ga_pop`).
2. CMA generates a bonus batch `Q_cma` of `n_cma` candidates from the learned distribution
   (default: `n_cma = max(10, ga_pop/5)` = 24 for pop=120).
3. The elitist merge becomes: `R = P ∪ Q ∪ Q_cma` (size = 2×pop + n_cma).
4. Standard NSGA-II fast_nondominated_sort + crowding truncation reduces R back to `ga_pop`.
5. CMA state updates from the **front-0 feasible** candidates of the current P.

The CMA candidates simply compete on equal footing with SBX offspring. No special treatment —
if a CMA candidate is dominated, it's discarded at truncation. The overhead is evaluating
`n_cma` extra candidates per generation: at pop=120, that's 24 extra evaluations out of 240
total offspring (10% cost increase). With 12 threads this is absorbed in the same wall-clock
slot.

### CMA state

```cpp
struct CMAState {
    int    n = 0;            // genome dimension (52)
    int    mu = 0;           // number of selected individuals for update
    double sigma = 0.3;      // global step size (initialized to 0.3 * avg gene range)
    Eigen::VectorXd  mean;   // distribution mean
    Eigen::MatrixXd  C;      // covariance matrix (n×n); initialized to I
    Eigen::VectorXd  p_c;    // evolution path for rank-1 C update
    Eigen::VectorXd  p_sigma;// evolution path for step size (CSA)
    Eigen::MatrixXd  L;      // Cholesky factor of C; recomputed each gen
    double mu_eff  = 0.0;    // variance effective selection mass
    // CMA constants (set once at initialization based on n, mu):
    double c_sigma, d_sigma, c_c, c_1, c_mu;
    bool initialized = false;
};
```

**Update rule (called with front-0 feasible candidates, sorted by crowding then drag):**

```cpp
void cma_update(CMAState& cma, const std::vector<Candidate>& selected,
                const geom::GenomeSpec& spec) {
    // Normalize genes to [0,1] for the covariance (avoids bound scale differences)
    // ... compute weighted mean, update p_sigma, p_c, C, sigma ...
    // Recompute L = chol(C) for sampling next generation
}
```

**Sampling:**

```cpp
Candidate cma_sample(CMAState& cma, const geom::GenomeSpec& spec, Rng& rng) {
    std::normal_distribution<double> N(0.0, 1.0);
    Eigen::VectorXd z(cma.n);
    for (int i = 0; i < cma.n; ++i) z(i) = N(rng.gen);
    Eigen::VectorXd x_norm = cma.mean + cma.sigma * cma.L * z;
    // Map back from [0,1] to gene bounds and clamp
    Candidate c;
    for (int i = 0; i < cma.n; ++i)
        c.genes[i] = std::min(spec.hi[i], std::max(spec.lo[i],
                     spec.lo[i] + x_norm(i) * (spec.hi[i] - spec.lo[i])));
    return c;
}
```

**CMA constants** (Hansen defaults for n=52):
- `mu = n_cma / 2 = 12` (half the CMA batch size is used for the update)
- `mu_eff = 1 / Σwi²` where `wi = ln(mu+0.5) - ln(i)` (normalized)
- `c_sigma ≈ (mu_eff+2)/(n+mu_eff+5)` ≈ 0.065
- `d_sigma ≈ 1 + c_sigma + 2*max(0, sqrt((mu_eff-1)/(n+1))-1)` ≈ 1.10
- `c_c ≈ (4+mu_eff/n)/(n+4+2*mu_eff/n)` ≈ 0.092
- `c_1 ≈ 2/(n+1.3)^2` ≈ 0.00073
- `c_mu ≈ 2*(mu_eff-2+1/mu_eff)/((n+2)^2+mu_eff)` ≈ 0.0035

The Cholesky recomputation at n=52 costs ~52³/3 ≈ 47K flops (~50 ns). Negligible.

> [!IMPORTANT]
> **Determinism is preserved.** CMA uses the same `Rng` (seeded `mt19937`) as SBX. The CMA
> update is deterministic given the population. The evaluation of `Q_cma` is included in the
> same OpenMP parallel loop as Q — the `thread_local` panel cache already handles this. The
> ordering of CMA candidates within the merged R is deterministic (appended after Q).

> [!NOTE]
> **Cold-start behavior:** For gen 0 and the first few generations before enough feasible
> front-0 candidates exist (< mu feasible), CMA falls back to sampling from `N(mean, σ²I)`
> (identity covariance) until at least `mu` feasible front-0 candidates are available. This
> means CMA is effectively a wide random sampler early on and adapts into a correlated sampler
> as the front forms.

### Blast radius (P5C)

| File | Change |
|------|--------|
| **[MODIFY] [ga.h](file:///c:/Users/kadan/Downloads/aeroanalyzer/include/aeroanalyzer/ga.h)** | Add `CMAState` struct (forward-decl or full). Add `int n_cma_helpers = -1;` (−1 → pop/5) to `Params`. |
| **[MODIFY] [ga.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/ga.cpp)** | Add `CMAState` full definition. Add `cma_update()` and `cma_sample()` free functions. In `NSGA2::run()`: initialize CMA at gen 0. After generating `Q`, generate `Q_cma` batch. Evaluate `Q_cma` in the same parallel `eval_pop` call (extend the loop range). Merge `R = P ∪ Q ∪ Q_cma`. |
| **[MODIFY] [tests/test_ga.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/tests/test_ga.cpp)** | Add: CMA update is deterministic; sampling produces vectors within gene bounds; CMA batch candidates are present in R after merge. |

---

## P6 — UX & Analysis Features

### P6-A: Checkpoint / Restart (with CMA state)

**File format** (plaintext, human-readable):

```
# aeroanalyzer checkpoint v2
# genome_hash = <fnv1a hex>
# generation = <N>
# n_genes = 52
# pop_size = 120
# cma_sigma = 0.187234
# cma_mean = 0.23 0.45 ... (52 values)
# cma_p_sigma = ... (52 values)
# cma_p_c = ... (52 values)
# cma_C = ... (52*52 = 2704 values, row-major, one row per line with prefix "# cma_C_row_<i> = ")
# ---
genes[0],genes[1],...,genes[51],cv,rank
...one row per candidate...
```

The CMA C matrix is 52×52 = 2704 doubles = ~21 KB — negligible. Writing it once per
generation (to a temp file then atomic rename) adds < 1 ms of I/O.

**Genome hash:** FNV-1a of all `spec.names` concatenated. Detects genome changes (e.g., if
K changes from 5 to 6 between runs). A hash mismatch prints a warning and starts fresh.

#### [MODIFY] [ga.h](file:///c:/Users/kadan/Downloads/aeroanalyzer/include/aeroanalyzer/ga.h)
- Add `std::string checkpoint_path` to `Params`.

#### [MODIFY] [ga.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/ga.cpp)
- `save_checkpoint()`: write to `<path>.tmp`, then `rename()` — atomic on Windows (POSIX
  rename-over is atomic; on Windows use `MoveFileExW` with `MOVEFILE_REPLACE_EXISTING`).
- `load_checkpoint()`: parse header, verify hash, restore P, restore CMA state.

#### [MODIFY] [main.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/app/main.cpp)
- `checkpoint_path = cfg.gets("checkpoint_path", "out/checkpoint.csv")`.

---

### P6-B: Pareto Visualization

#### [NEW] [tools/plot_pareto.py](file:///c:/Users/kadan/Downloads/aeroanalyzer/tools/plot_pareto.py)

Reads `out/pareto.csv`. Produces four figures:

1. **`pareto.png`** — 2×2 subplot: 3D scatter (drag/mass/SM, colour = static margin), plus
   three 2D projections.
2. **`pareto_genes.png`** — gene heatmap: rows = Pareto designs, columns = 52 genes,
   normalized [0,1]. Annotated with gene names from CSV header. Colour bands reveal which
   genes are converged (narrow) vs diverse (wide spectrum).
3. **`profiles.png`** — For the 3 incumbents: plots all 5 section airfoil shapes (upper +
   lower surface) overlaid, colour-coded by η. Plus a planform outline showing the organic
   LE/TE bow (if non-zero).
4. **`cma_evolution.png`** (optional, if checkpoint CSV is present) — σ vs generation, showing
   CMA step size adaptation. Useful for diagnosing premature convergence.

---

### P6-C: Structured Export (Per-Incumbent Folders + Per-Section .dat)

```
out/
  min_drag/
    wing.avl
    section_s0_root.dat
    section_s1_eta50.dat
    section_s2_eta75.dat
    section_s3_eta875.dat
    section_s4_tip.dat
    panel_ref.txt
  min_mass/ ...
  knee/ ...
  checkpoint.csv
  pareto.csv
```

AVL SECTION entries use the η-interpolated geometry (`x_le` from the curved planform):

```
SECTION
<x_le_at_eta> <y_at_eta> 0.0 <chord_at_eta> <twist_at_eta_deg>
AFILE
section_s2_eta75.dat
```

#### [MODIFY] [avl_export.h](file:///c:/Users/kadan/Downloads/aeroanalyzer/include/aeroanalyzer/avl_export.h)
- Change: `bool write_case(const std::string& out_dir, ...)`.

#### [MODIFY] [avl_export.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/src/avl_export.cpp)
- `create_directories(out_dir)`.
- Write `section_sk_<label>.dat` for each of the K control sections.
- `wing.avl`: K SECTION blocks with relative AFILE paths and curved planform x_le values.
- `panel_ref.txt`: existing sidecar content (unchanged internally).

#### [MODIFY] [main.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/app/main.cpp)
- `avl::write_case("out/" + name + "/", ...)`.

#### [MODIFY] [validate_avl.ps1](file:///c:/Users/kadan/Downloads/aeroanalyzer/validate_avl.ps1)
- Path: `out/<name>/wing.avl`, `out/<name>/panel_ref.txt`.

---

### P6-D: Sensitivity Sweeps

#### [MODIFY] [main.cpp](file:///c:/Users/kadan/Downloads/aeroanalyzer/app/main.cpp)

```powershell
# Usage:
aeroanalyzer.exe config/baseline.cfg --sweep knee le_bow 30
# → varies le_bow across its bounds at 30 points, knee as base
# Output: out/sweep_le_bow.csv
```

Parses `--sweep <incumbent_name> <gene_name> [N=30]`:
1. Reads incumbent gene vector from `out/<name>/panel_ref.txt` (add gene values to sidecar
   in P6-C).
2. Finds gene index by matching `gene_name` against `spec.names`.
3. Evaluates N points via `eval.run()` (single-threaded sweep, fast).
4. Writes `out/sweep_<gene_name>.csv`: columns `gene_val, drag_N, mass_kg, sm_dev, sm_pct, cv`.

#### [NEW] [tools/plot_sweep.py](file:///c:/Users/kadan/Downloads/aeroanalyzer/tools/plot_sweep.py)
- Reads `out/sweep_<gene>.csv`.
- 3-panel plot: drag / mass / SM vs gene value; incumbent value marked with a vertical line.
- Saves `out/sweep_<gene>.png`.

---

## Execution Order (Dependency-Ordered)

| Step | Scope | Blocking dependency |
|------|-------|---------------------|
| 1 | **P6-C** export folder restructure | — (first so avl_export API changes once) |
| 2 | **P6-A** checkpoint/restart (no CMA yet) | P6-C (`out/` layout fixed) |
| 3 | **P5A** K=5 genome + loft | P6-A (so checkpoint works with new genome) |
| 4 | **P5B** organic wing bow genes | P5A (adds genes to existing P5A genome block) |
| 5 | **P5C** CMA helper | P5A (needs genome dimension n=52 fixed) |
| 6 | **P6-A v2** CMA state in checkpoint | P5C (CMAState struct must exist) |
| 7 | **P6-B** visualization | P5A (uses K-section gene columns from pareto.csv) |
| 8 | **P6-D** sweep | P6-C (uses `panel_ref.txt` gene storage) |

Build and test after **every step**:
```powershell
powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test
```

---

## Full File Change Index

| File | Steps that touch it |
|------|-------------------|
| `include/aeroanalyzer/geom.h` | P5A, P5B |
| `include/aeroanalyzer/engine_core.h` | P5A, P5B |
| `include/aeroanalyzer/ga.h` | P5C, P6-A |
| `include/aeroanalyzer/avl_export.h` | P6-C |
| `src/geom.cpp` | P5A, P5B |
| `src/seeds.cpp` | P5A |
| `src/evaluate.cpp` | P5A (thickness loop), P5B (chord-collapse gate) |
| `src/avl_export.cpp` | P6-C (folder + K sections) |
| `src/ga.cpp` | P5C (CMA), P6-A (checkpoint) |
| `app/main.cpp` | P6-C (folder call), P6-A (checkpoint cfg), P6-D (sweep mode) |
| `validate_avl.ps1` | P6-C (path update) |
| `tests/test_geom.cpp` | P5A, P5B |
| `tests/test_seeds.cpp` | P5A |
| `tests/test_ga.cpp` | P5C |
| `tools/plot_pareto.py` | P6-B (new file) |
| `tools/plot_sweep.py` | P6-D (new file) |

---

## Verification Plan

### Automated (must stay green at every step)

```powershell
powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test
```

- **test_geom.cpp** additions:
  - K=5 decode: `w.sections.size() == 5`
  - Station at η=0.875 blends s3↔s4, not s2↔s3
  - Station at η=0.75 exactly equals s2 (at the knot point)
  - `le_bow = 0.03`: midspan x_le = 0.5*tan(sweep)*semi_span + 0.03 (within 1e-9)
  - `le_bow = 0`, `te_bow < 0`: chord at midspan < nominal; triggers chord-collapse gate
- **test_seeds.cpp:** `spec.size() == 52`
- **test_ga.cpp:** CMA update deterministic; sampled genes within bounds; R.size() = 2*pop + n_cma after merge

### Calibration runs (manual)

| Run | Config | Success criterion |
|-----|--------|-------------------|
| Baseline K=5, 50 gen, no CMA | `n_cma_helpers=0` | ≥ 10 feasible in front-0; no geometry crashes |
| K=5 + CMA helper, 50 gen | `n_cma_helpers=24` | ≥ 10 feasible; CMA σ decreases from gen 0→50 (converging, not exploding) |
| Checkpoint test | Kill at gen 40, restart with same cfg | Resumes at gen 40; gen 41 population identical to uninterrupted run |
| AVL cross-check | `validate_avl.ps1` (new paths) | PASS: e ≤ 3.5%, Xnp ≤ 3 mm for all 3 incumbents |
| Organic wing | `le_bow` in [−0.05,0.05]; 50 gen | No NaN in panel solver; `out/min_drag/wing.avl` shows curved LE visually in AVL |

### Manual inspection

- `out/knee/wing.avl`: 5 SECTION entries with η-correct x_le and chord values
- `tools/plot_pareto.py` → `out/profiles.png`: 5 overlaid airfoil shapes per incumbent, tip cluster showing differentiated s2/s3/s4 shapes
- CMA σ trajectory in `out/pareto_genes.png` or a log print each 10 gens showing σ adapting
