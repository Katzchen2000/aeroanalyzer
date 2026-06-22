Paste this whole file as the opening message of a fresh AeroAnalyzer session.
Deep background is in the milestone3-status memory -- read it first (the "SESSION 2"
block at the end is the current state of record).

---
You are continuing AeroAnalyzer Milestone 3: making the Morino panel solver the default and
demoting the VLM to a documented fallback. The previous session finished the correctness work
(steps 1-3) and staged the cross-check harness; this session does the robustness + perf gates
(steps 5, 6), supports the AVL pass (step 4, I drive AVL), and THEN flips the default + docs
(step 7). Do NOT flip until 5+6 pass and I have run AVL.

## State of record (already done last session -- verify, don't redo)

- **Step 1 (span efficiency e) DONE.** `panel::span_efficiency` (src/aero_panel.cpp ~508): the
  bogus 0.30 physical floor is replaced by `E_DEGEN=0.05` (numerical guard only). The fit is
  sound; the old "e=0.30 vs AVL 0.41" was a diagnostic-alpha artifact (basic/washout loading at
  low alpha). At trim, e is healthy. Probe: scratch/e_probe.cpp.
- **Step 2 (panel-mode analytic gates) DONE.** tests/test_aero.cpp gained panel-mode gates
  (lift slope vs Prandtl, elliptic e~1, induced-drag consistency, trim converges, x_np at
  quarter chord). **All 42 tests pass** (`build_mingw.ps1 -Test`).
- **Step 3 (x_np) DONE -- premise overturned.** The quarter-chord proxy in solve() is ALREADY
  within ~1.5-2% MAC of AVL (rect unswept -> exactly 0.25c). A chordwise load-centre method was
  tried and REJECTED (worse on reflexed sections). solve() unchanged; `neutral_point_load` /
  `panel_xnp_debug` kept as diagnostics only. Probe: scratch/xnp_probe.cpp.
- **Cross-check harness STAGED.** scratch/knee_xcheck.cpp now PARSES the current out/*.avl decks
  (knee/min_drag/min_mass), reconstructs each wing, runs panel+VLM, and prints CLa/e/x_np/Cm/SM
  with blank AVL rows. Panel and VLM agree (CLa <3%, x_np <0.4% MAC, e a few %).

## Build / run (MinGW + Eigen is primary; keep build*.ps1 pure ASCII)

    powershell -ExecutionPolicy Bypass -File build_mingw.ps1 -Test     # build + 42 gates
Scratch probes: compile src/*.cpp + the probe with the SAME flags build_mingw.ps1 uses:
    g++ -std=c++17 -O2 -fno-math-errno -fopenmp -DEIGEN_DONT_PARALLELIZE
        -static-libgcc -static-libstdc++ -I include -I "<eigen>" src/*.cpp scratch/x.cpp -o build/x.exe
with C:\msys64\mingw64\bin on PATH. **Caveats:** (1) Eigen 5.0.0 FAILS TO LINK at -O0 here --
use -O2 for anything that links Eigen. (2) Do NOT run the optimizer (`build_mingw.ps1 -Run`)
casually: it overwrites out/*.avl/*.dat, which the AVL cross-check (step 4) depends on. (3) The
OLD knee airfoil is gone (out/ was regenerated); the baked AVL oracle e=0.4148/Xnp=0.07306 was
for a sweep-8 planform that no longer matches the current decks (sweep 24-30).

## THE perf crux (read before steps 5/6)

Every panel `solve()` rebuilds a dense ~N*N AIC from transcendental kernels (~10ms; N=2*nc*
strips). `stability::trim`'s finite-difference Jacobian rebuilds AND refactorizes on EVERY
alpha perturbation, because the panel cache key includes `wake_alpha` -- a ~50-100x cost
multiplier. This is why naive fuzz/timing harnesses are slow, and it is the single biggest
lever for the GA hot path.

## Step 5 -- Robustness fuzz (do first; pure code). scratch/panel_fuzz.cpp exists but is too slow.

The panel must never NaN/Inf/diverge on any in-bounds genome (one poisoned candidate corrupts
the Pareto sort). The current fuzz calls solve thousands of times at production resolution, so
it is too slow. Redesign for speed -- RECOMMENDED: **targeted corner cases at production
resolution.** Enumerate combinations of the ~6 risky genes (G_SWEEP, G_TAPER, G_WASHOUT, G_TE,
G_SEMISPAN, plus extreme CST reflex) at their lo/hi bounds + random fill for the rest -> a few
hundred cases, real nc=10, ~2-3s. Run decode->loft->massprops->panel::solve at a few
(alpha,delta_e) and assert every output finite/bounded; run trim() on ~20 of them. (Alt: coarse
mesh nc=3/n_stations=6, thousands of random genomes, also ~3s -- finiteness is resolution-
independent.) **Accept: zero non-finite/divergent results; fix any guard/division that trips.**

## Step 6 -- Perf + determinism. scratch/panel_timing.cpp exists but runs a GA (too slow); redesign.

1. **Microbenchmark, do NOT run a GA for timing.** Time ONE panel solve, ONE vlm solve, ONE
   panel trim; extrapolate to pop120 x gen150 arithmetically. Report the panel/vlm cost ratio
   and the projected wall-clock.
2. **Fix the perf crux (biggest win).** Make trim's FD Jacobian reuse a single factorization --
   freeze the wake across the alpha perturbations (decouple wake_alpha from the FD step, or add
   a frozen-wake mode the Jacobian uses). Re-validate trim_converges + all gates afterwards.
   Recommend a default panel_chordwise if the production cost needs it (do NOT change the
   default config yet).
3. **Determinism.** Confirm same-seed result is identical serial vs OpenMP (tiny GA, e.g.
   pop 12 gen 4, is enough; the panel cache is thread_local keyed on geom_sig_panel).
   **Accept: acceptable wall-clock + byte-identical reproducibility.**

## Step 4 -- AVL cross-check (I drive AVL; you prep + interpret)

Decks are in out/{knee,min_drag,min_mass}.avl. Run tools/bin/avl352.exe headless on each (pipe
a .cmd to stdin, OPER -> x, parse the ST dump; pattern in ~/Downloads/testfitness/sweep_avl.cpp),
read CLa/e/Xnp/Cm. Fill the blank AVL rows printed by `build/knee_xcheck.exe` and compare.
**Accept: panel CLa/e/x_np/Cm within a few % of AVL on all THREE decks.**

## Step 7 -- Flip the default + reconcile docs (DEFERRED until 5+6 pass AND AVL is clean)

Change the dispatch default to panel (src/aero_potential.cpp ~247 `cfg.gets("aero_model","vlm")`
and/or add `aero_model = panel` to config/baseline.cfg), regenerate out/, retarget the decks,
and demote the VLM to the documented fallback in: CLAUDE.md (milestone table + "analytic
reference model" language), app/main.cpp banner (~lines 42-44), and the milestone3-status
memory. Keep all 42+ gates green and build*.ps1 pure ASCII.

## Files touched last session
- src/aero_panel.cpp (e floor -> E_DEGEN guard; added neutral_point_load + panel_xnp_debug diagnostics)
- include/aeroanalyzer/aero_panel.h (XnpDebug struct + panel_xnp_debug decl)
- tests/test_aero.cpp (5 new panel-mode gates; 42 total)
- scratch/: e_probe.cpp, xnp_probe.cpp, knee_xcheck.cpp (rebuilt), panel_fuzz.cpp + panel_timing.cpp (need redesign)

## Order of attack
Step 5 (fast fuzz) -> step 6 (microbench timing + the trim wake-freeze perf fix + determinism)
-> step 4 (I run AVL on the 3 decks) -> only then step 7 (flip + docs). Verify with
`build_mingw.ps1 -Test` (all gates green) at each stage. Update the milestone3-status memory at the end.
