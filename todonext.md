# TODO next

Current priority list as of 2026-06-24. Build/test after each change with ./build_mingw.ps1 -Test.

## Release hygiene
- [x] Remove the per-evaluation stdout diagnostic in src/evaluate.cpp.
- [x] Keep panel_relaxed_wake off in the GA hot path.

## Safety gates
- [x] Re-enable spar/OML clearance enforcement.
- [x] Re-enable hinge-moment enforcement.
- [x] Decide whether hardware keep-out is a hard gate or a warning.

## Trust gate
- [x] Thread NeuralFoil confidence into the evaluator and gate low-confidence results.

## Validation
- [x] Automate the AVL cross-check for min_drag, min_mass, and knee.
      validate_avl.ps1 drives avl_exe over the three incumbents; out/<name>_panel.txt sidecar written by aeroanalyzer.exe.
      Re-run on improved-baseline incumbents (2026-06-25): e within 2.7%, Xnp within 2mm — PASS all cases.
      CL 21% delta INFO-only (expected trim-alpha mismatch, not gated).
- [x] Before/after GA comparison: sm_floor_penalty=0 vs =30.
      gen 10: 27/120 feasible (sm=0) vs 2/120 (sm=30). sm=30 also collapses planform: all 7
      Pareto designs had identical root_c/tip_c/span. Root cause: SM objective already enforces
      the band; penalty at weight 30 double-counts and kills diversity.
      Fix applied: sm_floor_penalty 30->5, ga_eta_cx 15->10, ga_eta_mut 20->10.

## Optimizer quality fixes (applied 2026-06-25)
- [x] SM floor penalty reduced from 30 to 5 (prevents planform collapse; gen10 feasibility
      improves 13x; SM still enforced as cv floor, just less dominant).
- [x] GA eta values lowered: eta_cx 15->10, eta_mut 20->10 (wider crossover/mutation
      distribution for more planform and CST diversity).
- [x] pareto.csv now includes all 26 raw gene values (diagnosis column for future convergence checks).
- [x] Thread tuning: OMP=12 beats OMP=6 by 1.63x (257s vs 419s, 50-gen sm0.cfg run).
      HT helps: AIC+NeuralFoil has enough memory latency to benefit from hyperthreading.
      Use 12 threads. 150-gen runtime ~12.9 min. README "~6-7 min" too optimistic but 12T is correct.

## Next features
- [x] Add dynamic stability metrics (Dutch-roll + phugoid, Roadmap #2).
- [ ] Expand the geometry search space (Roadmap #3: mid-span CST breakpoint, 3 profiles).
- [ ] Add visualization, checkpoint/restart, sensitivity sweeps, and export helpers (Roadmap #4).
