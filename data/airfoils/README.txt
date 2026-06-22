Drop airfoil coordinate files (*.dat) here to seed the optimizer.

- Selig format (single TE -> LE -> TE loop) and Lednicer format (count header +
  upper block + lower block) are both auto-detected.
- Good flying-wing seeds are reflexed sections (e.g. MH60, MH45, EH 1.5/2.0).
- Coordinates may be at any chord scale; they are normalized to [0,1] on load.

NACA 4-digit seeds are generated automatically from the `seed_naca` list in
config/baseline.cfg (no file needed).

Each seed is fit to CST (4 upper + 4 lower weights), injected into the initial
GA population, and the CST gene bounds are widened to contain it so the GA can
morph around it.
