Third-party executables (vendored) — drop these here:

  tools/bin/xfoil.exe   standalone XFOIL (Fortran; NO Python needed)
                        https://web.mit.edu/drela/Public/web/xfoil/
                        Used OFFLINE by build_surrogate.exe to calibrate the
                        viscous polar tables. The runtime aeroanalyzer.exe never
                        calls it.

  tools/bin/avl.exe     standalone AVL (vortex lattice)
                        https://web.mit.edu/drela/Public/web/avl/
                        The cross-check oracle for the reference aero model and,
                        later, the Milestone-3 panel solver. Open out/*.avl in it.

Config (config/baseline.cfg) already points here:
  xfoil_exe = tools/bin/xfoil.exe
  avl_exe   = tools/bin/avl.exe

Notes:
- Run the tools from the PROJECT ROOT so these relative paths resolve
  (build.ps1 and the executables already use the project root as cwd).
- Paths with spaces are fine to override to an absolute path in the config.
- These are large binaries you download separately; they are not part of the
  source tree (don't commit them).
