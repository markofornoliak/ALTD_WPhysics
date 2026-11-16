# ALTD_WPhysics

Browser wind tunnel, live at [markofornoliak.github.io/ALTD_WPhysics](https://markofornoliak.github.io/ALTD_WPhysics/).

A two-dimensional section of the AeroLab 3D fractional-step MAC scheme: staggered faces, tunnel boundaries, smoothstep inlet ramp from rest, and a density-weighted pressure Poisson equation with \(p = 0\) on the outlet face.

Air is \(\rho = 1.204\,\mathrm{kg/m^3}\), \(\mu = 1.81\times 10^{-5}\,\mathrm{Pa\cdot s}\). The nominal cylinder is the AeroLab point, \(D = 20\,\mathrm{mm}\), \(U = 0.113\,\mathrm{m/s}\), \(\mathrm{Re} \approx 150\).

## What is adapted, and what is not

AeroLab 3D Beta (the v0.8.0 source archive; notes inside still say v0.7.0) is a 3D C++/OpenGL solver with PCG, optional geometric multigrid, and a particle renderer. Its own acceptance note says steady \(C_D\), wake convergence, and the Mac frontend were not passed.

This repository is the browser section, not that native app.

- Advection is semi-Lagrangian, so a 96×48 mesh can still grow a wake. The native code uses conservative upwind or MUSCL.
- Pressure is Jacobi-preconditioned CG with a capped iteration count, not the native \(10^{-7}\) PCG gate.
- \(C_D\) and \(C_L\) are section forces per unit depth (pressure plus a wall-shear estimate). They are not a certified drag coefficient.

Bodies follow the AeroLab set as tunnel silhouettes: cylinder, cube, NACA 0012, funnel, torus section (two disks), helix (four alternating disks).

## Layout

- `src/lib/wind/solver.ts` — MAC step, masks, coefficients
- `src/lib/wind/colormap.ts` — AeroLab color stops
- `src/components/wind-lab.tsx` — tunnel view and controls
- `pages/` — static entry built for GitHub Pages (`vite.pages.config.ts`)
- `site/` — the built page that Actions publishes
