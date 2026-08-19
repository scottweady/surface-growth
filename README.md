# surface-growth

Code and data associated with the manuscript "Geometric regulation of stress-limited growth on surfaces".

## Dependencies

- [chebfun](https://www.chebfun.org/) — see the site for installation instructions.
- [surfacefun](https://github.com/scottweady/surfacefun) (fork of danfortunato/surfacefun) — clone it and add it to your MATLAB path.

Both must be on the MATLAB path before running anything in this repo. Run `setup.m` to add this repo's own `utils` and `src` folder to the path.

## Structure

- `src/` — core solvers (`growth`, `growthHomogenized`, `growthSphere`, `roughnessTensor`)
- `scripts/` — higher-level simulations built on `src/`: `growthEvolve` (moving-boundary growth on a flat/homogenized surface), `sphereProfiles`/`sphereGrowth`/`sphereGrowthRate` (spherical-cap growth), `roughnessGrowthRate`/`movingBoundaryRoughness` (roughness homogenization)
- `data/` — figure data; `data/makedata.m` runs the `scripts/` solves and saves their output as `fig2.mat`, ..., `fig6/isotropic.mat`, `fig6/anisotropic.mat` (each run only if its output isn't already on disk)
- `test/` — convergence studies, run via `test/runtests.m`
- `figs/` — figure-generating scripts, run via `figs/makefigs.m` (each `figN.m` loads its precomputed data from `data/`)
- `utils/` — shared plotting and numerical helpers