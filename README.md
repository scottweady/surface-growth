<p align="center">
  <img src="figs/hero.png" alt="A colony spreading over a sphere, with growth fronts at equally spaced times shaded from light (early) to dark (late)" width="320">
</p>

# Geometric regulation of stress-limited growth on surfaces

We analyze a model for microbial communities growing on surfaces in which individual growth is inhibited by collective growth-induced pressure. Asymptotic and numerical solutions on the sphere and homogenization of periodic surfaces show surface geometry can be tuned to either slow or accelerate growth. We use these insights to identify classes of surfaces for which the growth rate can be made arbitrarily small or arbitrarily close to the stress-free growth rate. These results show stress-sensing is coupled to surface curvature and topography at scales much larger than the individual, and provide design principles for controlling the growth of stress-sensitive microbial communities.

## Dependencies

- MATLAB (tested with R2022a)
- [chebfun](https://www.chebfun.org/): see the site for installation instructions.
- [surfacefun](https://github.com/scottweady/surfacefun) (fork of danfortunato/surfacefun): clone it and add it to your MATLAB path.

Both must be on the MATLAB path before running anything in this repo. Run `setup.m` to add this repo's own `utils` and `src` folders to the path.

## Quick start

From the repository root in MATLAB:

```matlab
setup                          % add src/ and utils/ to the path

cd data,  makedata, cd ..      % solve and save figure data (skips existing .mat files)
cd figs,  makefigs, cd ..      % write figs/output/fig2.pdf ... fig6.pdf
cd test,  runtests, cd ..      % convergence studies, written to test/output/
```

The precomputed data in `data/` is already committed, so `makefigs` can be run directly without re-solving. To force a re-solve, delete the corresponding `.mat` file before running `makedata`.

## Structure

- `src/`: core solvers
  - `growth`: surface-sensing growth flow on a general `surfacemesh`, optionally anisotropic
  - `growthHomogenized`: homogenized problem on a disk (polar spectral discretization)
  - `growthSphere`: axisymmetric growth on a spherical cap
  - `roughnessTensor`: effective tensor $K$ and area element $J$ of a periodic rough surface
- `scripts/`: higher-level simulations built on `src/`
  - `growthEvolve`: moving-boundary growth on a flat or homogenized surface
  - `sphereProfiles`, `sphereGrowth`, `sphereGrowthRate`: spherical-cap growth
  - `roughnessGrowthRate`, `movingBoundaryRoughness`: roughness homogenization
  - `planarStability`: linear stability of the growing disk to non-axisymmetric perturbations
- `data/`: figure data. `data/makedata.m` runs the `scripts/` solves and saves their output as `fig2.mat`, ..., `fig6/isotropic.mat`, `fig6/anisotropic.mat`.
- `figs/`: figure scripts, run via `figs/makefigs.m`. Each `figN.m` loads its precomputed data from `data/`. `figs/scripts/hero.m` makes the image above.
- `test/`: convergence studies, run via `test/runtests.m`
- `utils/`: shared plotting and numerical helpers
