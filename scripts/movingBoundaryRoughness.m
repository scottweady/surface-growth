
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Runs the moving-boundary growth simulation (growthEvolve) with the
% homogenized permeability tensor K of two roughness patterns: a 2D
% (checkerboard) pattern, whose symmetry makes K isotropic, and a 1D (wave)
% pattern, whose K is anisotropic.
%
% Run via data/makedata.m; output is used by fig6.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 1;
a = 1;   % roughness amplitude

doPlot = true;
dt = 0.05;
tf = 50;

% ── Isotropic (2D roughness) ─────────────────────────────────────────────────

fprintf('Moving boundary growth: isotropic (2D roughness)...\n')

h2 = @(x, y) a * sin(2*pi*x) .* sin(2*pi*y);
K = roughnessTensor(h2);

growthEvolve

iso.Xb = Xb;  iso.Yb = Yb;  iso.AA = AA;  iso.A0 = A0;  iso.J = J;

% ── Anisotropic (1D roughness) ───────────────────────────────────────────────

fprintf('Moving boundary growth: anisotropic (1D roughness)...\n')

h1 = @(x, y) a * sin(2*pi*x);
K = roughnessTensor(h1);

growthEvolve

ani.Xb = Xb;  ani.Yb = Yb;  ani.AA = AA;  ani.A0 = A0;  ani.J = J;
