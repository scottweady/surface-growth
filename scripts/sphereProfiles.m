
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solves for the pressure and velocity profiles on a spherical cap of fixed
% polar angle H = pi/2, for a range of sphere radii R.
%
% Run via data/makedata.m; output is used by fig2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 1;
N = 100;
H = pi / 2;

Rplt = 10 + [0, 2.^(1:5)];

hAll = cell(length(Rplt), 1);
pAll = cell(length(Rplt), 1);
uAll = cell(length(Rplt), 1);

for k = 1:length(Rplt)
  R = Rplt(k);
  clc
  fprintf('Sphere pressure/velocity profiles...\n')
  fprintf('%s R = %g\n', progressBar(k / length(Rplt)), R)
  [~, p, h, u] = growthSphere(H, lambda, R, N);
  hAll{k} = h;
  pAll{k} = p;
  uAll{k} = u;
end
