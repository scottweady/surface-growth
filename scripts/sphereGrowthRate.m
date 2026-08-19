
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solves the growth rate as a function of polar angle H on a sphere for a
% range of sensitivities lambda, and the angle of maximal growth H* as a
% function of lambda.
%
% Run via data/makedata.m; output is used by fig4.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lspan = [0 0.01 0.1 1];
Lspan = linspace(0, 1, 100);

hmax = pi/2;

H0 = 1e-2;
Hf = pi-H0;

R = 10;
N = 100;

% ── Growth rate vs. polar angle ──────────────────────────────────────────────

AsaveAll = cell(length(lspan), 1);

for k = 1:length(lspan)
  lambda = lspan(k);

  clc
  fprintf('Growth rate vs. polar angle...\n')
  fprintf('%s lambda = %g\n', progressBar(k / length(lspan)), lambda)

  F = @(H) growthSphere(H, lambda, R, N);
  Asave = [];

  for H = linspace(H0, Hf, 100)
    dAdt = 2 * pi * R^2 * sin(H) * F(H);
    Asave = [Asave; H dAdt];
  end

  AsaveAll{k} = Asave;
end

% ── Angle of maximal growth vs. lambda ───────────────────────────────────────

Amax = [0 pi];

for k = 2:length(Lspan)
  lambda = Lspan(k);

  clc
  fprintf('Angle of maximal growth vs. lambda...\n')
  fprintf('%s lambda = %g\n', progressBar(k / length(Lspan)), lambda)

  H = H0;
  F = @(H) growthSphere(H, lambda, R, N);
  dAdt = @(H) -2 * pi * sin(H) * F(H);

  hmax = fminsearch(dAdt, hmax);

  Amax = [Amax; lambda hmax];

end
