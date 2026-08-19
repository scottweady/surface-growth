
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solves growth trajectories on spherical caps for a range of radii R, from
% an initial area A0 to a final area Af, plus the flat (R -> infinity) limit
% solved in closed form.
%
% Run via data/makedata.m; output is used by fig3.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 1;
N = 100;

Rspan = 10 + [0, 2.^(1:5)];

A0 = 4 * pi;
Af = 4 * pi * 10^2;

% ── Growth trajectories ──────────────────────────────────────────────────────

ts = [];
AAall = cell(length(Rspan), 1);

for k = 1:length(Rspan)
  R = Rspan(k);

  clc
  fprintf('Sphere growth trajectories...\n')
  fprintf('%s R = %g\n', progressBar(k / length(Rspan)), R)

  dt = 0.1;
  tf = 100;
  t = 0;

  H0 = acos(1 - A0 / (2 * pi * R^2));
  H = H0;

  Hf = acos(1 - Af / (2 * pi * R^2));

  A = 2 * pi * R^2 * (1 - cos(H0));
  AA = [0, A];

  F = @(Theta) growthSphere(Theta, lambda, R, N);

  Hsave = [0, H];

  while t < tf && H < Hf

    H = rk4(H, F, dt, 1);
    t = t + dt;

    Hsave = [Hsave; t, H];

    A = 2 * pi * R^2 * (1 - cos(H));
    AA = [AA; t, A];

  end

  if R == 10
    Hplt = Hsave;
  end

  % Reverse at high res to find final time
  while H > Hf
    H = rk4(H, F, -0.01*dt, 1);
    t = t - 0.01*dt;
  end

  ts = [ts; R, t];

  AAall{k} = AA;

end

% ── Flat limit (R -> infinity, solved in closed form) ────────────────────────

t = 0;
R = sqrt(A0 / pi);
A = pi * R^2;
AAinf = [0, A];

while A < Af
  c = @(R) (sqrt(1 + lambda * R^2 / 2) - 1) / (2 * lambda * R^2);
  F = @(R) 32 * c(R)^2 * R / (1 - 8 * c(R)^2 * lambda * R^2);
  R = rk4(R, F, dt, 1);
  t = t + dt;
  A = pi * R^2;
  AAinf = [AAinf; t, A];
end

while A > Af
  R = rk4(R, F, -0.01*dt, 1);
  A = pi * R^2;
  t = t - 0.01*dt;
end
ts = [ts; inf, t];
