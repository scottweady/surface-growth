
addpath('../utils')
addpath('../src')
addpath('../scripts')

%% Convergence of the direct (rough-surface) solve to the homogenized solve
% as the roughness wavelength eps -> 0.
%
% The homogenized reference (calcK -> K, J -> growthFlowHomogenized) matches
% exactly how the manuscript's "hom" growth-rate curve was computed (see
% figs/scripts/fig_homogenization.m). For the direct solve, build the
% actual rough surface and find its radius R via fzero so the true
% surface area is exactly pi, matching the area convention used there.

lambda = 1;
a      = 0.1;
h = @(xi, eta) a * sin(2 * pi * xi) .* sin(2 * pi * eta);
% Roughness cell problem: 2D wave, unit-period surface shape -> effective
% tensor
N = 31;
xi  = trigpts(N, [0 1]);
eta = trigpts(N, [0 1]);
[eta, xi] = meshgrid(eta, xi);
Hcell = h(xi, eta);
[K, J] = calcK(Hcell);

% Fixed mesh (unit disk), reused for both the homogenized reference below
% and the direct solve
P    = 12;
Nref = 12;
D    = surfacemesh.disk(P, Nref);

R_hom = 1/sqrt(J);
Omega_hom = surfacemesh.image(@(u, v) [R_hom * u, R_hom * v, 0 * u], D);
[~, p_hom] = growthFlow(Omega_hom, lambda, K);
dOmegadt_hom = J * integral(exp(-lambda * p_hom));

% Roughness wavelength sweep
epsspan = [2.0 1.0 0.5 0.25 0.125];
dOmegadt    = zeros(size(epsspan));

for k = 1:length(epsspan)
  veps = epsspan(k);
  clc
  fprintf('Direct-to-homogenized convergence...\n')
  fprintf('%s eps = %g\n', progressBar(k / length(epsspan)), veps)

  Omega = surfacemesh.image(@(u, v) [R_hom * u, R_hom * v, veps * h(R_hom * u/veps, R_hom * v/veps)], D);
  [~, p] = growthFlow(Omega, lambda);
  dOmegadt(k) = integral(exp(-lambda * p));
end

%% Plot eps convergence

journal_figure([3.75 2.5], 2);

err = abs(dOmegadt - dOmegadt_hom);

loglog(epsspan, err, 'ko-', 'DisplayName', 'error'), hold on
loglog(epsspan, 10 * err(end) * (epsspan / epsspan(end)).^1, '--', 'Color', [0.2 0.2 0.2], 'DisplayName', '$O(\varepsilon)$')

xlim([0.5 * min(epsspan) 2 * max(epsspan)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('roughness wavelength, $\varepsilon$')
ylabel('error, $|\dot A_\varepsilon - \dot\Omega_{\rm hom}|$')
legend('Location', 'best')

print('output/fig-homogenization-eps-convergence', '-dpdf')

function [K, J] = calcK(h)

  N  = length(h);
  D  = diffmat(N, 1, 'periodic', [0 1]);
  I  = speye(N);
  Dx = kron(I, D);
  Dy = kron(D, I);

  h  = h(:);
  hx = Dx * h;
  hy = Dy * h;

  g11 = 1 + hx.^2;
  g12 = hx .* hy;
  g21 = hy .* hx;
  g22 = 1 + hy.^2;
  detg = g11 .* g22 - g12 .* g21;

  ig11 =  g22 ./ detg;
  ig12 = -g12 ./ detg;
  ig21 = -g21 ./ detg;
  ig22 =  g11 ./ detg;

  J = sqrt(detg);

  JiG11 = J .* ig11;
  JiG12 = J .* ig12;
  JiG21 = J .* ig21;
  JiG22 = J .* ig22;

  L = -Dx*(JiG11.*Dx) - Dx*(JiG12.*Dy) - Dy*(JiG21.*Dx) - Dy*(JiG22.*Dy);
  e = ones(length(L), 1);
  L = [L e; e' 0];

  f1 = Dx*JiG11 + Dy*JiG12;
  f2 = Dx*JiG12 + Dy*JiG22;

  w1 = L \ [f1; 0];  w1 = w1(1:end-1);
  w2 = L \ [f2; 0];  w2 = w2(1:end-1);

  k11 = sum(JiG11 + JiG11.*(Dx*w1) + JiG12.*(Dy*w1)) / sum(J);
  k12 = sum(JiG12 + JiG11.*(Dx*w2) + JiG12.*(Dy*w2)) / sum(J);
  k21 = sum(JiG21 + JiG21.*(Dx*w1) + JiG22.*(Dy*w1)) / sum(J);
  k22 = sum(JiG22 + JiG21.*(Dx*w2) + JiG22.*(Dy*w2)) / sum(J);

  K = [k11 k12; k21 k22];
  J = sum(J) / N^2;

end
