
addpath('../scripts')

%% Convergence of the direct (rough-surface) solve to the homogenized solve
% as the roughness wavelength eps -> 0.
%
% The homogenized reference (roughnessTensor -> K, J -> growthHomogenized) matches
% exactly how the manuscript's "hom" growth-rate curve was computed (see
% figs/scripts/fig_homogenization.m). For the direct solve, build the
% actual rough surface and find its radius R via fzero so the true
% surface area is exactly pi, matching the area convention used there.

lambda = 1;
a      = 0.1;
h = @(xi, eta) a * sin(2 * pi * xi) .* sin(2 * pi * eta);
[K, J] = roughnessTensor(h);

% Fixed mesh (unit disk), reused for both the homogenized reference below
% and the direct solve
P    = 12;
Nref = 12;
D    = surfacemesh.disk(P, Nref);

R_hom = 1/sqrt(J);
Omega_hom = surfacemesh.image(@(u, v) [R_hom * u, R_hom * v, 0 * u], D);
[~, p_hom] = growth(Omega_hom, lambda, K);
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
  [~, p] = growth(Omega, lambda);
  dOmegadt(k) = integral(exp(-lambda * p));
end

%% Plot eps convergence

fig = journalFig([3.75 2.5], 2);

err = abs(dOmegadt - dOmegadt_hom);

loglog(epsspan, err, 'ko-', 'DisplayName', 'error'), hold on
loglog(epsspan, 10 * err(end) * (epsspan / epsspan(end)).^1, '--', 'Color', [0.2 0.2 0.2], 'DisplayName', '$O(\varepsilon)$')

xlim([0.5 * min(epsspan) 2 * max(epsspan)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('roughness wavelength, $\varepsilon$')
ylabel('error, $|\dot A_\varepsilon - \dot\Omega_{\rm hom}|$')
legend('Location', 'best')
title('Direct-to-homogenized: $\varepsilon$ convergence')
