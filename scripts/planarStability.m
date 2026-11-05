%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear and numerical stability of the axisymmetric (disk) solution on the
% plane, Eqs. (30)-(31), to non-axisymmetric perturbations of the interface.
%
% Linear theory (Appendix A): perturbing the boundary of a disk of radius R
% to r = R + delta * cos(k*theta), the amplitude evolves as
% d(delta)/dt = sigma_k delta with
%   sigma_k = -16 c^2 (k^2 - 1) / (1 + 4 c (k - 1)),   c = c(lambda R^2),
% and the base radius as dR/dt = 4 c R, so the shape (relative amplitude
% delta/R) evolves at the rate sigma_k - 4c.
%
% The script
%   1. evaluates sigma_k(Lambda), Lambda = lambda R^2,
%   2. validates it against the full 2D solver (src/growth.m) on slightly
%      perturbed disks.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

% ── Parameters ───────────────────────────────────────────────────────────────

LamSpan = [0 1 10 100 1000];    % Lambda = lambda R^2 for dispersion relations
kSpan   = 1 : 16;               % azimuthal wavenumbers
LamGrid = logspace(-1, 4, 61);  % Lambda grid for sigma_k(Lambda)
kGrid   = [2 3 4 8 16];

LamCheck = [1 10 100];          % 2D solver validation
kCheck   = 2 : 2 : 8;
epsCheck = 1e-3;

% ── 1. Linear theory ─────────────────────────────────────────────────────────

sigmaLin = linearGrowthRate(LamSpan(:), kSpan);
sigmaGrid = linearGrowthRate(LamGrid(:), kGrid);

% ── 2. Validation with the 2D solver ─────────────────────────────────────────
% Solve the full nonlinear problem on r = 1 + eps cos(k theta) with lambda =
% Lambda and project the boundary radial velocity onto cos(k theta).

fprintf('Linear stability: 2D solver validation...\n')

sigma2D  = zeros(length(LamCheck), length(kCheck));
sigmaRef = zeros(length(LamCheck), length(kCheck));

flat = @(u, v) [u v 0 * u];
Nref = 2;
Mb   = 4 * Nref;
P    = 12;
Dmesh = surfacemesh.disk(P, Nref);

for iL = 1 : length(LamCheck)
  for ik = 1 : length(kCheck)

    kk = kCheck(ik);
    [xbc, ybc] = polarBoundary(@(s) 1 + epsCheck * cos(kk * s), Mb, P);
    OmegaC = surfacemesh.image(flat, surfacemesh.ellipticgrid(xbc, ybc, Dmesh));
    UC = growth(OmegaC, LamCheck(iL));

    urc = cell(Mb, 1);
    for mb = 1 : Mb
      ux = UC.components{1}.vals{mb}(:, 1);
      uy = UC.components{2}.vals{mb}(:, 1);
      urc{mb} = (xbc{mb} .* ux + ybc{mb} .* uy) ./ hypot(xbc{mb}, ybc{mb});
    end

    delta = modeCoefficient(cellfun(@hypot, xbc, ybc, 'UniformOutput', false), xbc, ybc, kk);
    sigma2D(iL, ik)  = modeCoefficient(urc, xbc, ybc, kk) / delta;
    sigmaRef(iL, ik) = linearGrowthRate(LamCheck(iL), kk);

    fprintf('  Lambda = %5g, k = %2d: sigma (2D) = %+.5f, sigma (linear) = %+.5f\n', ...
            LamCheck(iL), kk, sigma2D(iL, ik), sigmaRef(iL, ik))
  end
end

% ── Figure ───────────────────────────────────────────────────────────────────

fig = journalFig([6.75 2.5], 2);
cmap = cmocean('algae', length(LamSpan) + 1);
cmap = cmap(2 : end, :);

% (a) Dispersion relation sigma_k vs k
subplot(1, 2, 1)
for iL = 1 : length(LamSpan)
  plot(kSpan, sigmaLin(iL, :), '-', 'Color', cmap(iL, :), ...
       'DisplayName', sprintf('$\\Lambda = %g$', LamSpan(iL))), hold on
  iC = find(LamCheck == LamSpan(iL));
  if ~isempty(iC)
    plot(kCheck, sigma2D(iC, :), 'o', 'Color', cmap(iL, :), 'HandleVisibility', 'off')
  end
end
plot(kSpan, 0 * kSpan, 'k:', 'HandleVisibility', 'off')
xlabel('wavenumber, $k$')
ylabel('growth rate, $\sigma_k$')
legend('Location', 'southwest')
title('(a) lines: linear theory, circles: 2D solver')

% (b) Decay rate -sigma_k vs Lambda
subplot(1, 2, 2)
cmapk = cmocean('dense', length(kGrid) + 1);
cmapk = cmapk(2 : end, :);
for ik = 1 : length(kGrid)
  loglog(LamGrid, -sigmaGrid(:, ik), 'Color', cmapk(ik, :), ...
         'DisplayName', sprintf('$k = %d$', kGrid(ik))), hold on
end
xlim([LamGrid(1) LamGrid(end)])
xlabel('$\Lambda = \lambda R^2$')
ylabel('decay rate, $-\sigma_k$')
legend('Location', 'southwest')
title('(b)')

% ── Local functions ──────────────────────────────────────────────────────────

% Growth rate sigma_k of boundary perturbation mode k for a disk with
% Lambda = lambda R^2 (Appendix A). Lambda is a column vector and k a row
% vector; sigma is length(Lambda) x length(k).
function sigma = linearGrowthRate(Lam, k)

  c = (sqrt(1 + Lam / 2) - 1) ./ (2 * Lam);
  c(Lam == 0) = 1/8;

  sigma = -16 * c.^2 .* (k.^2 - 1) ./ (1 + 4 * c .* (k - 1));

end

% Boundary segments (Chebyshev points on each of Mb arcs) of the polar curve
% r = rfun(s), matching the boundary discretization used in growthEvolve.
function [xb, yb] = polarBoundary(rfun, Mb, P)

  xb = cell(Mb, 1);
  yb = cell(Mb, 1);
  for mb = 1 : Mb
    s = chebpts(P, [0 2*pi/Mb]) + 2*pi*(mb-1)/Mb;
    xb{mb} = rfun(s) .* cos(s);
    yb{mb} = rfun(s) .* sin(s);
  end

end

% Fourier coefficient (1/pi) * int g cos(k theta) d(theta) of a boundary
% function g, where theta is the polar angle. Integrates spectrally in the
% segment parameter, so it is independent of the boundary parametrization.
function a = modeCoefficient(g, xb, yb, k)

  P = length(xb{1});
  [~, w] = chebpts(P, [0 1]);
  D = diffmat(P, [0 1]);

  a = 0;
  for mb = 1 : length(xb)
    theta = unwrap(atan2(yb{mb}, xb{mb}));
    a = a + w * (g{mb} .* cos(k * theta) .* (D * theta)) / pi;
  end

end
