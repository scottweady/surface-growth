%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear and numerical stability of the axisymmetric (disk) solution on the
% plane, Eqs. (30)-(31), to non-axisymmetric perturbations of the interface.
%
% Linear theory (Appendix A): perturbing the boundary of a disk of radius R
% to r = R (1 + eps cos(m theta)), the relative amplitude evolves as
% d(eps)/dt = sigma_m eps with
%   sigma_m = -4c [1 + 4c (m - 1)(m + 2)] / [1 + 4c (m - 1)],
% where c = c(lambda R^2) and the base radius grows as dR/dt = 4 c R.
% See planarStabilitySymbolic.m for a symbolic derivation.
%
% The script
%   1. evaluates sigma_m for m = 1..10 at several Lambda = lambda R^2,
%   2. validates it against the full 2D solver (src/growth.m) on slightly
%      perturbed disks.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

% ── Parameters ───────────────────────────────────────────────────────────────

LamSpan = [0 1 10 100];         % Lambda = lambda R^2 for dispersion relations
mSpan   = 1 : 10;               % azimuthal wavenumbers

LamCheck = [1 10 100];          % 2D solver validation
mCheck   = [1 2 4 6 8 10];
epsCheck = 1e-3;

% ── 1. Linear theory ─────────────────────────────────────────────────────────

sigmaLin = linearGrowthRate(LamSpan(:), mSpan);

% ── 2. Validation with the 2D solver ─────────────────────────────────────────
% Solve the full nonlinear problem on r = 1 + eps cos(m theta) with lambda =
% Lambda and project the boundary radial velocity u_r onto 1 and cos(m theta).
% Writing r = R (1 + eps cos(m theta)) with absolute amplitude delta = R eps,
%   d(eps)/dt / eps = d(delta)/dt / delta - dR/dt / R,
% where dR/dt is the mean of u_r and d(delta)/dt its cos(m theta) coefficient.

fprintf('Linear stability: 2D solver validation...\n')

sigma2D  = zeros(length(LamCheck), length(mCheck));
sigmaRef = zeros(length(LamCheck), length(mCheck));

flat = @(u, v) [u v 0 * u];
Nref = 2;
Mb   = 4 * Nref;
P    = 12;
Dmesh = surfacemesh.disk(P, Nref);

for iL = 1 : length(LamCheck)
  for im = 1 : length(mCheck)

    mm = mCheck(im);
    [xbc, ybc] = polarBoundary(@(s) 1 + epsCheck * cos(mm * s), Mb, P);
    OmegaC = surfacemesh.image(flat, surfacemesh.ellipticgrid(xbc, ybc, Dmesh));
    UC = growth(OmegaC, LamCheck(iL));

    urc = cell(Mb, 1);
    for mb = 1 : Mb
      ux = UC.components{1}.vals{mb}(:, 1);
      uy = UC.components{2}.vals{mb}(:, 1);
      urc{mb} = (xbc{mb} .* ux + ybc{mb} .* uy) ./ hypot(xbc{mb}, ybc{mb});
    end

    rbc   = cellfun(@hypot, xbc, ybc, 'UniformOutput', false);
    R     = modeCoefficient(rbc, xbc, ybc, 0) / 2;    % mean radius
    Rdot  = modeCoefficient(urc, xbc, ybc, 0) / 2;    % mean radial velocity
    delta = modeCoefficient(rbc, xbc, ybc, mm);
    ddot  = modeCoefficient(urc, xbc, ybc, mm);

    sigma2D(iL, im)  = ddot / delta - Rdot / R;
    sigmaRef(iL, im) = linearGrowthRate(LamCheck(iL), mm);

    fprintf('  Lambda = %5g, m = %2d: sigma (2D) = %+.5f, sigma (linear) = %+.5f\n', ...
            LamCheck(iL), mm, sigma2D(iL, im), sigmaRef(iL, im))
  end
end

% ── Figure ───────────────────────────────────────────────────────────────────

fig = journalFig([3.375 2.5], 2);
cmap = cmocean('algae', length(LamSpan) + 1);
cmap = cmap(2 : end, :);

% Dispersion relation sigma_m vs m
for iL = 1 : length(LamSpan)
  plot(mSpan, sigmaLin(iL, :), '-', 'Color', cmap(iL, :), ...
       'DisplayName', sprintf('$\\Lambda = %g$', LamSpan(iL))), hold on
  iC = find(LamCheck == LamSpan(iL));
  if ~isempty(iC)
    plot(mCheck, sigma2D(iC, :), 'o', 'Color', cmap(iL, :), 'HandleVisibility', 'off')
  end
end
plot([0 mSpan(end)+1], [0 0], 'k:', 'HandleVisibility', 'off')
xlim([0 mSpan(end)+1])
xlabel('wavenumber, $m$')
ylabel('growth rate, $\sigma_m$')
legend('Location', 'southwest')


% ── Local functions ──────────────────────────────────────────────────────────

% Growth rate sigma_m of the relative amplitude of boundary perturbation mode
% m for a disk with Lambda = lambda R^2 (Appendix A). Lambda is a column
% vector and m a row vector; sigma is length(Lambda) x length(m).
function sigma = linearGrowthRate(Lam, m)

  c = (sqrt(1 + Lam / 2) - 1) ./ (2 * Lam);
  c(Lam == 0) = 1/8;

  sigma = -4 * c .* (1 + 4 * c .* (m - 1) .* (m + 2)) ./ (1 + 4 * c .* (m - 1));

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

% Fourier coefficient (1/pi) * int g cos(m theta) d(theta) of a boundary
% function g, where theta is the polar angle. Integrates spectrally in the
% segment parameter, so it is independent of the boundary parametrization.
% For m = 0 this is twice the mean of g.
function a = modeCoefficient(g, xb, yb, m)

  P = length(xb{1});
  [~, w] = chebpts(P, [0 1]);
  D = diffmat(P, [0 1]);

  a = 0;
  for mb = 1 : length(xb)
    theta = unwrap(atan2(yb{mb}, xb{mb}));
    a = a + w * (g{mb} .* cos(m * theta) .* (D * theta)) / pi;
  end

end
