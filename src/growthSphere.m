function [U, p, h, u] = growthSphere(H, lambda, R, N)
% GROWTHSPHERE Solve the growth flow problem on a spherical cap
%
% Solves the nonlinear growth flow problem on a spherical cap of polar
% angle H and radius R using a 1D Chebyshev discretization of the
% Laplace-Beltrami operator and Newton's method.
%
% INPUTS:
%   H      - Polar angle of the cap boundary
%   lambda - Non-dimensional parameter controlling the exponential nonlinearity
%   R      - Sphere radius
%   N      - Number of Chebyshev points
%
% OUTPUTS:
%   U - Interface (boundary) velocity
%   p - Pressure field on the polar grid
%   h - Polar angle grid
%   u - Velocity field on the polar grid

  f = @(p) [0; exp(-lambda * p(2 : end - 1)); 0];

  % ── Grid and differential operators ──────────────────────────────────────

  h = chebpts(N, [0 H]);
  D = diffmat(N, [0 H]);

  L = -1 ./ (R^2 * sin(h)) .* (D * (sin(h) .* D));

  % Dirichlet/Neumann boundary conditions
  L(1, :) = D(1, :);
  L(end, :) = [zeros(1, N - 1), 1];

  % ── Newton iteration ─────────────────────────────────────────────────────
  % L becomes (nearly) singular as H -> 0 (a vanishing cap) -- expected and
  % harmless, so warnings from the linear solves below are suppressed.

  warnState = warning('off', 'all');

  p = L \ f(zeros(N, 1));
  err = max(abs(L * p - f(p)));

  newton_tol = 1e-6;
  newton_iter = 100;
  iter = 0;

  while err > newton_tol && iter < newton_iter

    F = L * p - f(p);
    DF = L + lambda * f(p) .* speye(N);
    dp = DF \ F;
    p = p - dp;

    err = max(abs(L * p - f(p)));
    iter = iter + 1;

  end

  warning(warnState);

  % ── Velocity field ───────────────────────────────────────────────────────

  U = -D(end, :) * p / R^2;
  u = -D * p / R^2;

end
