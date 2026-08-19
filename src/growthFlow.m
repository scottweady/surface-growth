function [U, p] = growthFlow(Omega, lambda, K)
% GROWTHFLOW Solve the (anisotropic) surface-sensing growth flow problem
%
% Solves the nonlinear anisotropic growth flow problem
%   -div(K grad p) - exp(-λp) = 0
% on a surfacemesh Omega using Newton's method, then recovers the velocity
% field from the pressure gradient.
%
% INPUTS:
%   Omega  - Surface domain structure (surfacemesh)
%   lambda - Non-dimensional parameter controlling the exponential nonlinearity
%   K      - 2x2 stiffness (permeability) tensor, K = [K11 K12; K21 K22]
%            (optional, defaults to the identity, i.e. the isotropic problem)
%
% OUTPUTS:
%   U - Cartesian velocity field (computed as U = -K*grad(p))
%   p - Pressure field

  if nargin < 3
    K = eye(2);
    pdo = struct('lap', -1);
    divKgradp = @(p) lap(p);
  else
    pdo = struct('dxx', -K(1,1), 'dxy', -K(1,2), 'dyx', -K(2,1), 'dyy', -K(2,2));
    divKgradp = @(p) K(1,1) * diffx(p, 2) + (K(1,2) + K(2,1)) * diffx(diffy(p)) + K(2,2) * diffy(p, 2);
  end

  % ── Initial guess ────────────────────────────────────────────────────────
  % Linear anisotropic problem -div(K grad p) = 1

  L = surfaceop(Omega, pdo, 1);
  p = L.solve(0);

  % ── Newton iteration ─────────────────────────────────────────────────────

  tol = 1e-6;
  err = 10 * tol;
  maxIter = 100;
  iter = 0;

  while err > tol && iter < maxIter
    % Linearized problem: -div(K grad q) + λ*exp(-λp)*q = -(-div(K grad p) - exp(-λp))
    % where q is the Newton correction
    pdo.b = lambda * exp(-lambda * p);
    L = surfaceop(Omega, pdo);
    L.rhs = -divKgradp(p) - exp(-lambda * p);

    q = L.solve(0);
    p = p - q;

    err = norm(q, 'inf');
    iter = iter + 1;
  end

  if err > tol
    fprintf(' growthFlow failed with tolerance %1.4e. Error is %1.4e\n', tol, err)
  end

  % ── Velocity field ───────────────────────────────────────────────────────
  % U = -K*grad(p)

  [px, py, ~] = grad(p);
  u = -(K(1,1) * px + K(1,2) * py);
  v = -(K(2,1) * px + K(2,2) * py);
  U = surfacefunv(u, v, 0 * u);

end
