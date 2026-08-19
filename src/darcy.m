function [U, p] = darcy(Omega, lambda, K)
% DARCY Solve (anisotropic) Darcy flow problem on a surface using Newton's method
%
% This function solves the nonlinear anisotropic Darcy flow problem:
%   -div(K grad p) - exp(-λp) = 0
% on a surfacemesh Omega using Newton's method.
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
%
% ALGORITHM:
%   Uses Newton's method to solve the nonlinear equation iteratively:
%   1. Initialize with linear anisotropic problem
%   2. Iterate: solve linearized system for correction
%   3. Update pressure and check convergence
%   4. Compute velocity from pressure gradient
%

  if nargin < 3
    K = eye(2);
    pdo = struct('lap', -1);
    divKgradp = @(p) lap(p);
  else
    pdo = struct('dxx', -K(1,1), 'dxy', -K(1,2), 'dyx', -K(2,1), 'dyy', -K(2,2));
    divKgradp = @(p) K(1,1) * diffx(p, 2) + (K(1,2) + K(2,1)) * diffx(diffy(p)) + K(2,2) * diffy(p, 2);
  end

  % Solve initial linear problem for pressure
  % -div(K grad p) = 1 with appropriate boundary conditions
  

  L = surfaceop(Omega, pdo, 1);
  p = L.solve(0);

  % Set convergence parameters for Newton iteration
  tol = 1e-6;        % Convergence tolerance
  err = 10 * tol;    % Initialize error above tolerance
  maxIter = 100;     % Maximum number of iterations
  iter = 0;          % Iteration counter

  % Newton's method for nonlinear anisotropic Darcy equation
  while err > tol && iter < maxIter
    % Set up linearized problem: -div(K grad q) + λ*exp(-λp)*q = -(-div(K grad p) - exp(-λp))
    % where q is the Newton correction
    pdo.b = lambda * exp(-lambda * p);
    L = surfaceop(Omega, pdo);
    L.rhs = -divKgradp(p) - exp(-lambda * p);  % Residual of nonlinear equation

    % Solve for Newton correction
    q = L.solve(0);

    % Update pressure with Newton step
    p = p - q;

    % Check convergence using infinity norm of correction
    err = norm(q, 'inf');
    iter = iter + 1;
  end

  % Check if Newton iteration converged
  if err > tol
    fprintf(' Darcy failed with tolerance %1.4e. Error is %1.4e\n', tol, err)
  end

  % Compute the Cartesian velocity field from pressure gradient
  % Darcy's law: U = -K*grad(p)
  [px, py, ~] = grad(p);
  u = -(K(1,1) * px + K(1,2) * py);
  v = -(K(2,1) * px + K(2,2) * py);
  U = surfacefunv(u, v, 0 * u);

end
