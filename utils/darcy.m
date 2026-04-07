function [U, p] = darcy(Omega, lambda)
% DARCY Solve Darcy flow problem on a surface using Newton's method
%
% This function solves the nonlinear Darcy flow problem:
%   -Δp - exp(-λp) = 0
% on a surfacemesh Omega using Newton's method.
%
% INPUTS:
%   Omega  - Surface domain structure (surfacemesh)
%   lambda - Non-dimensional parameter controlling the exponential nonlinearity
%
% OUTPUTS:
%   U - Cartesian velocity field (computed as U = -∇p)
%   p - Pressure field
%
% ALGORITHM:
%   Uses Newton's method to solve the nonlinear equation iteratively:
%   1. Initialize with linear Laplace solution
%   2. Iterate: solve linearized system for correction
%   3. Update pressure and check convergence
%   4. Compute velocity from pressure gradient
%

  % Solve initial linear problem for pressure
  % -Δp = 1 with appropriate boundary conditions
  pdo = struct('lap', -1);
  L = surfaceop(Omega, pdo, 1);
  p = L.solve(0);

  % Set convergence parameters for Newton iteration
  tol = 1e-8;        % Convergence tolerance
  err = 10 * tol;    % Initialize error above tolerance
  maxIter = 100;     % Maximum number of iterations
  iter = 0;          % Iteration counter
  
  % Newton's method for nonlinear Darcy equation
  while err > tol && iter < maxIter
    % Set up linearized problem: -Δq + λ*exp(-λp)*q = -(-Δp - exp(-λp))
    % where q is the Newton correction
    pdo = struct('lap', -1, 'c', lambda * exp(-lambda * p));
    L = surfaceop(Omega, pdo);
    L.rhs = -lap(p) - exp(-lambda * p);  % Residual of nonlinear equation
    
    % Solve for Newton correction
    q = L.solve(0);
    
    % Update pressure with Newton step
    p = p - q;
    
    % Check convergence using infinity norm of correction
    err = norm(q, 'inf');
    iter = iter + 1;
    fprintf(' Iteration %d, error = %1.4e\n', iter, err)
  end

  % Check if Newton iteration converged
  if err > tol
    fprintf(' Darcy failed with tolerance %1.4e. Error is %1.4e\n', tol, err)
  end

  % Compute the Cartesian velocity field from pressure gradient
  % Darcy's law: U = -∇p
  U = grad(-p);

end