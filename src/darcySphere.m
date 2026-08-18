
function [U, p, h, u] = darcySphere(H, lambda, R, N)

  % Growth function
  f = @(p) [0; exp(-lambda * p(2 : end - 1)); 0];
  
  % Grid
  h = chebpts(N, [0 H]);

  % Derivative operator
  D = diffmat(N, [0 H]);

  % Laplace-Beltrami operator
  L = -1 ./ (R^2 * sin(h)) .* (D * (sin(h) .* D));

  % Dirichlet boundary conditions
  L(1, :) = D(1, :);
  L(end, :) = [zeros(1, N - 1), 1];
  
  % Initial guess for pressure
  p = L \ f(zeros(N, 1));

  % Initial error
  err = max(abs(L * p - f(p)));

  % Newton solve parameters
  newton_tol = 1e-6;
  newton_iter = 100;
  iter = 0;

  % Newton solve 
  while err > newton_tol && iter < newton_iter
    
    % Residual
    F = L * p - f(p);

    % Jacobian
    DF = L + lambda * f(p) .* speye(N);

    % Newton update
    dp = DF \ F;

    % Update
    p = p - dp;

    % L-inf error
    err = max(abs(L * p - f(p)));

    % Update counter
    iter = iter + 1;

  end

  % Compute interface velocity
  U = -D(end, :) * p / R^2;
  u = -D * p / R^2;

end