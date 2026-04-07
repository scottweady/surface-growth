
function [dOmegadt, u, v, p, x, y] = anisotropicDarcy(lambda, K, J)

  Nr = 32;
  Ntheta = 31;

  R = 1 / sqrt(J);

  [r, dr] = chebpts(Nr, [-1 1]);
  [theta, dtheta] = trigpts(Ntheta, [0 2 * pi]);
  dr = dr';
  theta = theta';

  nx = cos(theta);
  ny = sin(theta);

  x = R * r .* cos(theta);
  y = R * r .* sin(theta);

  [theta, r] = meshgrid(theta, r);

  Dr = diffmat(Nr, 1, [-1 1], 'chebkind2');
  Dtheta = diffmat(Ntheta, 1, 'periodic', [0 2 * pi]);

  Dtheta = kron(Dtheta, speye(Nr));
  Dr = kron(speye(Ntheta), Dr);

  Dx = (1 / R) * cos(theta(:)) .* Dr - (1 ./ (R * r(:)).^2) .* sin(theta(:)) .* (R * r(:)) .* Dtheta;
  Dy = (1 / R) * sin(theta(:)) .* Dr + (1 ./ (R * r(:)).^2) .* cos(theta(:)) .* (R * r(:)) .* Dtheta;

  L = -(Dx * (K(1, 1) * Dx + K(1, 2) * Dy) + Dy * (K(2, 1) * Dx + K(2, 2) * Dy));

  % Dirichlet boundary conditions at r = 1, -1
  L(1 : Nr : end, :) = 0;
  L(1 : Nr : end, 1 : Nr : end) = eye(Ntheta);

  L(Nr : Nr : end, :) = 0;
  L(Nr : Nr : end, Nr : Nr : end) = eye(Ntheta);

  g = ones(Nr * Ntheta, 1);
  ib = [1 : Nr : Nr * Ntheta, Nr : Nr : Nr * Ntheta];
  g(ib) = 0;
  p = L \ g;

  tol = 1e-6;
  err = max(abs(L * p - exp(-lambda * p)));
  maxit = 25;
  it = 0;

  while err > tol && it < maxit
    g = exp(-lambda * p);
    g(ib) = 0;
    DF = L - lambda * spdiags(g, 0, Nr * Ntheta, Nr * Ntheta);
    F = L * p - g;
    dp = DF \ F;
    p = p - dp;
    err = norm(F);
    it = it + 1;
  end

  if err > tol
    warning("Darcy solver did not converge: err = %1.4e", err)
  else
    fprintf("Darcy solver converged in %d iterations with error %1.4e\n", it, err)
  end

  p = reshape(p, Nr, Ntheta);

  u = -(K(1, 1) * Dx + K(1, 2) * Dy) * p(:);
  v = -(K(2, 1) * Dx + K(2, 2) * Dy) * p(:);
  u = reshape(u, Nr, Ntheta);
  v = reshape(v, Nr, Ntheta);

  dOmegadt = J * R * sum(dtheta .* (u(end, :) .* nx + v(end, :) .* ny));

  x = x(Nr/2 + 1 : end, :);
  y = y(Nr/2 + 1 : end, :);
  u = u(Nr/2 + 1 : end, :);
  v = v(Nr/2 + 1 : end, :);
  p = p(Nr/2 + 1 : end, :);

end

