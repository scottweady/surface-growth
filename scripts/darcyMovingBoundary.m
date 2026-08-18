
addpath('../src')

% Parameters (may be pre-set in the workspace by a caller, e.g. a convergence
% study -- otherwise default to the values below)
if ~exist('lambda', 'var')
  lambda = 1;
end

% Stiffness (permeability) tensor -- set manually
if ~exist('K', 'var')
    K = [1.0 0.0;
         0.0 1.0];
end

J = 1.0 / K(1, 1);

% Set false to suppress plotting/saving, e.g. when called repeatedly from a
% convergence study
if ~exist('doPlot', 'var')
  doPlot = true;
end
if ~exist('doSave', 'var')
  doSave = true;
end

if doPlot
  close all
end

% Flat domain (K represents an effective/homogenized medium)
h = @(x, y) 0 * x;

f          = @(u, v) [u v h(u, v)];
discretize = @(xb, yb, D) surfacemesh.image(f, surfacemesh.ellipticgrid(xb, yb, D));

% ── Mesh ─────────────────────────────────────────────────────────────────────

Nref = 1;
Mb   = 4 * Nref;
P    = 12;
D    = surfacemesh.disk(P, Nref);

% Initial boundary: circle of unit physical area, R0 = 1/sqrt(J)
% (matches the R = 1/sqrt(J) convention in darcyHomogenized.m, since the
% homogenized cell's true area element is J times the flat domain's)
R0 = 1 / sqrt(J);
xb = cell(Mb, 1);
yb = cell(Mb, 1);
for mb = 1:Mb
  s        = chebpts(P, [0 2*pi/Mb]) + 2*pi*(mb-1)/Mb;
  xb{mb}   = R0 * cos(s);
  yb{mb}   = R0 * sin(s);
end

% ── Time stepping ─────────────────────────────────────────────────────────────

if ~exist('dt', 'var')
  dt = 0.05;
end
if ~exist('tf', 'var')
  tf = 50;
end
Nt = floor(tf / dt);

Omega = discretize(xb, yb, D);
A0    = surfacearea(Omega);

% ── Main loop (RK4 every step, remesh every stage) ────────────────────────────
% Stops once the domain area reaches 10x its initial value.

Xb = {};
Yb = {};
AA = [];

for nt = 1 : Nt

  t = dt * (nt - 1);

  Omega      = discretize(xb, yb, D);
  A          = surfacearea(Omega);
  [U1, p]    = darcy(Omega, lambda, K);
  [ub, vb]   = extractBoundaryVelocity(U1, Mb);

  % Store boundary and area
  Xb{nt} = xb;
  Yb{nt} = yb;
  AA     = [AA; t, A];

  if doPlot
    clc
    fprintf('t = %.3f / %.3f\n', t, tf)
    fprintf('|Omega| / |Omega0| = %.2f\n', A / A0)

    % Plot
    clf
    subplot(1, 2, 1)
%     for m = 1:length(Omega)
%       contourf(Omega.x{m}, Omega.y{m}, exp(-lambda * p.vals{m}), 'EdgeColor', 'none'), hold on
%     end
    xbcat = cellfun(@(x) x, xb, 'UniformOutput', false);
    ybcat = cellfun(@(y) y, yb, 'UniformOutput', false);
    plot([xbcat{:}], [ybcat{:}], 'k-')
    axis equal off
    colorbar
    title(sprintf('$t = %.2f$', t), 'Interpreter', 'latex')

    subplot(1, 2, 2)
    plot(AA(:,1), AA(:,2) / A0, 'k-')
    xlabel('$t$', 'Interpreter', 'latex')
    ylabel('$|\Omega| / |\Omega_0|$', 'Interpreter', 'latex')
    drawnow
  end

  % Stop once the domain reaches 10x its initial area
  if A >= 10 * A0
    break
  end

  % RK4 update (remesh the boundary curve before evaluating the velocity at
  % each intermediate stage)
  k1u = ub;
  k1v = vb;

  xb2 = cell(Mb, 1);  yb2 = cell(Mb, 1);
  for mb = 1:Mb
    xb2{mb} = xb{mb} + (dt/2) * k1u{mb};
    yb2{mb} = yb{mb} + (dt/2) * k1v{mb};
  end

  Omega2     = discretize(xb2, yb2, D);
  [U2, ~]    = darcy(Omega2, lambda, K);
  [k2u, k2v] = extractBoundaryVelocity(U2, Mb);

  xb3 = cell(Mb, 1);  yb3 = cell(Mb, 1);
  for mb = 1:Mb
    xb3{mb} = xb{mb} + (dt/2) * k2u{mb};
    yb3{mb} = yb{mb} + (dt/2) * k2v{mb};
  end

  Omega3     = discretize(xb3, yb3, D);
  [U3, ~]    = darcy(Omega3, lambda, K);
  [k3u, k3v] = extractBoundaryVelocity(U3, Mb);

  xb4 = cell(Mb, 1);  yb4 = cell(Mb, 1);
  for mb = 1:Mb
    xb4{mb} = xb{mb} + dt * k3u{mb};
    yb4{mb} = yb{mb} + dt * k3v{mb};
  end

  Omega4     = discretize(xb4, yb4, D);
  [U4, ~]    = darcy(Omega4, lambda, K);
  [k4u, k4v] = extractBoundaryVelocity(U4, Mb);

  for mb = 1:Mb
    xb{mb} = xb{mb} + (dt/6) * (k1u{mb} + 2*k2u{mb} + 2*k3u{mb} + k4u{mb});
    yb{mb} = yb{mb} + (dt/6) * (k1v{mb} + 2*k2v{mb} + 2*k3v{mb} + k4v{mb});
  end

  % Remesh (equal arclength) every step
  [xb, yb] = equalArclength(xb, yb, h);

end

Omega = discretize(xb, yb, D);
A = surfacearea(Omega);

% ── Plot ──────────────────────────────────────────────────────────────────────

if doPlot
  figure
  plot(AA(:,1), AA(:,2) / A0)
  xlabel('time, $t$', 'Interpreter', 'latex')
  ylabel('surface area, $|\Omega| / |\Omega_0|$', 'Interpreter', 'latex')
  title(sprintf('$\\lambda = %g$', lambda), 'Interpreter', 'latex')
end

% ── Save output ────────────────────────────────────────────────────────────

if doSave
  outdir = 'output';
  if ~exist(outdir, 'dir')
    mkdir(outdir)
  end
  save(fullfile(outdir, 'darcyMovingBoundary.mat'), 'Xb', 'Yb', 'AA', 'K', 'J', 'lambda', 'dt', 'A0')
end

% ── Local function ────────────────────────────────────────────────────────────

function [ub, vb] = extractBoundaryVelocity(U, Mb)

  ub = cell(Mb, 1);
  vb = cell(Mb, 1);

  for mb = 1:Mb
    ub{mb} = U.components{1}.vals{mb}(:, 1);
    vb{mb} = U.components{2}.vals{mb}(:, 1);
  end

  for mb = 1:Mb
    ub{mb}(end) = ub{mod(mb, Mb) + 1}(1);
    vb{mb}(end) = vb{mod(mb, Mb) + 1}(1);
  end

end

function [xb, yb] = equalArclength(xb, yb, h, varargin)

  if isempty(varargin)

    % Loop over segments
    for k = 1 : length(xb)
      [xb, yb] = equalArclength(xb, yb, h, k);
    end

  else

    k = varargin{1};

    % Get coordinates
    x = xb{k};
    y = yb{k};
    z = h(x, y);

    % Number of grid points
    N = length(x);

    % Create Chebyshev grid
    [t, w] = chebpts(N, [0 1]);
    w = w(:);
    Dmat = diffmat(N, 1, [0 1]);

    % Arclength spacing
    dt = t(2 : end) - t(1 : end - 1);

    % Save endpoint to check error
    xf = x(end);
    yf = y(end);
    zf = z(end);

    % Compute arclength
    ds = sqrt( (Dmat * x).^2 + (Dmat * y).^2 + (Dmat * z).^2 );
    L = sum(ds .* w);

    % Right hand side
    fphi = @(phi) equalArclengthF(phi, L, ds, t);

    % Number of subdivisions
    M = 1;

    % Maximum number of subdivisions
    Mmax = 32;

    % Endpoint error tolerance
    tol = 1e-8;

    err = 10 * tol;

    while err > tol && M <= Mmax

      % Initialize speed
      phi = zeros(N, 1);

      % Compute speed
      for n = 1 : N - 1
        phi(n + 1) = equalArclengthRK4(phi(n), fphi, dt(n) / M, dt(n));
      end

      % Update solution
      B = barymat(phi, t);
      xp = B * x;
      yp = B * y;
      zp = B * z;

      % Check error
      err = sqrt( (xp(end) - xf)^2 + (yp(end) - yf)^2 + (zp(end) - zf)^2);

      % Update subdivisions
      M = 2 * M;

    end

    % Reassign value
    B = barymat(phi, t);
    xb{k} = B * x;
    yb{k} = B * y;

  end

end

% Fourth order Runge-Kutta method (used by equalArclength's reparametrization)
function u = equalArclengthRK4(u, f, dt, tf)

  N = round(tf / dt);

  for n = 1 : N

    k1 = f(u);
    k2 = f(u + dt / 2 * k1);
    k3 = f(u + dt / 2 * k2);
    k4 = f(u + dt * k3);

    u = u + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

  end

end

function val = equalArclengthF(phi, L, ds, t)

  B = barymat(phi, t);
  val = L ./ (B * ds);

end
