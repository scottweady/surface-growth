
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulates a growing domain on a flat (possibly anisotropic/homogenized)
% surface. Starting from a circular boundary, the boundary is advanced with
% RK4 using the velocity field from growth, remeshed to equal arclength
% at every step, and the domain's surface area is tracked over time. Runs
% until the area reaches 10x its initial value.
%
% Designed to be called as a script with a shared workspace (e.g. by a
% convergence study), so that lambda, K, dt, tf, doPlot can be pre-set by
% the caller before it runs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ── Parameters ───────────────────────────────────────────────────────────────
% Fall back to defaults for anything not already set in the workspace.

if ~exist('lambda', 'var')
  lambda = 1;
end

if ~exist('K', 'var')
    K = [1.0 0.0;
         0.0 1.0];
end

J = 1.0 / K(1, 1);

if ~exist('doPlot', 'var')
  doPlot = true;
end

if doPlot
  close all
  fig = journalFig([6.75 3], 2);
end

% Flat domain (K represents an effective/homogenized medium)
h = @(x, y) 0 * x;

f          = @(u, v) [u v h(u, v)];
discretize = @(xb, yb, D) surfacemesh.image(f, surfacemesh.ellipticgrid(xb, yb, D));

% ── Mesh ─────────────────────────────────────────────────────────────────────
% Parametric disk with an initial boundary of unit physical area
% (R0 = 1/sqrt(J), matching the convention in growthHomogenized.m, since
% the homogenized cell's true area element is J times the flat domain's).

Nref = 2;
Mb   = 4 * Nref;
P    = 12;
D    = surfacemesh.disk(P, Nref);

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

% ── Main loop ────────────────────────────────────────────────────────────────
% RK4 in time; the boundary curve is remeshed to equal arclength after each
% full step, and the growth flow is re-solved at each RK4 stage since the
% domain shape changes at every stage.

Xb = {};
Yb = {};
AA = [];

for nt = 1 : Nt

  t = dt * (nt - 1);

  Omega      = discretize(xb, yb, D);
  A          = surfacearea(Omega);
  [U1, p]    = growth(Omega, lambda, K);
  [ub, vb]   = extractBoundaryVelocity(U1, Mb);

  Xb{nt} = xb;
  Yb{nt} = yb;
  AA     = [AA; t, A];

  if doPlot
    clc
    fprintf('t = %.3f / %.3f\n', t, tf)
    fprintf('|Omega| / |Omega0| = %.2f\n', A / A0)

    subplot(1, 2, 1)
    xbcat = cellfun(@(x) x, xb, 'UniformOutput', false);
    ybcat = cellfun(@(y) y, yb, 'UniformOutput', false);
    plot([xbcat{:}], [ybcat{:}], 'k'), hold on
    axis equal off
    title(sprintf('$t = %.2f$', t), 'Interpreter', 'latex')

    subplot(1, 2, 2)
    cla
    plot(AA(:,1), AA(:,2) / A0, 'k')
    xlabel('$t$', 'Interpreter', 'latex')
    ylabel('$|\Omega| / |\Omega_0|$', 'Interpreter', 'latex')
    drawnow
  end

  if A >= 10 * A0
    break
  end

  % RK4 stages (remesh the boundary curve before evaluating the velocity at
  % each intermediate stage)
  k1u = ub;
  k1v = vb;

  xb2 = cell(Mb, 1);  yb2 = cell(Mb, 1);
  for mb = 1:Mb
    xb2{mb} = xb{mb} + (dt/2) * k1u{mb};
    yb2{mb} = yb{mb} + (dt/2) * k1v{mb};
  end

  Omega2     = discretize(xb2, yb2, D);
  [U2, ~]    = growth(Omega2, lambda, K);
  [k2u, k2v] = extractBoundaryVelocity(U2, Mb);

  xb3 = cell(Mb, 1);  yb3 = cell(Mb, 1);
  for mb = 1:Mb
    xb3{mb} = xb{mb} + (dt/2) * k2u{mb};
    yb3{mb} = yb{mb} + (dt/2) * k2v{mb};
  end

  Omega3     = discretize(xb3, yb3, D);
  [U3, ~]    = growth(Omega3, lambda, K);
  [k3u, k3v] = extractBoundaryVelocity(U3, Mb);

  xb4 = cell(Mb, 1);  yb4 = cell(Mb, 1);
  for mb = 1:Mb
    xb4{mb} = xb{mb} + dt * k3u{mb};
    yb4{mb} = yb{mb} + dt * k3v{mb};
  end

  Omega4     = discretize(xb4, yb4, D);
  [U4, ~]    = growth(Omega4, lambda, K);
  [k4u, k4v] = extractBoundaryVelocity(U4, Mb);

  for mb = 1:Mb
    xb{mb} = xb{mb} + (dt/6) * (k1u{mb} + 2*k2u{mb} + 2*k3u{mb} + k4u{mb});
    yb{mb} = yb{mb} + (dt/6) * (k1v{mb} + 2*k2v{mb} + 2*k3v{mb} + k4v{mb});
  end

  [xb, yb] = equalArclength(xb, yb, h);

end

Omega = discretize(xb, yb, D);
A = surfacearea(Omega);

% ── Local functions ──────────────────────────────────────────────────────────

% Extracts the boundary velocity of each segment from the surfacefunv U,
% patching segment endpoints together so the boundary stays closed.
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

% Reparametrizes the boundary curve(s) to equal arclength spacing by solving
% for the speed profile phi that maps uniform Chebyshev points to uniform
% arclength, refining the number of subdivisions until the endpoint matches.
function [xb, yb] = equalArclength(xb, yb, h, varargin)

  if isempty(varargin)

    for k = 1 : length(xb)
      [xb, yb] = equalArclength(xb, yb, h, k);
    end

  else

    k = varargin{1};

    x = xb{k};
    y = yb{k};
    z = h(x, y);

    N = length(x);

    [t, w] = chebpts(N, [0 1]);
    w = w(:);
    Dmat = diffmat(N, 1, [0 1]);

    dt = t(2 : end) - t(1 : end - 1);

    xf = x(end);
    yf = y(end);
    zf = z(end);

    ds = sqrt( (Dmat * x).^2 + (Dmat * y).^2 + (Dmat * z).^2 );
    L = sum(ds .* w);

    fphi = @(phi) equalArclengthF(phi, L, ds, t);

    M = 1;
    Mmax = 32;
    tol = 1e-8;

    err = 10 * tol;

    while err > tol && M <= Mmax

      phi = zeros(N, 1);

      for n = 1 : N - 1
        phi(n + 1) = equalArclengthRK4(phi(n), fphi, dt(n) / M, dt(n));
      end

      B = barymat(phi, t);
      xp = B * x;
      yp = B * y;
      zp = B * z;

      err = sqrt( (xp(end) - xf)^2 + (yp(end) - yf)^2 + (zp(end) - zf)^2);

      M = 2 * M;

    end

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
