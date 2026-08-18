
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Renders the growth rate on the actual rough surface: the image of the
% disk of radius R = 1/sqrt(J*) under (xi,eta) -> (xi, eta, eps*h(xi/eps,eta/eps)).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

addpath('../../src')
addpath('../../utils')

lambda  = 1;
a       = 1;     % roughness amplitude
eps_val = 0.1;   % roughness period (number of bumps ~ 2R/eps_val)
ng      = 300;   % Cartesian grid resolution

h_fns  = { @(xi, eta) zeros(size(xi)), ...
           @(xi, eta) a * sin(2*pi*xi), ...
           @(xi, eta) a * sin(2*pi*xi) .* sin(2*pi*eta) };
labels = {'$h = 0$', ...
          '$h = a\sin(2\pi\xi'')$', ...
          '$h = a\sin(2\pi\xi'')\sin(2\pi\eta'')$'};

fig = figure('Units', 'inches', 'Position', 2*[1 1 10 4]);
clim_vals = [0.8 1.0];

for k = 1:3

  h_fn = h_fns{k};

  % ── Roughness tensor on unit cell ────────────────────────────────────────
  N = 31;
  [eta_cell, xi_cell] = meshgrid(trigpts(N, [0 1]), trigpts(N, [0 1]));
  [K, J] = calcK(h_fn(xi_cell, eta_cell));
  R = 1 / sqrt(J);

  sqrt(K(1,1) * K(2,2)) * J

  % ── Solve homogenized equations on the disk ───────────────────────────────
  [~, ~, ~, p, xc, yc] = darcyHomogenized(lambda, K, J);

  % ── Interpolate p onto a fine Cartesian grid ──────────────────────────────
  % darcyHomogenized uses a double cover (r in [-1,1]); use all points since
  % both halves sample the same physical disk.
  Fp = scatteredInterpolant(xc(:), yc(:), p(:), 'natural', 'none');

  xg   = linspace(-R, R, ng);
  [Xg, Yg] = meshgrid(xg, xg);
  mask = Xg.^2 + Yg.^2 > R^2;

  pg      = Fp(Xg, Yg);
  pg(mask) = NaN;

  % ── Lift to 3D surface ────────────────────────────────────────────────────
  Zg       = eps_val * h_fn(Xg / eps_val, Yg / eps_val);
  Zg(mask) = NaN;

  gamma = exp(-lambda * pg);

  % ── Plot ──────────────────────────────────────────────────────────────────
  ax = subplot(1, 3, k);
  surf(Xg, Yg, Zg, gamma, 'EdgeColor', 'none')
  axis equal
  axis([-1 1 -1 1 -eps_val*1.2 eps_val*1.2])
  colormap(ax, cmocean('algae'))
  clim(clim_vals)
  shading interp

  xlabel('$\xi$',  'Interpreter', 'latex')
  ylabel('$\eta$', 'Interpreter', 'latex')
  zlabel('$z$',    'Interpreter', 'latex')
  title(labels{k}, 'Interpreter', 'latex')

  view([-35 30])
  ax.Box = 'off';
  ax.TickLabelInterpreter = 'latex';

  if k == 3
    cb = colorbar;
    cb.Label.String      = 'growth rate, $\gamma(p)$';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    cb.Ticks = linspace(clim_vals(1), clim_vals(2), 5);
  end

end

% ── Local function ────────────────────────────────────────────────────────────

function [K, J] = calcK(h)

  N  = length(h);
  D  = diffmat(N, 1, 'periodic', [0 1]);
  I  = speye(N);
  Dx = kron(I, D);
  Dy = kron(D, I);

  h  = h(:);
  hx = Dx * h;
  hy = Dy * h;

  g11 = 1 + hx.^2;
  g12 = hx .* hy;
  g21 = hy .* hx;
  g22 = 1 + hy.^2;
  detg = g11 .* g22 - g12 .* g21;

  ig11 =  g22 ./ detg;
  ig12 = -g12 ./ detg;
  ig21 = -g21 ./ detg;
  ig22 =  g11 ./ detg;

  J = sqrt(detg);

  JiG11 = J .* ig11;
  JiG12 = J .* ig12;
  JiG21 = J .* ig21;
  JiG22 = J .* ig22;

  L = -Dx*(JiG11.*Dx) - Dx*(JiG12.*Dy) - Dy*(JiG21.*Dx) - Dy*(JiG22.*Dy);
  e = ones(length(L), 1);
  L = [L e; e' 0];

  f1 = Dx*JiG11 + Dy*JiG12;
  f2 = Dx*JiG12 + Dy*JiG22;

  w1 = L \ [f1; 0];  w1 = w1(1:end-1);
  w2 = L \ [f2; 0];  w2 = w2(1:end-1);

  k11 = sum(JiG11 + JiG11.*(Dx*w1) + JiG12.*(Dy*w1)) / sum(J);
  k12 = sum(JiG12 + JiG11.*(Dx*w2) + JiG12.*(Dy*w2)) / sum(J);
  k21 = sum(JiG21 + JiG21.*(Dx*w1) + JiG22.*(Dy*w1)) / sum(J);
  k22 = sum(JiG22 + JiG21.*(Dx*w2) + JiG22.*(Dy*w2)) / sum(J);

  K = [k11 k12; k21 k22];
  J = sum(J) / N^2;

end
