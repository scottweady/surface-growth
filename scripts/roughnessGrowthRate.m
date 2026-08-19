
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solves (a) the homogenized growth rate as a function of roughness
% amplitude for 1D and 2D roughness patterns, and (b) the flow fields (for
% a fixed roughness amplitude) on flat, 1D-wave, and 2D-wave surfaces from
% the homogenized growth equations.
%
% Run via data/makedata.m; output is used by fig5.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 1;
a      = 1;   % roughness amplitude for the flow field snapshots

% ── Growth rate vs amplitude ─────────────────────────────────────────────────

aspan = linspace(0, 1, 16);
Na    = length(aspan);
hom   = zeros(Na, 3);

for na = 1:Na
  ai = aspan(na);

  clc
  fprintf('Homogenized growth rate vs. roughness amplitude...\n')
  fprintf('%s a = %g\n', progressBar(na / Na), ai)

  h1 = @(x, y) ai * sin(2*pi*x);
  [K1, J1] = roughnessTensor(h1);
  dA1dt = growthHomogenized(lambda, K1, J1);

  h2 = @(x, y) ai * sin(2*pi*x) .* sin(2*pi*y);
  [K2, J2] = roughnessTensor(h2);
  dA2dt = growthHomogenized(lambda, K2, J2);

  hom(na, :) = [ai, dA1dt/pi, dA2dt/pi];
end

% ── Flow field snapshots ─────────────────────────────────────────────────────

h_list = {@(x, y) 0 * x, ...
          @(x, y) a * sin(2*pi*x), ...
          @(x, y) a * sin(2*pi*x) .* sin(2*pi*y)};

X = cell(3,1);  Y = cell(3,1);
U = cell(3,1);  V = cell(3,1);
G = cell(3,1);

for k = 1:3
  [K, J] = roughnessTensor(h_list{k});
  [~, u, v, p, xc, yc] = growthHomogenized(lambda, K, J);
  X{k} = [xc, xc(:,1)];
  Y{k} = [yc, yc(:,1)];
  U{k} = [u,  u(:,1)];
  V{k} = [v,  v(:,1)];
  G{k} = [exp(-lambda*p), exp(-lambda*p(:,1))];
end
