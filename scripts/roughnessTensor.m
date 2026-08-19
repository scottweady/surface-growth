
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Homogenizes surface roughness into an effective permeability tensor K and
% area element J by solving the periodic cell problem on a unit-period
% roughness patch, then uses the homogenized growth flow to get the resulting
% growth rate as a function of roughness amplitude, for both a 1D (wave) and
% a 2D (checkerboard) roughness pattern.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../src')

close all

lambda = 1;

% ── Roughness amplitude sweep ──────────────────────────────────────────────

aspan = linspace(0, 1, 16);
Na = length(aspan);
hom = zeros(Na, 3);

% Cell problem discretization (periodic, unit box)
N = 31;
x = trigpts(N, [0 1]);
y = trigpts(N, [0 1]);
[y, x] = meshgrid(y, x);

for na = 1 : Na

  a = aspan(na);

  % 1D roughness: sinusoidal wave in one direction
  h1 = a * sin(2 * pi * x);
  [K, J] = calcK(h1);
  dA1dt = growthFlowHomogenized(lambda, K, J);

  % 2D roughness: checkerboard pattern
  h2 = a * sin(2 * pi * x) .* sin(2 * pi * y);
  [K, J] = calcK(h2);
  dA2dt = growthFlowHomogenized(lambda, K, J);

  hom(na, :) = [a dA1dt / pi dA2dt / pi];

end

% ── Cell problem ─────────────────────────────────────────────────────────────
% Solves the periodic elliptic cell equations for a rough surface graph
% h(x,y) and returns the homogenized permeability tensor K and area
% element J.
function [K, J] = calcK(h)

  N = length(h);
  D = diffmat(N, 1, 'periodic', [0 1]);
  I = speye(N);
  Dx = kron(I, D);
  Dy = kron(D, I);

  h  = h(:);
  hx = Dx * h;
  hy = Dy * h;

  % Metric tensor of the surface graph and its inverse
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

  % Cell equations with a mean-zero constraint
  L = -Dx * (JiG11 .* Dx) - Dx * (JiG12 .* Dy) - Dy * (JiG21 .* Dx) - Dy * (JiG22 .* Dy);
  e = ones(length(L), 1);
  L = [L e; e' 0];

  f1 = Dx * JiG11 + Dy * JiG12;
  f2 = Dx * JiG12 + Dy * JiG22;

  w1 = L \ [f1; 0];
  w2 = L \ [f2; 0];

  w1 = w1(1 : end - 1);
  w2 = w2(1 : end - 1);

  % Integrate cell solutions to get the effective (homogenized) tensor
  k11 = sum(JiG11 + JiG11 .* (Dx * w1) + JiG12 .* (Dy * w1)) / sum(J);
  k12 = sum(JiG12 + JiG11 .* (Dx * w2) + JiG12 .* (Dy * w2)) / sum(J);
  k21 = sum(JiG21 + JiG21 .* (Dx * w1) + JiG22 .* (Dy * w1)) / sum(J);
  k22 = sum(JiG22 + JiG21 .* (Dx * w2) + JiG22 .* (Dy * w2)) / sum(J);

  K = [k11 k12; k21 k22];
  J = sum(J) / N^2;

end
