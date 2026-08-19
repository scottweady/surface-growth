function [K, J] = roughnessTensor(h, N)
% ROUGHNESSTENSOR Homogenize a rough surface graph into an effective tensor
%
% Solves the periodic elliptic cell problem for a rough surface graph
% z = h(x,y) on the unit-period box [0,1]^2, and returns the homogenized
% permeability tensor K and area element J.
%
% INPUTS:
%   h - Surface height, as a function handle h(x,y), periodic on [0,1]^2
%   N - Grid resolution (number of points per dimension, optional,
%       defaults to 31; should be odd to avoid mismatched mode)
%
% OUTPUTS:
%   K - 2x2 homogenized permeability tensor
%   J - Homogenized area element

  if nargin < 2
    N = 31;
  end

  if mod(N, 2) == 0
    warning('roughnessTensor:evenN', 'N = %d is even; must be odd to avoid mismatched mode.', N)
  end

  x = trigpts(N, [0 1]);
  y = trigpts(N, [0 1]);
  [y, x] = meshgrid(y, x);

  D  = diffmat(N, 1, 'periodic', [0 1]);
  I  = speye(N);
  Dx = kron(I, D);
  Dy = kron(D, I);

  hx = Dx * h(x(:), y(:));
  hy = Dy * h(x(:), y(:));

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
