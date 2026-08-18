
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solves the cell equations for the homogenized roughness tensor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../src')

close all

lambda = 1;

% Range of gamma
aspan = linspace(0, 1, 16);
% aspan = [1.0];

% Number of values
Na = length(aspan);

% Preallocate
hom = zeros(Na, 3);

% Number of modes
N = 31;

% Discretize box
x = trigpts(N, [0 1]);
y = trigpts(N, [0 1]);
[y, x] = meshgrid(y, x);

for na = 1 : Na

  a = aspan(na);
  
  % Domain
  dom = [0 1 0 1];
  
  % Compute roughness tensor and solve
  h1 = a * sin(2 * pi * x);
  [K, J] = calcK(h1);
  dA1dt = darcyHomogenized(lambda, K, J);

  display(K)
  display(J)

  %%

  h2 = a * sin(2 * pi * x) .* sin(2 * pi * y);
  [K, J] = calcK(h2);
  dA2dt = darcyHomogenized(lambda, K, J);

  % Store value
  hom(na, :) = [a dA1dt / pi dA2dt / pi];

end

%%
c_growth = [0.1 0.6 0.3];

% Initialize figure
close all
fig = journal_figure([3.25 2.5], 2);
plot(hom(:, 1), hom(:, 3), '--', 'Color', 'k', 'DisplayName', '$h(\xi'',\eta'') = 0$'), hold on
plot(hom(:, 1), hom(:, 2), 'o-.', 'Color', c_growth, 'DisplayName', '$h(\xi'',\eta'') = a\sin(2 \pi \xi'')$'), hold on
plot(hom(:, 1), hom(:, 3), 's:', 'Color', 1-c_growth, 'DisplayName', '$h(\xi'',\eta'') = a\sin(2 \pi \xi'') \sin(2 \pi \eta'')$')
xlabel('roughness amplitude, $a$')
ylabel('growth rate, $|\dot\Omega|/|\Omega|$')
legend('Location', 'southwest', 'EdgeColor', 'none')
ylim([0.86 0.96])

h1 = 10 * sin(2 * pi * x);
[K, J] = calcK(h1);
[dA2dt, u, v, p, xi, eta] = darcyHomogenized(lambda, K, J);


function [K, J] = calcK(h)

  % Number of discretization points in each dimension
  N = length(h);
  
  % Differentiation matrix
  D = diffmat(N, 1, 'periodic', [0 1]);

  % Form 3D operators
  I = speye(N);
  Dx = kron(I, D);
  Dy = kron(D, I);

  % Flatten
  h = h(:);
  hx = Dx * h;
  hy = Dy * h;

  % Metric tensor
  g11 = 1 + hx.^2;
  g12 = hx .* hy;
  g21 = hy .* hx;
  g22 = 1 + hy.^2;
  
  % Determinant of metric tensor
  detg = g11 .* g22 - g12 .* g21;
  
  % Inverse metric tensor
  ig11 =  g22 ./ detg;
  ig12 = -g12 ./ detg;
  ig21 = -g21 ./ detg;
  ig22 =  g11 ./ detg;
  
  % Jacobian
  J = sqrt(detg);
  
  % Coefficient matrix
  JiG11 = J .* ig11;
  JiG12 = J .* ig12;
  JiG21 = J .* ig21;
  JiG22 = J .* ig22;

  % Construct differential operator
  L = -Dx * (JiG11 .* Dx) - Dx * (JiG12 .* Dy) - Dy * (JiG21 .* Dx) - Dy * (JiG22 .* Dy);

  e = ones(length(L), 1);

  % Mean zero condition
  L = [L e; e' 0];

  % Right hand side of cell equations
  f1 = Dx * JiG11 + Dy * JiG12;
  f2 = Dx * JiG12 + Dy * JiG22;

  % Solve
  w1 = L \ [f1; 0];
  w2 = L \ [f2; 0];

  % Get grid values
  w1 = w1(1 : end - 1);
  w2 = w2(1 : end - 1);
  
  % Integrate to get roughness tensor
  k11 = sum(JiG11 + JiG11 .* (Dx * w1) + JiG12 .* (Dy * w1)) / sum(J);
  k12 = sum(JiG12 + JiG11 .* (Dx * w2) + JiG12 .* (Dy * w2)) / sum(J);
  k21 = sum(JiG21 + JiG21 .* (Dx * w1) + JiG22 .* (Dy * w1)) / sum(J);
  k22 = sum(JiG22 + JiG21 .* (Dx * w2) + JiG22 .* (Dy * w2)) / sum(J);

  % Store roughness tensor
  K = [k11 k12; k21 k22];

  % Area element
  J = sum(J) / N^2;

end
