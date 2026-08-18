
close all
addpath('../src')

% Stress sensitivity
lambda = 1;

% Perturbation parameters
veps = 0.25;

% Roughness amplitude span
aspan = linspace(0, 1, 10);

% Storage for growth rates
Gspan = zeros(length(aspan), 2);

% Mesh order
P = 12; 

% Construct parameteric domain
D = surfacemesh.disk(P, 4);

for na = length(aspan) : -1 : 1

  % Roughness amplitude
  a = aspan(na);

  % Surface parameterization
  k = 2 * pi / veps;
  h = @(x, y) a * veps * cos(k * x) .* cos(k * y);

  % Construct graph
  R = fzero(@(R) obj(R, h, D) - pi, 1); % Find R such that the area of the surface is pi
  Omega = surfacemesh.image(@(u, v) [R * u R * v h(R * u, R * v)], D);

  % Solve for velocity and pressure
  [U, p] = darcy(Omega, lambda);

  % Compute growth rate
  Adot = integral(exp(-lambda * p));

  % Store results
  Gspan(na, :) = [a Adot];

  cla
  surf(p)
  drawnow

  clc, fprintf('Solving %d of %d\n', na, length(aspan))

end

close all
plot(Gspan(:, 1), Gspan(:, 2), 'o-')

function A = obj(R, h, D)
  Omega = surfacemesh.image(@(u, v) [R * u R * v h(R * u, R * v)], D);
  A = surfacearea(Omega);
end