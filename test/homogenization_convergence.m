
addpath('../utils')
addpath('../src')
addpath('../scripts')

lambda = 1;
a = 1;

%% Spatial convergence (roughness tensor cell problem)

close all

NN = 5:10:65;
KK = [];

for k = 1:length(NN)
  N = NN(k);
  clc
  fprintf('Homogenization K-tensor convergence...\n')
  fprintf('%s N = %d\n', progressBar(k / length(NN)), N)
  xi  = trigpts(N, [0 1]);
  eta = trigpts(N, [0 1]);
  [eta, xi] = meshgrid(eta, xi);
  h = a * sin(2 * pi * xi);
  K = calcK(h);
  KK = [KK; K(1, 1)];
end

err = abs(KK(1 : end-1) - KK(end));

semilogy(NN(1 : end-1), err, 'ko-')

xlim([0 max(NN)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('resolution, $N$')
ylabel('error, $|K_{11,N} - K_{11,N_{\max}}|$')

print('output/fig-homogenization-K-convergence', '-dpdf')

%% Spatial convergence (homogenized disk solve)

close all

% Roughness tensor for a representative (anisotropic wave) surface
N  = 31;
xi  = trigpts(N, [0 1]);
eta = trigpts(N, [0 1]);
[eta, xi] = meshgrid(eta, xi);
h = a * sin(2 * pi * xi);
[K, J] = calcK(h);

NN = 4:4:40;
UU = [];

for k = 1:length(NN)
  Nr = NN(k);
  clc
  fprintf('Homogenization disk-solve convergence...\n')
  fprintf('%s Nr = %d\n', progressBar(k / length(NN)), Nr)
  UU = [UU; darcyHomogenized(lambda, K, J, Nr, Nr - 1)];
end

err = abs(UU(1 : end-1) - UU(end));

semilogy(NN(1 : end-1), err, 'ko-')

xlim([0 max(NN)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('resolution, $N$')
ylabel('error, $|\dot\Omega_N - \dot\Omega_{N_{\max}}|$')

print('output/fig-homogenization-spatial-convergence', '-dpdf')

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