
addpath('../scripts')

lambda = 1;
a = 1;

%% Spatial convergence (roughness tensor cell problem)

fig = journalFig([3.75 2.5], 2);

NN = 5:10:65;
KK = [];

for k = 1:length(NN)
  N = NN(k);
  clc
  fprintf('Homogenization K-tensor convergence...\n')
  fprintf('%s N = %d\n', progressBar(k / length(NN)), N)
  h = @(x, y) a * sin(2 * pi * x);
  K = roughnessTensor(h, N);
  KK = [KK; K(1, 1)];
end

err = abs(KK(1 : end-1) - KK(end));

semilogy(NN(1 : end-1), err, 'ko-', 'DisplayName', 'error')

xlim([0 max(NN)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('resolution, $N$')
ylabel('error, $|K_{11,N} - K_{11,N_{\max}}|$')
legend('Location', 'best')
title('Roughness tensor: spatial convergence')