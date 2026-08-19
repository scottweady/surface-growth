
addpath('../scripts')

lambda = 1;
R = 1;

%% Spatial convergence

figSpatial = journalFig([3.75 2.5], 2);
UU = [];
NN = 4 : 2 : 32;

for k = 1:length(NN)
  N = NN(k);
  clc
  fprintf('Spatial convergence for sphere...\n')
  fprintf('%s N = %d\n', progressBar(k / length(NN)), N)
  UU = [UU; growthSphere(pi/2, lambda, R, N)];
end

err = abs(diff(UU));

semilogy(NN(2 : end), err, 'ko-', 'DisplayName', 'error')

xlim([0 max(NN)+(NN(end)-NN(end-1))])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('resolution, $N$')
ylabel('error, $|U_N - U_{N/2}|$')
legend('Location', 'best')
title('Sphere: spatial convergence')

%% Temporal convergence

figTemporal = journalFig([3.75 2.5], 2);
HH = [];

N = 100;
NNt = 2.^(1 : 8);
F = @(Theta) growthSphere(Theta, lambda, R, N);

for k = 1:length(NNt)
  Nt = NNt(k);
  dt = 1 / Nt;
  clc
  fprintf('Temporal convergence for sphere...\n')
  fprintf('%s dt = %g\n', progressBar(k / length(NNt)), dt)
  HH = [HH; rk4(pi/4, F, dt, Nt)];
end

dt = 1 ./ NNt(2 : end);

err = abs(diff(HH));

loglog(dt, err, 'ko-', 'DisplayName', 'error'), hold on
loglog(dt, 10 * err(1) * (dt/dt(1)).^(4), '--', 'Color', [0.2 0.2 0.2], 'DisplayName', '$O(\Delta t^4)$')

xlim([0.5 * min(dt) 2 * max(dt)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('time step size, $\Delta t$')
ylabel('error, $|U_{\Delta t} - U_{\Delta t/2}|$')
legend('Location', 'best')
title('Sphere: temporal convergence')

