
addpath('../utils')
addpath('../src')
addpath('../scripts')

lambda = 1;
R = 1;

%% Spatial convergence

close all
UU = [];
NN = 4 : 2 : 32;

for N = NN
  fprintf('N = %d\n', N)
  UU = [UU; darcySphere(pi/2, lambda, R, N)];
end

err = abs(diff(UU));

semilogy(NN(2 : end), err, 'ko-')

xlim([0 max(NN)+(NN(end)-NN(end-1))])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('resolution, $N$')
ylabel('error, $|U_N - U_{N/2}|$')

print('output/fig-sphere-spatial-convergence', '-dpdf')

%% Temporal convergence

close all
HH = [];

N = 100;
NNt = 2.^(1 : 8);
F = @(Theta) darcySphere(Theta, lambda, R, N);

for Nt = NNt
  dt = 1 / Nt;
  fprintf('dt = %g\n', dt)
  HH = [HH; rk4(pi/4, F, dt, Nt)];
end

dt = 1 ./ NNt(2 : end);

err = abs(diff(HH));

loglog(dt, err, 'ko-'), hold on
loglog(dt, 10 * err(1) * (dt/dt(1)).^(4), '--', 'Color', [0.2 0.2 0.2])

xlim([0.5 * min(dt) 2 * max(dt)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('time step size, $\Delta t$')
ylabel('error, $|U_{\Delta t} - U_{\Delta t/2}|$')

print('output/fig-sphere-temporal-convergence', '-dpdf')

