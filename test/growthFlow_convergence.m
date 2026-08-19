
addpath('../utils')
addpath('../src')
addpath('../scripts')

%% Temporal convergence (RK4 moving-boundary stepping)
% Runs the growthFlowMovingBoundary.m script itself (shared workspace) at
% decreasing time step sizes over a short, fixed time window.

lambda = 1;
K      = [0.1 0.0; 0.0 1.0];

doPlot = false;
doSave = false;

T    = 0.2;
NNt  = 2.^(0 : 4);
Afin = [];

for k = 1:length(NNt)
  Nsteps = NNt(k);
  dt = T / Nsteps;
  tf = T;
  clc
  fprintf('Moving boundary solve temporal convergence...\n')
  fprintf('%s dt = %g\n', progressBar(k / length(NNt)), dt)
  growthFlowMovingBoundary
  Afin = [Afin; A];
end

err = abs(diff(Afin));
dt  = T ./ NNt(2 : end);

%% Plot temporal convergence

journal_figure([3.75 2.5], 2);

loglog(dt, err, 'ko-', 'DisplayName', 'error'), hold on
loglog(dt, 10 * err(1) * (dt/dt(1)).^(4), '--', 'Color', [0.2 0.2 0.2], 'DisplayName', '$O(\Delta t^4)$')

xlim([0.5 * min(dt) 2 * max(dt)])
ylim([0.1 * min(err) 10 * max(err)])

xlabel('time step size, $\Delta t$')
ylabel('error, $||\Omega|_{\Delta t} - |\Omega|_{\Delta t/2}|$')
legend('Location', 'best')

print('output/fig-anisotropic-temporal-convergence', '-dpdf')