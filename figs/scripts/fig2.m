
% Pressure and velocity profiles on a spherical cap for a range of radii R.

addpath('../../utils')
addpath('../../src')

c_pressure = [0.1 0.1 0.6];
c_velocity = [0.6 0.1 0.6];

N = 100;

fig = journalFig([6.75 2.5], 2);

Rplt = 10 + [0, 2.^(1:5)];
lambda = 1;

cmap_pressure = 1 - linearrgbmap(1 - c_pressure, length(Rplt)+1);
cmap_pressure(1, :) = [];

cmap_velocity = 1 - linearrgbmap(1 - c_velocity, length(Rplt)+1);
cmap_velocity(1, :) = [];

% ── Solve ────────────────────────────────────────────────────────────────────

cache = '../../data/fig2.mat';

if isfile(cache)
  load(cache, 'hAll', 'pAll', 'uAll')
else
  H = pi / 2;

  hAll = cell(length(Rplt), 1);
  pAll = cell(length(Rplt), 1);
  uAll = cell(length(Rplt), 1);

  for k = 1:length(Rplt)
    R = Rplt(k);
    [~, p, h, u] = growthFlowSphere(H, lambda, R, N);
    hAll{k} = h;
    pAll{k} = p;
    uAll{k} = u;
  end

  save(cache, 'hAll', 'pAll', 'uAll')
end

% ── Figure ───────────────────────────────────────────────────────────────────

nc = 1;
for k = 1:length(Rplt)

  h = hAll{k};
  p = pAll{k};
  u = uAll{k};

  subplot(1, 2, 1)
  plot(h, p, 'Color', cmap_pressure(nc, :), 'DisplayName', sprintf('$R$ = %2d', Rplt(k))), hold on
  plot(-flip(h), flip(p), 'Color', cmap_pressure(nc, :), 'HandleVisibility', 'off'), hold on

  subplot(1, 2, 2)
  plot(h, u, 'Color', cmap_velocity(nc, :), 'DisplayName', sprintf('$R$ = %2d', Rplt(k))), hold on
  plot(-flip(h), -flip(u), 'Color', cmap_velocity(nc, :), 'HandleVisibility', 'off'), hold on

  nc = nc + 1;

end

sp1 = subplot(1, 2, 1);
xlabel("polar angle, $\eta$")
ylabel("pressure, $p$")
xlim(1 * [-pi/2 pi/2])
xticks([-pi/2 0 pi/2])
xticklabels(["$-\pi/2$", "$0$", "$\pi/2$"])
legend('location', 'south')

sp2 = subplot(1, 2, 2);
xlabel("polar angle, $\eta$")
ylabel("velocity, $u$")
xlim(1 * [-pi/2 pi/2])
xticks([-pi/2 0 pi/2])
xticklabels(["$-\pi/2$", "$0$", "$\pi/2$"])
legend('location', 'south')

sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position(1) = 0.85;
sp1.Position(2) = 0.7;
sp1.Position(3) = 0.4 * fig.PaperSize(1);
sp1.Position(4) = 0.8 * fig.PaperSize(2);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.2;
sp2.Position(2 : 4) = sp1.Position(2 : 4);

ax = gca;

tba = annotation('textbox', 'string', '(a)', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tba.Position(1) = sp1.Position(1);
tba.Position(2) = sum(sp1.Position([2 4])) - tba.Position(4);

tbb = annotation('textbox', 'string', '(b)', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tbb.Position(1) = sp2.Position(1);
tbb.Position(2) = sum(sp2.Position([2 4])) - tbb.Position(4);

tba.Position(2) = tbb.Position(2);
