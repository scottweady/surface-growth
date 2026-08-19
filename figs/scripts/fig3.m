
% Growth trajectories on spherical caps for a range of radii R, compared
% against the flat (R -> infinity) limit.

c_growth = [0.1 0.6 0.3];

N = 100;

fig = journalFig([6.75 2.5], 2);

Rplt = 10 + [0, 2.^(1:5)];
Rspan = Rplt;

cmap = 1 - linearrgbmap(1 - c_growth, length(Rplt) + 1);
cmap = cmap(2 : end, :);

lambda = 1;

A0 = 4 * pi;
Af = 4 * pi * 10^2;

% ── Load ─────────────────────────────────────────────────────────────────────
% Solved by scripts/sphereGrowth.m, run via data/makedata.m

load('../../data/fig3.mat', 'Hplt', 'ts', 'AAall', 'AAinf')

% ── Figure ───────────────────────────────────────────────────────────────────

sp2 = subplot(1, 2, 2);
for k = 1:length(Rplt)
  plot(AAall{k}(:, 1), AAall{k}(:, 2)/(4*pi), 'DisplayName', sprintf('$R$ = %2d', Rplt(k)), 'Color', cmap(k, :)), hold on
end
plot(AAinf(:, 1), AAinf(:, 2) / (4 * pi), 'k--', 'DisplayName', '$R \rightarrow \infty$'), hold on

legend('location', 'southeast')
xlabel('time, $t$')
ylabel('surface area, $|\Omega| / |\Omega_0|$')
ylim([0 Af/(4*pi)])

sp1 = subplot(1, 2, 1);

nspan = round(linspace(1, length(Hplt)-1, 8));
cmap = cmocean('algae', length(nspan)+1);
cmap = cmap(2:end,:);

hm1 = 0;
for n = 1 : length(nspan)
  xi = linspace(0, 2 * pi, 32);
  h = Hplt(nspan(n), 2);
  eta = linspace(hm1, h, 32);
  [eta, xi] = meshgrid(eta, xi);
  surf(10 * cos(xi) .* sin(eta), 10 * sin(xi) .* sin(eta), 10 * cos(eta), 'FaceColor', cmap(n, :), 'EdgeColor', 'none'), hold on

  xi = linspace(0, 2 * pi, 32);
  h = Hplt(nspan(n), 2);
  a = 10;
  plot3(a * cos(xi) .* sin(h), a * sin(xi) .* sin(h), a * cos(h) + 0 * xi, 'Color', cmap(n, :).^1.5), hold on
  hm1 = h;
end

colormap(cmap)
hb = colorbar;
clim([0 1])

axis equal tight off
box on

ax = gca;
hb.TickLabelInterpreter = 'latex';
hb.Ticks = [0 1];
hb.FontSize = ax.FontSize;
hb.TickLabels = {'0', '$t_f$'};
ylabel(hb, '$t$', 'interpreter', 'latex', 'rotation', 0, 'FontSize', ax.FontSize)

sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position([1 2]) = [-0.25 -0.25];
sp1.Position([3 4]) = 1.1 * fig.PaperSize(2);

hb.Units = 'inches';

hb.Position(1) = sum(sp1.Position([1 3]))-0.75;
hb.Position(4) = 0.5 * sp1.Position(4);
hb.Position(2) = sp1.Position(2) + 0.5 * hb.Position(4);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.5;
sp2.Position(2) = 0.75;
sp2.Position([3 4]) = [0.4 * fig.PaperSize(1) 0.75 * fig.PaperSize(2)];

tba = annotation('textbox', 'string', '(a)', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tba.Position(1) = sp1.Position(1);
tba.Position(2) = sum(sp1.Position([2 4])) - tba.Position(4);

tbb = annotation('textbox', 'string', '(b)', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tbb.Position(1) = sp2.Position(1);
tbb.Position(2) = sum(sp2.Position([2 4])) - tbb.Position(4);

tba.Position(2) = tbb.Position(2);
tba.Position(1) = 0.6;
