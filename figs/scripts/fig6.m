
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plots (a)-(b) boundary snapshots (equally spaced in time) for the isotropic
% and anisotropic moving-boundary simulations, and (c) domain area growth
% over time, normalized by the initial area.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
addpath('../../utils')
addpath('../../src')

Niso = load('../../data/fig6/isotropic.mat');
Nani = load('../../data/fig6/anisotropic.mat');

Nsnap = 5;

% ── Figure ────────────────────────────────────────────────────────────────────

fig = journal_figure([6.75 2.25], 2);
ax_a = subplot(1, 3, 1);
ax_b = subplot(1, 3, 2);
ax_c = subplot(1, 3, 3);

c_growth = [0.1 0.6 0.3];
cmap = cmocean('algae', Nsnap + 1);
cmap = cmap(2:end, :);

[Ra, tRangeA] = plotSnapshots(ax_a, Niso, Nsnap, cmap);
[Rb, tRangeB] = plotSnapshots(ax_b, Nani, Nsnap, cmap);

R = 1.05 * max(Ra, Rb);
xlim(ax_a, [-R R]); ylim(ax_a, [-R R]);
xlim(ax_b, [-R R]); ylim(ax_b, [-R R]);

addTimeColorbar(ax_a, cmap, tRangeA)
addTimeColorbar(ax_b, cmap, tRangeB)

axes(ax_c)
plot(Niso.AA(:,1), Niso.AA(:,2) / Niso.A0, 'Color', c_growth([1 3 2]), ...
    'DisplayName', '$h(\xi'',\eta'') = a\sin(2\pi \xi'')\sin(2\pi\eta'')$'), hold on
plot(Nani.AA(:,1), Nani.AA(:,2) / Nani.A0, 'Color', c_growth, ...
    'DisplayName', '$h(\xi'',\eta'') = a\sin(2\pi \xi'')$')
xlabel('time, $t$', 'Interpreter', 'latex')
ylabel('$|\Omega| / |\Omega_0|$', 'Interpreter', 'latex')
legend('Location', 'southeast', 'EdgeColor', 'none', 'Color', 'none')
ylim([0 10])

ax_a.Units = 'inches';
ax_b.Units = 'inches';
ax_c.Units = 'inches';

ax_a.Position(1) = -0.25;
ax_a.Position(2) = 0.25;
ax_a.Position(3) = 0.25 * fig.Position(3);
ax_a.Position(4) = fig.Position(4);

ax_b.Position(2:4) = ax_a.Position(2:4);

ax_c.Position(2) = 0.8;
ax_c.Position(3) = 0.35 * fig.Position(3);
ax_c.Position(4) = 0.75 * fig.Position(4);

ax_b.Position(1) = ax_a.Position(1) + ax_a.Position(3) + 0.5;
ax_c.Position(1) = ax_b.Position(1) + ax_b.Position(3) + 1.5;


tb_a = annotation('textbox', 'String', '(a)', ...
    'EdgeColor', 'none', 'Interpreter', 'latex', 'FontSize', ax_a.FontSize, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Units', 'inches');
tb_b = annotation('textbox', 'String', '(b)', ...
    'EdgeColor', 'none', 'Interpreter', 'latex', 'FontSize', ax_b.FontSize, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Units', 'inches');
tb_c = annotation('textbox', 'String', '(c)', ...
    'EdgeColor', 'none', 'Interpreter', 'latex', 'FontSize', ax_c.FontSize, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Units', 'inches');

tb_a.Position(1) = ax_a.Position(1)+0.5;
tb_a.Position(2) = ax_a.Position(2) + ax_a.Position(4) - 1;

tb_b.Position(1) = ax_b.Position(1)+0.5;
tb_b.Position(2) = tb_a.Position(2);
tb_c.Position(1) = ax_c.Position(1);
tb_c.Position(2) = tb_a.Position(2);

function [R, tRange] = plotSnapshots(ax, data, Nsnap, cmap)

  axes(ax)
  idx = round(linspace(1, size(data.AA, 1), Nsnap));
  tRange = data.AA(idx([1 end]), 1)';

  % Convert flat/homogenized coordinates to physical scale (area ~ J * flat area)
  s = sqrt(data.J);

  % Draw largest (latest) snapshot first so smaller/earlier ones stay visible
  R = 0;
  for k = Nsnap:-1:1
    xb = data.Xb{idx(k)};
    yb = data.Yb{idx(k)};
    xfull = s * vertcat(xb{:});
    yfull = s * vertcat(yb{:});
    fill(xfull, yfull, cmap(k, :), 'EdgeColor', 'none'), hold on
    R = max(R, max(sqrt(xfull.^2 + yfull.^2), [], 'all'));
  end
  axis equal off

end

function addTimeColorbar(ax, cmap, tRange)

  pos = ax.Position;
  colormap(ax, cmap)
  clim(ax, tRange)
  cb = colorbar(ax);
  ax.Position = pos;  % restore axes size

  cb.TickLabelInterpreter = 'latex';
  cb.Ticks = tRange;
  cb.TickLabels = {sprintf('%.0f', tRange(1)), sprintf('%.0f', tRange(2))};
  ylabel(cb, '$t$', 'Interpreter', 'latex', 'Rotation', 0, 'FontSize', ax.FontSize)

end
