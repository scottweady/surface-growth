
% Plots (a) growth rate vs roughness amplitude and (b)-(d) flow snapshots
% for flat, 1D-wave, and 2D-wave surfaces from the homogenized growth flow
% equations.

% ── Load ─────────────────────────────────────────────────────────────────────
% Solved by scripts/roughnessGrowthRate.m, run via data/makedata.m

load('../../data/fig5.mat', 'hom', 'X', 'Y', 'U', 'V', 'G')

gr_lim = [min(cellfun(@(g) min(g(:)), G)), 1];

% ── Figure ────────────────────────────────────────────────────────────────────

fig = journalFig([6.75 2], 2);

ml = 0.08;  mr = 0.08;  mb = 0.18;  mt = 0.06;
gap = 0.015;
h_ax = 1 - mb - mt;
wa   = 0.25;
wb   = (1 - ml - mr - wa - 3*gap) / 3;

ax_a = axes('Position', [ml, mb, wa, h_ax]);

c_growth = [0.1 0.6 0.3];
plot(hom(:,1), hom(1,3) * ones(size(hom(:,1))), '-',  'Color', 'k', ...
     'DisplayName', '$h = 0$'), hold on
plot(hom(:,1), hom(:,2), 'o-.', 'Color', c_growth, ...
     'DisplayName', '$h = a\sin(2\pi\xi'')$')
plot(hom(:,1), hom(:,3), 's:',  'Color', c_growth([1 3 2]), ...
     'DisplayName', '$h = a\sin(2\pi\xi'')\sin(2\pi\eta'')$')
xlabel('roughness amplitude, $a$')
ylabel('growth rate, $|\dot\Omega|/|\Omega|$')
legend('Location', 'southwest', 'EdgeColor', 'none')
ylim([0.86 0.96])

% Panels (b)-(d)
ng = 80;   % uniform grid resolution for streamline interpolation
Ns = 24;   % seeds per panel

ax_snap = gobjects(3,1);
for k = 1:3
  xl = ml + wa + k*gap + (k-1)*wb;
  ax = axes('Position', [xl, mb, wb, h_ax]);
  ax_snap(k) = ax;
  pcolor(X{k}, Y{k}, G{k}); shading interp; axis equal off; hold on

  colormap(ax, cmocean('algae'))
  clim([0.8 1.0]);

  R  = max(abs(X{k}(:)));
  xg = linspace(-R, R, ng);
  yg = linspace(-R, R, ng);
  [Xg, Yg] = meshgrid(xg, yg);

  warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
  Fu = scatteredInterpolant(X{k}(:), Y{k}(:), U{k}(:));
  Fv = scatteredInterpolant(X{k}(:), Y{k}(:), V{k}(:));
  warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

  Ug = Fu(Xg, Yg);
  Vg = Fv(Xg, Yg);
  Ug(Xg.^2 + Yg.^2 > R^2) = NaN;
  Vg(Xg.^2 + Yg.^2 > R^2) = NaN;

  s0 = 1e-3;
  if k == 2
    xs = [linspace(-R, R, Ns), linspace(-R, R, Ns)];
    ys = [s0*ones(1,Ns), -s0*ones(1,Ns)];
  else
    th = linspace(0, 2*pi, Ns+1); th = th(1:end-1);
    xs = s0*cos(th);
    ys = s0*sin(th);
  end

  verts = stream2(xg, yg, Ug, Vg, xs, ys);
  hw = 0.02*R;
  hl = 0.04*R;
  for i = 1:length(verts)
    v = verts{i};
    v = v(~any(isnan(v), 2), :);
    if size(v,1) < 4, continue; end
    plot(v(:,1), v(:,2), 'Color', 0.99*[1 1 1], 'LineWidth', 0.6)
    n = size(v,1);
    for frac = 0.2:0.2:1
      mid  = max(1, min(round(frac*n), n-1));
      d    = v(mid+1,:) - v(mid,:);
      if norm(d) == 0, continue; end
      d    = d / norm(d);
      perp = [-d(2), d(1)];
      tip   = v(mid,:) + d    * hl;
      base1 = v(mid,:) + perp * hw;
      base2 = v(mid,:) - perp * hw;
      patch([tip(1) base1(1) base2(1)], [tip(2) base1(2) base2(2)], ...
            0.99*[1 1 1], 'EdgeColor', 'none')
    end
  end

  if k == 3
    pos_d = ax.Position;
    cb = colorbar; cb.Label.Interpreter = 'latex';
    cb.Label.String = 'growth rate, $\gamma(p)$';
    cb.Ticks = linspace(0.8, 1, 5);
    cb.TickLabelInterpreter = 'latex';
    ax.Position = pos_d;  % restore axes size
    cb_h = 0.8 * pos_d(4);
    cb.Position = [pos_d(1)+pos_d(3)+0.005, pos_d(2)+0.5*(pos_d(4)-cb_h), 0.01, cb_h];
  end

end

% Panel labels at a shared vertical position
axs  = [{ax_a}, num2cell(ax_snap)'];

ytop = max(cellfun(@(ax) ax.Position(2) + ax.Position(4), axs)) - 0.05;
for k = 1:4
  pos = axs{k}.Position;
  lbl = char('a' + k - 1);
  tb = annotation('textbox', [pos(1) ytop 0.05 0.05], 'String', ['(' lbl ')'], ...
      'EdgeColor', 'none', 'Interpreter', 'latex', 'FontSize', axs{k}.FontSize, ...
      'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
if k == 1
    tb.Position(1) = 0;
end
end

% Height function labels at the bottom of panels (b)-(d)
h_labels = {'$h = 0$', ...
            '$h = a\sin(2\pi\xi'')$', ...
            '$h = a\sin(2\pi\xi'')\sin(2\pi\eta'')$'};
for k = 1:3
  ax = axs{k+1};
  text(ax, 0.5, -0.05, h_labels{k}, 'Units', 'normalized', ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
      'Interpreter', 'latex', 'FontSize', ax.FontSize, 'Clipping', 'off')
end
