% HERO  README hero image: a colony spreading over a sphere (as in Fig. 3a),
% with fronts at equally spaced times shaded from early (light) to late (dark).
%
% Saves a transparent PNG to figs/hero.png. Run from figs/scripts.

load('../../data/fig3.mat', 'Hplt')

R = 10;
nspan = round(linspace(1, length(Hplt)-1, 8));
cmap = cmocean('algae', length(nspan)+1);
cmap = cmap(2:end, :);

fig = figure('Units', 'inches', 'Position', [2 2 4 4], 'InvertHardcopy', 'off');
axes(fig, 'Position', [0 0 1 1]);

xi = linspace(0, 2 * pi, 256);
hm1 = 0;
for n = 1:length(nspan)
  h = Hplt(nspan(n), 2);
  [eta, x] = meshgrid(linspace(hm1, h, 32), xi);
  surf(R * cos(x) .* sin(eta), R * sin(x) .* sin(eta), R * cos(eta), ...
       'FaceColor', cmap(n, :), 'EdgeColor', 'none'), hold on
  plot3(R * cos(xi) * sin(h), R * sin(xi) * sin(h), R * cos(h) + 0 * xi, ...
        'Color', cmap(n, :).^1.5, 'LineWidth', 1.5)
  hm1 = h;
end

axis equal tight off

% ── Transparent export ───────────────────────────────────────────────────────
% Render on white and on black; the difference gives the (anti-aliased)
% alpha channel exactly, and the black render gives premultiplied color.

dpi = 300;
fig.Color = 'w';
Iw = double(print(fig, '-RGBImage', sprintf('-r%d', dpi))) / 255;
fig.Color = 'k';
Ib = double(print(fig, '-RGBImage', sprintf('-r%d', dpi))) / 255;

alpha = min(max(1 - mean(Iw - Ib, 3), 0), 1);
rgb   = min(Ib ./ max(alpha, eps), 1);

rows = find(any(alpha > 0, 2));
cols = find(any(alpha > 0, 1));
rows = rows(1):rows(end);
cols = cols(1):cols(end);

imwrite(rgb(rows, cols, :), '../hero.png', 'Alpha', alpha(rows, cols))
