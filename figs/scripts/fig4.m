
addpath('../../utils')
addpath('../../src')

close all

fig = journal_figure([6.75 2.25], 2);

lspan = [0 0.01 0.1 1];
Lspan = linspace(0, 1, 100);
cmap = cmocean('algae', length(lspan) + 1);
cmap = cmap(2 : end, :);

cache = '../../data/fig4.mat';

if isfile(cache)
  load(cache, 'AsaveAll', 'Amax')
else
  hmax = pi/2;

  H0 = 1e-2;
  Hf = pi-H0;

  R = 10;
  N = 100;

  AsaveAll = cell(length(lspan), 1);

  for k = 1:length(lspan)
    lambda = lspan(k);

    F = @(H) darcySphere(H, lambda, R, N);
    Asave = [];

    % Main loop
    for H = linspace(H0, Hf, 100)
      dAdt = 2 * pi * R^2 * sin(H) * F(H);
      Asave = [Asave; H dAdt];
    end

    AsaveAll{k} = Asave;
  end

  Amax = [0 pi];

  for lambda = Lspan(2 : end)

    t = 0;

    H = H0;
    F = @(H) darcySphere(H, lambda, R, N);
    dAdt = @(H) -2 * pi * sin(H) * F(H);

    hmax = fminsearch(dAdt, hmax);

    Amax = [Amax; lambda hmax];

  end

  save(cache, 'AsaveAll', 'Amax')
end

sp1 = subplot(1, 2, 1);
for k = 1:length(lspan)
  Asave = AsaveAll{k};
  plot(Asave(:, 1), Asave(:, 2) / max(Asave(:, 2)), 'Color', cmap(k, :), 'DisplayName', strcat('$\lambda = ~$', num2str(lspan(k))))
  hold on
  plot([pi/2 pi/2], [0 1], 'k--', 'HandleVisibility', 'off')
end
xlim([0 pi])
ylim([0 1])
xlabel('polar angle, $H$')
ylabel('growth rate, $|\dot\Omega| / |\dot\Omega|_{\rm max}$')
xticks([0 pi/2 pi])
xticklabels(["$0$", "$\pi/2$", "$\pi$"])

leg = legend('location', 'southeast');
leg.Units = 'inches';

sp2 = subplot(1, 2, 2);
plot(Amax(:, 1), Amax(:, 2), 'Color', [0.2 0.2 0.6]), hold on
plot(Amax(:, 1), pi/2 + 0 * Amax(:, 1), 'k--')

ylim([0.9 * pi/2 pi])
yticks([pi/2 3*pi/4 pi])
yticklabels(["$\pi/2$", "$3\pi/4$", "$\pi$"])
xlabel('stress sensitivity, $\lambda$')
ylabel('angle of maximal growth, $H^*$')

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

leg.Position(1) = leg.Position(1) - 0.5;
