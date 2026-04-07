
close all

fig = journal_figure([6.5 2.45], 2);

source = '../../data/homogenization-evolution';

sp1 = subplot(1, 3, 1);
render(strcat(source,'/frames-eps4-lambda1'))

sp2 = subplot(1, 3, 2);
render(strcat(source,'/frames-eps2-lambda1'))
ylabel('')
yticklabels('')
sp3 = subplot(1, 3, 3);
render(strcat(source,'/frames-eps1-lambda1'))
ylabel('')
yticklabels('')

cmap = cmocean('algae');
colormap(cmap)
clim([0 1])
hb = colorbar;

sp1.Units = 'inches';
sp2.Units = 'inches';
sp3.Units = 'inches';

sp1.Position(1 : 2) = [0.6 0.7];
sp1.Position(3 : 4) = fig.PaperSize(2) * 0.75;

sp2.Position(2 : 4) = sp1.Position(2 : 4);
sp3.Position(2 : 4) = sp1.Position(2 : 4);

sp2.Position(1) = sum(sp1.Position([1 3])) + 0.3;
sp3.Position(1) = sum(sp2.Position([1 3])) + 0.3;

ax = gca;

hb.TickLabelInterpreter = 'latex';
hb.Ticks = [0 1];
hb.FontSize = ax.FontSize;
hb.TickLabels = {'0', '$t_f$'};
ylabel(hb, '$t$', 'interpreter', 'latex', 'rotation', 0, 'FontSize', ax.FontSize)


tb = annotation('textbox', 'string', '(a) $\varepsilon = 4$', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tb.Position(1) = sp1.Position(1);
tb.Position(2) = sum(sp1.Position([2 4])) - tb.Position(4);

tb = annotation('textbox', 'string', '(b) $\varepsilon = 2$', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tb.Position(1) = sp2.Position(1);
tb.Position(2) = sum(sp2.Position([2 4])) - tb.Position(4);

tb = annotation('textbox', 'string', '(c) $\varepsilon = 1$', 'interpreter', 'latex', 'FontSize', ax.FontSize', 'EdgeColor', 'none', 'Units', 'inches');
tb.Position(1) = sp3.Position(1);
tb.Position(2) = sum(sp3.Position([2 4])) - tb.Position(4);

print('../output/fig6', '-dpdf')

function render(src)

  Nf = length(dir(strcat(src, '/*.mat')));
  cmap = cmocean('algae', Nf + 1);
  Nspan = flip(unique(floor(linspace(2, Nf, 8))));
  
  for nf = Nspan
    
    load(strcat(src, sprintf('/f-%d.mat', nf)), 'xb', 'yb', 't')
    x = [];
    y = [];
  
    for mb = 1 : length(xb)
      x = [x; xb{mb}];
      y = [y; yb{mb}];
    end
  
    if nf == Nf
      ec = [0 0 0];
    else
      ec = cmap(nf + 1, :);
    end
  
    fill(x, y, cmap(nf + 1, :), 'LineWidth', 2, 'EdgeColor', ec), hold on
    
  end
  
  
  % s = linspace(0, 2 * pi);
  % fill(cos(s), sin(s), cmap(1, :), 'LineWidth', 2, 'EdgeColor', cmap(1, :)), hold on
  
  xlabel('$x$'), ylabel('$y$')
  axis equal
  axis(2.2 * [-1 1 -1 1])

end