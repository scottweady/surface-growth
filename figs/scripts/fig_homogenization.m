
close all

fig = journal_figure([6.5 2.25], 2);

eps = 0.5;

% Surface parameterization
k = 2 * pi / eps;
a = 0.5;
h = @(x, y) a * eps * sin(k * x) .* sin(k * y);

D = surfacemesh.disk(12, 4);

% Construct graph
x = D.x;
y = D.y;
z = cell(1, length(D));

s1 = linspace(-1, 1, 10);
s2 = linspace(-1, 1, 10);
[s2, s1] = meshgrid(s2, s1);

sp1 = subplot(1, 2, 1);

for m = 1 : length(D)
  z{m} = h(x{m}, y{m});
  surf(x{m}, y{m}, z{m}), hold on
  mesh(s1, s2, 0*s2, 'FaceAlpha', 0.25, 'EdgeColor', 'k')
end

colormap(flipud(cmocean('algae')))
shading interp

axis equal
xticks([-1 0 1])
yticks([-1 0 1])
axis(1.1 * [-1 1 -1 1 -a * eps a * eps])
grid off
box on

xlabel('$x$')
ylabel('$y$')
zlabel('$z$')

sp2 = subplot(1, 2, 2);

c_growth = [0.1 0.6 0.3];
cmap = 1 - linearrgbmap(1 - c_growth, 6);
cmap = cmap(2 : end, :);

source = '../../data/homogenization';
load(strcat(source, '/', 'eps1.mat'))
plot(dns(:, 1), dns(:, 2), 'o--', 'Color', cmap(1, :), 'DisplayName', 'DNS, $\epsilon=1$'), hold on

load(strcat(source, '/', 'eps0.5.mat'))
plot(dns(:, 1), dns(:, 2), 'o--', 'Color', cmap(2, :), 'DisplayName', 'DNS, $\epsilon=0.5$'), hold on

load(strcat(source, '/', 'eps0.25.mat'))
plot(dns(:, 1), dns(:, 2), 'o--', 'Color', cmap(3, :), 'DisplayName', 'DNS, $\epsilon=0.25$')

load(strcat(source, '/', 'eps0.125.mat'))
plot(dns(:, 1), dns(:, 2), 'o--', 'Color', cmap(4, :), 'DisplayName', 'DNS, $\epsilon=0.125$')

load(strcat(source, '/', 'homogenization.mat'))
plot(hom(:, 1), hom(:, 2), 'Color', cmap(4, :), 'DisplayName', 'homogenization, $\epsilon\rightarrow0$')

xlabel('roughness amplitude, $a$')
ylabel('growth rate, $|\dot\Omega| / |\Omega|$')
legend('location', 'southwest', 'color', 'none')


sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position(1) = 0.5;
sp1.Position(2) = 0.75;

sp1.Position(3) = 0.45 * fig.PaperSize(1);
sp2.Position(3) = 0.35 * fig.PaperSize(1);

sp1.Position(4) = 0.75 * fig.PaperSize(2);
sp2.Position(4) = 0.7 * fig.PaperSize(2);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.25;
sp2.Position(2) = 0.5 * (fig.PaperSize(2) - sp2.Position(4)) + 0.25;

print('../output/fig5', '-dpdf')