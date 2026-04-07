
addpath('utils')

c_growth = [0.1 0.6 0.3];
c_pressure = [0.1 0.1 0.6];
c_velocity = [0.6 0.1 0.6];

%%
close all

fig = journal_figure([6.5 2.5], 2);

% Range of sphere radii
Rplt = 10 + [0, 2.^(1:5)];
Rspan = Rplt;

% Colormap
cmap = 1 - linearrgbmap(1 - c_growth, length(Rplt) + 1);
cmap = cmap(2 : end, :);

% Color counter
nc = 1;

% Final time 
ts = [];

% Resolution
N = 100;

% Sensitivity
lambda = 1;
  
for R = Rspan

  clc, fprintf("%f of %f\n", R, Rspan(end))

  
  % Time step
  dt = 0.1;
  
  % Final time
  tf = 100;
  
  % Time
  t = 0;
  
  % Initial area
  A0 = 4 * pi;
  
  % Initial angle
  H0 = acos(1 - A0 / (2 * pi * R^2));
  H = H0;
  
  % Final location
  Af = 4 * pi * 10^2;
  Hf = acos(1 - Af / (2 * pi * R^2));
  
  % Area
  A = 2 * pi * R^2 * (1 - cos(H0));
  AA = [0, A];

  % RHS function
  F = @(Theta) darcy_sphere(Theta, lambda, R, N);
  
  Hsave = [0, H];

  % Main loop
  while t < tf && H < Hf
  
    % Update interface location
    H = rk4(H, F, dt, 1);
    t = t + dt;

    Hsave = [Hsave; t, H];

    % Area
    A = 2 * pi * R^2 * (1 - cos(H));
    AA = [AA; t, A];  

  end
  
  if R == 10
    Hplt = Hsave;
  end

  % Reverse at high res to find final time
  while H > Hf
    H = rk4(H, F, -0.01*dt, 1);
    t = t - 0.01*dt;
  end

  % Get final time
  ts = [ts; R, t];

  subplot(1, 2, 2)
  if any(R == Rplt)
    plot(AA(:, 1), AA(:, 2)/(4*pi), 'DisplayName', sprintf('$R$ = %2d', R), 'Color', cmap(nc, :)), hold on
    nc = nc + 1;
  end
  
end

t = 0;
R = sqrt(A0 / pi);
A = pi * R^2;
AA = [0, A];

while A < Af
  % Flat interface
  c = @(R) (sqrt(1 + lambda * R^2 / 2) - 1) / (2 * lambda * R^2);
  F = @(R) 32 * c(R)^2 * R / (1 - 8 * c(R)^2 * lambda * R^2);
  R = rk4(R, F, dt, 1);
  t = t + dt;
  A = pi * R^2;
  AA = [AA; t, A];
  
end

% Reverse at high res to find final time
while A > Af
  R = rk4(R, F, -0.01*dt, 1);
  A = pi * R^2;
  t = t - 0.01*dt;
end
ts = [ts; inf t];

sp2 = subplot(1, 2, 2);
plot(AA(:, 1), AA(:, 2) / (4 * pi), 'k--', 'DisplayName', '$R \rightarrow \infty$'), hold on

legend('location', 'southeast')
xlabel('time, $t$')
ylabel('surface area, $|\Omega| / |\Omega_0|$')
ylim([0 Af/(4*pi)])

sp1 = subplot(1, 2, 1);

nspan = round(linspace(1, length(Hplt)-1, 8));
cmap = flipud(cmocean('algae', length(nspan)+1));

hm1 = 0;
for n = 1 : length(nspan)
  xi = linspace(0, 2 * pi, 256);
  h = Hplt(nspan(n), 2);
  eta = linspace(hm1, h, 256);
  [eta, xi] = meshgrid(eta, xi);
  surf(10 * cos(xi) .* sin(eta), 10 * sin(xi) .* sin(eta), 10 * cos(eta), 'FaceColor', cmap(n, :), 'EdgeColor', 'none'), hold on

  xi = linspace(0, 2 * pi, 256);
  h = Hplt(nspan(n), 2);
  a = 10;
  plot3(a * cos(xi) .* sin(h), a * sin(xi) .* sin(h), a * cos(h) + 0 * xi, 'Color', cmap(n, :).^1.5), hold on
  hm1 = h;
end

% Colormap
colormap(cmap)
hb = colorbar;
clim([0 1])

% Axis format
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

print('../output/fig3', '-dpdf')


%%
close all
fig = journal_figure([6.5 2.5], 2);

Rplt = 10 + [0, 2.^(1:5)];
lambda = 1;

cmap_pressure = 1 - linearrgbmap(1 - c_pressure, length(Rplt)+1);
cmap_pressure(1, :) = [];

cmap_velocity = 1 - linearrgbmap(1 - c_velocity, length(Rplt)+1);
cmap_velocity(1, :) = [];

nc = 1;
H = pi / 2;

for R = Rplt

  [U, p, h, u] = darcy_sphere(H, lambda, R, N);

  subplot(1, 2, 1)
  plot(h, p, 'Color', cmap_pressure(nc, :), 'DisplayName', sprintf('$R$ = %2d', R)), hold on
  plot(-flip(h), flip(p), 'Color', cmap_pressure(nc, :), 'HandleVisibility', 'off'), hold on

  subplot(1, 2, 2)
  plot(h, u, 'Color', cmap_velocity(nc, :), 'DisplayName', sprintf('$R$ = %2d', R)), hold on
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

print('../output/fig2', '-dpdf')

%%
close all

fig = journal_figure([6.5 2.25], 2);

lspan = [0 0.01 0.1 1];

Lspan = linspace(0, 1, 100);
cmap = cmocean('algae', length(lspan) + 1);
cmap = cmap(2 : end, :);
nl = 1;

hmax = pi/2;

H0 = 1e-2;
Hf = pi-H0;

R = 10;
N = 100;

for lambda = lspan

  F = @(H) darcy_sphere(H, lambda, R, N);
  Asave = [];

  % Main loop
  for H = linspace(H0, Hf, 100)
    dAdt = 2 * pi * R^2 * sin(H) * F(H);
    Asave = [Asave; H dAdt];
  end
  
  sp1 = subplot(1, 2, 1);
  plot(Asave(:, 1), Asave(:, 2) / max(Asave(:, 2)), 'Color', cmap(nl, :), 'DisplayName', strcat('$\lambda = ~$', num2str(lambda)))
  nl = nl + 1;
  hold on
  plot([pi/2 pi/2], [0 1], 'k--', 'HandleVisibility', 'off')
  xlim([0 pi])
  ylim([0 1])
  xlabel('polar angle, $H$')
  ylabel('growth rate, $|\dot\Omega| / |\dot\Omega|_{\rm max}$')
  xticks([0 pi/2 pi])
  xticklabels(["$0$", "$\pi/2$", "$\pi$"])
  
end

leg = legend('location', 'southeast');
leg.Units = 'inches';

Amax = [0 pi];

for lambda = Lspan(2 : end)

  t = 0;
  
  H = H0;
  F = @(H) darcy_sphere(H, lambda, R, N);
  dAdt = @(H) -2 * pi * sin(H) * F(H);
  
  hmax = fminsearch(dAdt, hmax);
  
  Amax = [Amax; lambda hmax];
  
end

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

print('../output/fig4', '-dpdf')

