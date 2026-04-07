
close all

lambda = 1;
D = surfacemesh.disk(12, 4);
[~, p0] = darcy(D, lambda);
G0 = integral(exp(-p0));

mspan = 1 : 5;
G = zeros(length(mspan), 1);

for nm = 1 : length(mspan)

  m = mspan(nm);
  Omega = surfacemesh.star(12, a, m, 8);  
  [~, p] = darcy(Omega, lambda);
  G(nm) = integral(exp(-p));

  cla
  plot(p)
  drawnow

end

plot(mspan, G, 'o-')
hold on
plot(mspan, G0 + 0 * mspan, '--')
plot(mspan, pi + 0 * mspan, '--')
