
function u = rk4(u, f, dt, Nt)

  for nt = 1 : Nt
    k1 = f(u);
    k2 = f(u + dt * k1 / 2);
    k3 = f(u + dt * k2 / 2);
    k4 = f(u + dt * k3);
    u = u + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
  end

end