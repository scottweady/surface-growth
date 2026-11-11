%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Symbolic verification of the linear stability of the axisymmetric (disk)
% solution on the plane (Appendix A).
%
%   -Lap p = exp(-lambda p)  in Omega,   p = 0 on dOmega,
%   V_n = -dp/dn             on dOmega.
%
% The boundary is r = R(theta, t). Two perturbation ansatzes are compared:
%   (abs) R = H(t) + eps(t) cos(m theta)       (absolute amplitude)
%   (rel) R = H(t) (1 + eps(t) cos(m theta))   (relative amplitude)
% and in each case d(eps)/dt = sigma eps is derived from scratch.
%
% Expected results, with c = c(lambda H^2):
%   sigma_abs = -16 c^2 (m^2 - 1) / (1 + 4 c (m - 1))
%   sigma_rel = sigma_abs - 4 c
%             = -4 c [1 + 4 c (m - 1)(m + 2)] / (1 + 4 c (m - 1))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

syms r theta H c lam positive
syms m integer positive
syms s epsl Hdot epsdot real

% Lambda = lambda H^2 enters only via c. Impose p0(H) = 0, i.e.
% 8 lambda H^2 c^2 + 8 c - 1 = 0, by eliminating lambda.
lamH = (1 - 8*c) / (8 * c^2 * H^2);

% ── 1. Base state from the Liouville solution ────────────────────────────────
% p = (1/lambda) log[ lambda (1 - |f|^2)^2 / (8 |f'|^2) ],  z = r e^{i theta}

a  = 8 * lam * c^2;
zk = @(k) r^k * exp( 1i*k*theta);   % z^k
zbk = @(k) r^k * exp(-1i*k*theta);  % conj(z)^k

liouville = @(F2, Fp2) (1/lam) * log(lam * (1 - F2)^2 / (8 * Fp2));

% Perturbed map f = sqrt(a) z + s z^(m+1); conjugates written explicitly
f   = sqrt(a) * zk(1)  + s * zk(m+1);
fb  = sqrt(a) * zbk(1) + s * zbk(m+1);
fp  = sqrt(a) + s * (m+1) * zk(m);
fpb = sqrt(a) + s * (m+1) * zbk(m);

pLiou = liouville(f*fb, fp*fpb);

p0 = simplify(subs(pLiou, s, 0));
fprintf('p0(r) = %s\n', char(p0))

lap = @(q) diff(q, r, 2) + diff(q, r)/r + diff(q, theta, 2)/r^2;

res0 = simplify(-lap(p0) - exp(-lam*p0));
fprintf('Base PDE residual:              %s\n', char(res0))
fprintf('p0(H) with lambda eliminated:   %s\n', char(simplify(subs(subs(p0, r, H), lam, lamH))))
fprintf('p0''(H)  = %s\n', char(simplify(subs(diff(p0, r), {r, lam}, {H, lamH}))))
fprintf('p0''''(H) = %s\n', char(simplify(subs(diff(p0, r, 2), {r, lam}, {H, lamH}))))

% ── 2. O(eps) pressure ───────────────────────────────────────────────────────
% d/ds at s = 0 of the Liouville solution; strip cos(m theta) to get the
% radial profile and compare with the form stated in Appendix A.

dps = subs(diff(pLiou, s), s, 0);
dps = simplify(rewrite(dps, 'cos'));
dps = simplify(subs(dps, exp(1i*m*theta), cos(m*theta) + 1i*sin(m*theta)));

p1hat = r^m * ((m+1) - (m-1)*a*r^2) / (1 - a*r^2);  % Appendix A, A = 1
ratio = simplify(dps / (p1hat * cos(m*theta)));
fprintf('d p / d s |_{s=0} / [p1hat cos(m theta)] = %s  (r, theta independent)\n', ...
        char(ratio))

% Independent check: p1hat cos(m theta) solves the linearized PDE
%   -Lap P1 = -lambda exp(-lambda p0) P1
P1 = p1hat * cos(m*theta);
res1 = simplify(-lap(P1) + lam*exp(-lam*p0)*P1);
fprintf('Linearized PDE residual:        %s\n', char(res1))

% ── 3. Linearize the boundary conditions ─────────────────────────────────────
% Exact kinematic condition for r = R(theta, t):
%   R_t = -( p_r - R_theta p_theta / R^2 ) |_{r = R}

syms A real
p = p0 + epsl * A * p1hat * cos(m*theta);

sigma_abs_expected = -16*c^2*(m^2 - 1) / (1 + 4*c*(m - 1));
sigma_rel_expected = -4*c*(1 + 4*c*(m - 1)*(m + 2)) / (1 + 4*c*(m - 1));

cases = {'abs', H + epsl*cos(m*theta), Hdot + epsdot*cos(m*theta), ...
                sigma_abs_expected;
         'rel', H*(1 + epsl*cos(m*theta)), ...
                Hdot*(1 + epsl*cos(m*theta)) + H*epsdot*cos(m*theta), ...
                sigma_rel_expected};

sigma = struct();
for ic = 1 : size(cases, 1)

  name = cases{ic, 1};
  R    = cases{ic, 2};
  Rt   = cases{ic, 3};
  fprintf('\n── Ansatz (%s): R = %s\n', name, char(R))

  % Dirichlet condition p(R, theta) = 0 to O(eps) fixes A
  bc  = subs(p, r, R);
  bc1 = simplify(subs(diff(bc, epsl), epsl, 0) / cos(m*theta));
  Asol = solve(bc1 == 0, A, 'IgnoreAnalyticConstraints', true);
  fprintf('A = %s\n', char(simplify(subs(Asol, lam, lamH))))

  % Kinematic condition, expanded to O(eps)
  Rth  = diff(R, theta);
  kin  = Rt + subs(diff(p, r) - Rth*diff(p, theta)/r^2, r, R);
  kin  = subs(kin, A, Asol);
  kin0 = subs(kin, {epsl, epsdot}, {0, 0});   % epsdot = O(eps)

  % O(1): base growth law
  Hdot_sol = simplify(subs(solve(kin0 == 0, Hdot), lam, lamH));
  fprintf('O(1):   dH/dt = %s\n', char(Hdot_sol))

  % O(eps): epsdot = sigma eps, with epsdot = O(eps) so set epsdot -> sigma*eps
  syms sig real
  kin1 = subs(diff(subs(kin, {epsdot, Hdot}, {sig*epsl, Hdot_sol}), epsl), epsl, 0);
  kin1 = simplify(subs(kin1, lam, lamH) / cos(m*theta));
  sig_sol = simplify(solve(kin1 == 0, sig));
  sigma.(name) = sig_sol;

  fprintf('O(eps): sigma = %s\n', char(sig_sol))
  fprintf('        sigma - expected = %s\n', ...
          char(simplify(sig_sol - cases{ic, 4})))
end

fprintf('\nsigma_rel - sigma_abs = %s\n', char(simplify(sigma.rel - sigma.abs)))

% ── 4. Explicit holomorphic function for the perturbed domain ────────────────
% Claim: for R = H (1 + eps cos(m theta)),
%   f(z) = (1 - 8c)^(1/2) (z/H) [1 - eps beta (z/H)^m] + O(eps^2),
%   beta = (1 - 8c) / (1 + 4c(m - 1)),
% Check by substituting directly into the Liouville formula (lambda
% eliminated) and verifying p = 0 on r = R and the growth rate, to O(eps).

beta = (1 - 8*c) / (1 + 4*c*(m - 1));
g  = @(k) (1 - 8*c)^(1/2) * (zk(1)/H)  * (1 - epsl*beta*(zk(m)/H^m));
gb = @(k) (1 - 8*c)^(1/2) * (zbk(1)/H) * (1 - epsl*beta*(zbk(m)/H^m));
gp  = (1 - 8*c)^(1/2) / H * (1 - epsl*beta*(m+1)*zk(m)/H^m);
gpb = (1 - 8*c)^(1/2) / H * (1 - epsl*beta*(m+1)*zbk(m)/H^m);

pf = subs(liouville(g(1)*gb(1), gp*gpb), lam, lamH);
pf0 = simplify(subs(pf, epsl, 0));
pf1 = subs(diff(pf, epsl), epsl, 0);
pf1 = simplify(rewrite(pf1, 'cos'));
pf1 = simplify(subs(pf1, exp(1i*m*theta), cos(m*theta) + 1i*sin(m*theta)));

fprintf('\n── Explicit f\n')
fprintf('O(1) pressure - p0:            %s\n', char(simplify(pf0 - subs(p0, lam, lamH))))

Rrel = H*(1 + epsl*cos(m*theta));
pR   = subs(pf0 + epsl*pf1, r, Rrel);
fprintf('O(eps) of p on r = R:          %s\n', ...
        char(simplify(subs(diff(pR, epsl), epsl, 0))))

% Growth rate from the explicit pressure
syms sig real
kinf = Hdot*(1 + epsl*cos(m*theta)) + H*sig*epsl*cos(m*theta) ...
       + subs(diff(pf0 + epsl*pf1, r) ...
              - diff(Rrel, theta)*diff(pf0 + epsl*pf1, theta)/r^2, r, Rrel);
kinf = subs(kinf, Hdot, 4*c*H);
kinf1 = simplify(subs(diff(kinf, epsl), epsl, 0) / cos(m*theta));
sig_f = simplify(solve(kinf1 == 0, sig));
fprintf('sigma from explicit f - sigma_rel: %s\n', char(simplify(sig_f - sigma.rel)))

% Explicit p1(r) (coefficient of eps cos(m theta))
p1explicit = 16*c^2*H^2/(1 + 4*c*(m-1)) * (r/H)^m ...
             * ((m+1) - (m-1)*(1-8*c)*(r/H)^2) / (1 - (1-8*c)*(r/H)^2);
fprintf('p1 from f - p1 (explicit):     %s\n', ...
        char(simplify(pf1/cos(m*theta) - p1explicit)))

% ── 5. Limits ────────────────────────────────────────────────────────────────

fprintf('\nLambda = 0 (c = 1/8):\n')
fprintf('  sigma_abs = %s\n', char(simplify(subs(sigma.abs, c, sym(1)/8))))
fprintf('  sigma_rel = %s\n', char(simplify(subs(sigma.rel, c, sym(1)/8))))

fprintf('m = 1:\n')
fprintf('  sigma_abs = %s\n', char(simplify(subs(sigma.abs, m, 1))))
fprintf('  sigma_rel = %s\n', char(simplify(subs(sigma.rel, m, 1))))

fprintf('Lambda >> 1 (c -> 0), leading order:\n')
fprintf('  sigma_abs ~ %s\n', char(taylor(sigma.abs, c, 0, 'Order', 3)))
fprintf('  sigma_rel ~ %s\n', char(taylor(sigma.rel, c, 0, 'Order', 2)))
