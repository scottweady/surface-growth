
addpath('../utils')
addpath('../src')
addpath('../scripts')

fprintf('Running sphere convergence\n')
sphere_convergence

fprintf('Running homogenization convergence\n')
homogenization_convergence

fprintf('Running anisotropic convergence\n')
anisotropic_convergence

close all