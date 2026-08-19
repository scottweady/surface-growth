
% MAKEDATA  Runs each solve script and saves its output for the figures.
%
% Each solve script is run only if its output isn't already on disk; delete
% the corresponding figN.mat to force a re-solve.

addpath('../scripts')

if ~isfile('fig2.mat')
  sphereProfiles
  save('fig2.mat', 'hAll', 'pAll', 'uAll')
end

if ~isfile('fig3.mat')
  sphereGrowth
  save('fig3.mat', 'Hplt', 'ts', 'AAall', 'AAinf')
end

if ~isfile('fig4.mat')
  sphereGrowthRate
  save('fig4.mat', 'AsaveAll', 'Amax')
end

if ~isfile('fig5.mat')
  roughnessGrowthRate
  save('fig5.mat', 'hom', 'X', 'Y', 'U', 'V', 'G')
end

if ~isfile('fig6/isotropic.mat') || ~isfile('fig6/anisotropic.mat')
  if ~exist('fig6', 'dir')
    mkdir('fig6')
  end
  movingBoundaryRoughness
  save('fig6/isotropic.mat', '-struct', 'iso')
  save('fig6/anisotropic.mat', '-struct', 'ani')
end
