
addpath('../scripts')

sphereConvergence
print(figSpatial, 'output/fig-sphere-spatial-convergence', '-dpdf')
print(figTemporal, 'output/fig-sphere-temporal-convergence', '-dpdf')

roughnessConvergence
print(fig, 'output/fig-roughness-convergence', '-dpdf')

growthConvergence
print(fig, 'output/fig-anisotropic-temporal-convergence', '-dpdf')

homogenizationConvergence
print(fig, 'output/fig-homogenization-eps-convergence', '-dpdf')
