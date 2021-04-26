function PD=circle_elasticity_example(lambdaInput)
% Example for how to solve for free vibration with a circle

% circle of radius 1 centered at origin
PD.DistFunc = @(r)(r(:,1).^2/(1^2) + r(:,2).^2/(1^2) - 1.0);

PD.InitEdgeLen = 0.15;
PD.BBox = [-2,-2;2,2];

mu = 1;
lambda = lambdaInput*mu;

%plain strain for isotropic material
PD.Material.Stiffness = ...
[lambda+2*mu,      lambda,   0;
      lambda, lambda+2*mu,   0;
           0,           0,  mu];
PD.Material.Density = 1.0;
PD.Material.Thickness = 1.0;

PD=PD_elasticity_modes(PD,0);
