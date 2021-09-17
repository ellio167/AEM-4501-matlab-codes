lambdaInput = 1;
% function PD=sphere_elasticity_example(lambdaInput)
% Example for how to solve for free vibration with a 3d sphere domain

% sphere of radius 1 centered at origin
% PD.DistFunc = @(r)(r(:,1).^2/(1^2) + r(:,2).^2/(1^2) + r(:,3).^2/(1^2) - 1.0);
% PD.DistFunc = inline('sqrt(sum(p.^2,2))-1','p');
PD.DistFunc = @(r)(sqrt(sum(r.^2,2))-1);

PD.InitEdgeLen = 0.2;
PD.BBox = [-1,-1,-1; 1,1,1];

mu = 1;
lambda = lambdaInput*mu;

% 3D for isotropic material
PD.Material.Stiffness = ...
[lambda+2*mu,      lambda,      lambda,   0,   0,   0;
      lambda, lambda+2*mu,      lambda,   0,   0,   0;
      lambda,      lambda, lambda+2*mu,   0,   0,   0;
           0,           0,           0,   mu,  0,   0;
           0,           0,           0,    0, mu,   0;
           0,           0,           0,    0,  0,  mu];

PD.Material.Density = 1.0;

PD=PD_3d_elasticity_modes(PD,1);
