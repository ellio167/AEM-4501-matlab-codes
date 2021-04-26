function PD=rectangle_elasticity_example(lambda)
% Example for how to solve for free vibration with a rectangle shape

 VertList = [-.5,-1;
             .5,-1;
             .5,1;
             -.5,1;
             -.5,-1];

PD.DistFunc = @(p)(dpoly(p,VertList));

PD.InitEdgeLen = 0.125;
PD.BBox = [-2,-2;2,2];

% use lambda from user argument instead of below value
% lambda = 1.0;
mu = 1.0;

%plain strain for isotropic material
PD.Material.Stiffness = ...
[lambda+2*mu,      lambda,   0;
      lambda, lambda+2*mu,   0;
           0,           0,  mu];
PD.Material.Density = 1.0;
PD.Material.Thickness = 1.0;

PD=PD_elasticity_modes(PD,1);
