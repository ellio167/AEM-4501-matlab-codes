function PD=butterfly_elasticity_example()
% Example for how to solve for free vibration with a butterfly shape


PD.DistFunc = @butterfly_func;

PD.InitEdgeLen = 0.5;
PD.BBox = [-10,-10;10,10];

lambda = 1.0;
mu = 1.0;

%plain strain for isotropic material
PD.Material.Stiffness = ...
[lambda+2*mu,      lambda,   0;
      lambda, lambda+2*mu,   0;
           0,           0,  mu];
PD.Material.Density = 1.0;
PD.Material.Thickness = 1.0;

PD=PD_elasticity_modes(PD,1);

function fd = butterfly_func(p)
[TH,R] = cart2pol(p(:,1), p(:,2));

rth = 0.5 * (12.0 - 0.5 * sin(TH) + 2.5 * sin(3.0 * TH) + ...
      2.0 * sin(5.0 * TH) - 1.7 * sin(7.0 * TH) + 3.0 * cos(2.0 * TH) - ...
      2.0 * cos(4.0 * TH) - 0.4 * cos(16 * TH));

fd = R - rth;
return
