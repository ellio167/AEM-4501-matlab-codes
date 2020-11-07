function PD=churro_elasticity_example()
% Example for how to solve for free vibration with a churro shape


PD.DistFunc = @churro_func;

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

function fd = churro_func(p)
[TH,R] = cart2pol(p(:,1), p(:,2));

rth = 7-cos(TH).*sin(7.*TH);

fd = R - rth;
