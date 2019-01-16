function PD=rectangle_elasticity_example()

VertList = [0,0;
            1,0;
            1,2;
            0,2;
            0,0];
PD.DistFunc = @(p)(dpoly(p,VertList));
PD.InitEdgeLen = 0.5;
PD.BBox = [-1,-1;3,3];

lambda = 1.0;
mu = 1.0;
PD.Material.Stiffness = ...
[lambda+2*mu, lambda, 0;
 lambda, lambda+2*mu, 0;
 0, 0, mu];
PD.Material.Density = 1.0;
PD.Material.Thickness = 1.0;

PD=PD_elasticity_modes(PD,10,1);
