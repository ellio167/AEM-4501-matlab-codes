function PD=ellipse_torsion_example()

a=4;
b=1;
G=70e9;
alpha=0.1;

PD.DistFunc = @(p)(p(:,1)/a).^2 + (p(:,2)/b).^2 - 1;
PD.InitEdgeLen = 0.25*b;
PD.BBox = [-a,-b; a,b];
PD.RHS = -2*G*alpha;

PD=PD_torsion(PD,1);
