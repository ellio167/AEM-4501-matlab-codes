function PD=butterfly_torsion_example()

PD.DistFunc = @butterfly_func;
PD.InitEdgeLen = 0.25;
PD.BBox = [-15,-15;15,15];
PD.RHS = -2*70e9*(5*pi/180);

PD=PD_torsion(PD,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fd = butterfly_func(p)
[TH,R] = cart2pol(p(:,1), p(:,2));

rth = 0.5 * (12.0 - 0.5 * sin(TH) + 2.5 * sin(3.0 * TH) + ...
      2.0 * sin(5.0 * TH) - 1.7 * sin(7.0 * TH) + 3.0 * cos(2.0 * TH) - ...
      2.0 * cos(4.0 * TH) - 0.4 * cos(16 * TH));
      
fd = R - rth;  

return