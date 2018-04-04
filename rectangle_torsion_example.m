function PD=rectangle_torsion_example(a,b)

G = 70e9;
alpha = 0.1;

PD.VertexList = [-a, -b;
                 a, -b;
                 a, b;
                 -a, b;
                 -a, -b];
PD.InitEdgeLen = 0.125*min(a,b);
PD.BBox = [-a, -b; a, b];
PD.RHS = -2*G*alpha;

PD=PD_torsion_poly(PD,1);
