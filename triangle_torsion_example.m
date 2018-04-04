function PD=triangle_torsion_example()

a = 0.01;
G = 70e9;
alpha = 0.1;

PD.VertexList = [-a, -a;
                 a, -a;
                 0, 2*a;
                 -a, -a];
PD.InitEdgeLen = 0.5*a;
PD.BBox = [-a, -a; a, 2*a];
PD.RHS = -2*G*alpha;

PD=PD_torsion_poly(PD,1);
