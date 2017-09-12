function PD = create_PD_for_uniform_mesh(Nodes, L, E, G, A, J, Ix, Iy, Ixy, ...
                                         rho, MassMoment)
%
% Function to create PD mesh data structures for a uniform mesh
%
% Synopsis:
%     PD = create_PD_for_uniform_mesh(Nodes, L, E, G, A, J, Ix, Iy, Iz, ...
%                                     rho, MassMoment)
%
% Input:
%     Nodes      =   Number of Nodes (number of elements is N-1)
%     L          =   Length of bar
%     E          =   Young's Modulus
%     G          =   Shear Modulus
%     A          =   X-sec. area
%     J          =   X-sec. torsional constant
%     Ix         =   X-sec. moment of inertia about x axis
%     Iy         =   X-sec. moment of inertia about y axis
%     Ixy        =   X-sec. cross moment of inertia
%     rho        =   Mass density per unit length
%     MassMoment =   X-sec. mass moment of inertia per unit length
%

PD.N = Nodes;
PD.NE = PD.N - 1;

PD.LEndNode = 1;
PD.REndNode = Nodes;

PD.NodePos(1:PD.N,1) = [0:L/PD.NE:L]';

j = 1;
for i = 1:PD.NE
  PD.ElmConnect(i,:) = [j, j+1];
  j = j+1;
end;

Mat.E = E;
Mat.G = G;
Mat.A = A;
Mat.J = J;
Mat.Ix = Ix;
Mat.Iy = Iy;
Mat.Ixy = Ixy;
Mat.rho = rho;
Mat.MassMoment = MassMoment;

PD.NM = 1;
PD.MatsSets(1,1) = Mat;

PD.ElmMats(1:PD.NE,1) = 1;

PD.EqnNumbering = @(Node, NodeDOF)(6*(Node-1) + NodeDOF);