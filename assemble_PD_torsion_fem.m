function [K, F] = assemble_PD_torsion_fem(PD)
%
% Function to assemble FEM matrices of Prandtl torsion problems
% using a distmesh2d mesh specified with a PD data structure
%
% Synopsis:
%     [K, F]  =   assemble_PD_torsion_fem(PD)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        DistFunc   =   Distance Function handle for 2D x-section domain
%                       as @(x,y)().  The function must be zero on the
%                       boundary, negative inside the domain and positive
%                       outside the domain.
%        InitEdgeLen=   Initial mesh edge length
%        BBox       =   [xmin,ymin; xmax,ymax] Bounding box for domain
%        RHS        =   Right-hand side of torsion equation (-2*G*dTheta/dz)
%        N          =   The number of nodes in the generated mesh
%        NE         =   The number of elements in the generated mesh
%        NodePos    =   Nx2 array of the node positions
%        ElmConnect =   NEx3 array of the element connectivity
%
% Output:
%     K          =   NxN global "stiffness" matrix
%     F          =   Nx1 global "force" vector
%
% By: Ryan S. Elliott, Amartya Banerjee -- Apr. 2015
%

%
% The program PD_torsion.m is partially based on the FEM solver described by
% John Coady https://www.particleincell.com/2012/matlab-fem/
%


K = zeros(PD.N, PD.N);
F = zeros(PD.N, 1);

for i = 1:PD.NE
  % Set local element i node numbers
  Node1 = PD.ElmConnect(i, 1);
  Node2 = PD.ElmConnect(i, 2);
  Node3 = PD.ElmConnect(i, 3);
  % Set local element i node positions
  Pos1 = PD.NodePos(Node1,:);
  Pos2 = PD.NodePos(Node2,:);
  Pos3 = PD.NodePos(Node3,:);
  % compute element "stiffness" matrix and "force" vector
  [k, f] = local_stiffness(Pos1, Pos2, Pos3, PD.RHS);

  % Set global connectivity for element
  G1 = Node1;
  G2 = Node2;
  G3 = Node3;

  % Define Range variable for element
  Range = [G1, G2, G3];

  % add element i contribution to global "stiffness" matrix
  K(Range, Range) = K(Range, Range) + k;

  % add element i contribution to global "force" vector
  F(Range) = F(Range) + f;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k, f] = local_stiffness(Pos1, Pos2, Pos3, rhs)
%   function to compute local element "stiffness" matrix

Pe = [ones(3,1),[Pos1; Pos2; Pos3]]; % 3 by 3 with rows=[1 xcorner ycorner]
Area = abs(det(Pe))/2.0; % area of triangle e = half of parallelogram area

C = inv(Pe); % columns of C are coeffs in a+bx+cy to give s=1,0,0 at nodes

% Now compute 3 by 3 k
% element matrix from slopes b,c in grad
%
grad = C(2:3,:);
k = Area*(grad'*grad);

% Now compute 3 by 1 f
% element vector
%
f = - (Area/3.0) * rhs * ones(3,1); % minus sign comes from div-theorem

