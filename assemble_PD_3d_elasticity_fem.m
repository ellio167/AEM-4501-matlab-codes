function [M, K] = assemble_PD_3d_elasticity_fem(PD)
%
% Function to assemble FEM matrices of 3D elasticity free
% vibration problems using a general mesh specified with a
% PD data structure
%
% Synopsis:
%     [M, K]  =   assemble_PD_3d_elasticity_fem(PD)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx3 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx4 matrix of node numbers (for each element)
%        Material   =   Matlab structure containing (at least)
%                       the following fields
%           Stiffness  =   6x6 symmetric 3D voigt stiffness matrix
%           Density    =   Mass density per unit volume
%        EqnNumbering =   Node to global DOF numbering as @(Node, NodeDOF)()
%
% Output: (NTot = N * NodeDOFs)
%     M          =   (NTot)x(NTot) consistent global mass matrix
%     K          =   (NTot)x(NTot) global stiffness matrix
%
%
% By: Ryan S. Elliott -- Apr 2021
%

NodeDOFs = 3;

NTot = NodeDOFs * PD.N;

M = zeros(NTot, NTot);
K = zeros(NTot, NTot);

for i = 1:PD.NE
  % Set local element i node numbers
  Node1 = PD.ElmConnect(i, 1);
  Node2 = PD.ElmConnect(i, 2);
  Node3 = PD.ElmConnect(i, 3);
  Node4 = PD.ElmConnect(i, 4);
  % Set local element i node positions
  Pos1 = PD.NodePos(Node1,:);
  Pos2 = PD.NodePos(Node2,:);
  Pos3 = PD.NodePos(Node3,:);
  Pos4 = PD.NodePos(Node4,:);
  % compute element stiffness matrix
  [m, k] = local_matrices(Pos1, Pos2, Pos3, Pos4, PD.Material);

  % Set global connectivity for element
  G = zeros(NodeDOFs,4);
  for i = 1:NodeDOFs
    G(i,1) = PD.EqnNumbering(Node1, i);
    G(i,2) = PD.EqnNumbering(Node2, i);
    G(i,3) = PD.EqnNumbering(Node3, i);
    G(i,4) = PD.EqnNumbering(Node4, i);
  end;

  % Define Range variable for element
  Range = [G(:,1); G(:,2); G(:,3); G(:,4)];

  % add element i contribution to global mass matrix
  M(Range, Range) = M(Range, Range) + m;

  % add element i contribution to global stiffness matrix
  K(Range, Range) = K(Range, Range) + k;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m, k] = local_matrices(Pos1, Pos2, Pos3, Pos4, Material)
%   function to compute local element stiffness matrix
%

% NodeDOFs = 3; % [ux, uy, uz]

Xi = [ones(4,1),[Pos1; Pos2; Pos3; Pos4]]; % 4 by 4 with rows=[1 xcorner ycorner zcorner]
Volume = abs(det(Xi))/2.0; % volume of tetrahedral = half of parallepiped volume

% columns of XiInv are coeffs in a+bx+cy+dz to give s=1,0,0,0 at nodes
XiInv = inv(Xi);

% Compute B matrix from slopes b,c,d in gradN = [dN/dx; dN/dy; dN/dz]
%
% voigt strain: [e_xx; e_yy; e_zz; gamma_yz; gamma_zx; gamma_xy]
%
gradN = XiInv(2:4,:);
B(1,:) = [gradN(1,1), 0         , 0         , gradN(1,2), 0         , 0         , gradN(1,3), 0         , 0         , gradN(1,4), 0         , 0         ];
B(2,:) = [0         , gradN(2,1), 0         , 0         , gradN(2,2), 0         , 0         , gradN(2,3), 0         , 0         , gradN(2,4), 0         ];
B(3,:) = [0         , 0         , gradN(3,1), 0         , 0         , gradN(3,2), 0         , 0         , gradN(3,3), 0         , 0         , gradN(3,4)];
B(4,:) = [0         , gradN(3,1), gradN(2,1), 0         , gradN(3,2), gradN(2,2), 0         , gradN(3,3), gradN(2,3), 0         , gradN(3,4), gradN(2,4)];
B(5,:) = [gradN(3,1), 0         , gradN(1,1), gradN(3,2), 0         , gradN(1,2), gradN(3,3), 0         , gradN(1,3), gradN(3,4), 0         , gradN(1,4)];
B(6,:) = [gradN(2,1), gradN(1,1), 0         , gradN(2,2), gradN(1,2), 0         , gradN(2,3), gradN(1,3), 0         , gradN(2,4), gradN(1,4), 0         ];

BDB = B'*(Material.Stiffness)*B;

k = Volume*BDB;

k = 0.5*(k+k');

% See Chandrupatla & Belegundu, 2nd Ed, 1997; pg380 for m formula
m = ((Material.Density)*Volume/20.0) * ...
    [[2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0];
     [0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0];
     [0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1];
     [1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0];
     [0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0];
     [0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1];
     [1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0];
     [0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0];
     [0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1];
     [1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0];
     [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0];
     [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2]];

m = 0.5*(m+m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
