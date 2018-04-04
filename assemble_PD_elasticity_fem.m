function [M, K] = assemble_PD_elasticity_fem(PD)
%
% Function to assemble FEM matrices of 2D elasticity free
% vibration problems using a general mesh specified with a
% PD data structure
%
% Synopsis:
%     [M, K]  =   assemble_PD_elasticity_fem(PD)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx2 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx3 matrix of node numbers (for each element)
%        Material   =   Matlab structure containing (at least)
%                       the following fields
%           Stiffness  =   3x3 symmetric 2D stiffness matrix
%           Density    =   Mass density per unit volume
%           Thickness  =   Thickness in third-direction (1.0 for plane-strain)
%        EqnNumbering =   Node to global DOF numbering as @(Node, NodeDOF)()
%
% Output: (NTot = N * NodeDOFs)
%     M          =   (NTot)x(NTot) consistent global mass matrix
%     K          =   (NTot)x(NTot) global stiffness matrix
%
%
% By: Ryan S. Elliott -- Jan. 2015, Apr. 2018
%

NodeDOFs = 2;

NTot = NodeDOFs * PD.N;

M = zeros(NTot, NTot);
K = zeros(NTot, NTot);

for i = 1:PD.NE
  % Set local element i node numbers
  Node1 = PD.ElmConnect(i, 1);
  Node2 = PD.ElmConnect(i, 2);
  Node3 = PD.ElmConnect(i, 3);
  % Set local element i node positions
  Pos1 = PD.NodePos(Node1,:);
  Pos2 = PD.NodePos(Node2,:);
  Pos3 = PD.NodePos(Node3,:);
  % compute element stiffness matrix
  [m, k] = local_matrices(Pos1, Pos2, Pos3, PD.Material);

  % Set global connectivity for element
  G = zeros(NodeDOFs,3);
  for i = 1:NodeDOFs
    G(i,1) = PD.EqnNumbering(Node1, i);
    G(i,2) = PD.EqnNumbering(Node2, i);
    G(i,3) = PD.EqnNumbering(Node3, i);
  end;

  % Define Range variable for element
  Range = [G(:,1); G(:,2); G(:,3)];

  % add element i contribution to global mass matrix
  M(Range, Range) = M(Range, Range) + m;

  % add element i contribution to global stiffness matrix
  K(Range, Range) = K(Range, Range) + k;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m, k] = local_matrices(Pos1, Pos2, Pos3, Material)
%   function to compute local element stiffness matrix
%

% NodeDOFs = 2; % [ux, uy]

Xi = [ones(3,1),[Pos1; Pos2; Pos3]]; % 3 by 3 with rows=[1 xcorner ycorner]
Area = abs(det(Xi))/2.0; % area of triangle e = half of parallelogram area

% columns of XiInv are coeffs in a+bx+cy to give s=1,0,0 at nodes
XiInv = inv(Xi);

% Compute B matrix from slopes b,c in gradN = [dN/dx; dN/dy]
%
gradN = XiInv(2:3,:);
B = [gradN(1,1), 0,          gradN(1,2), 0,          gradN(1,3), 0;
     0,          gradN(2,1), 0,          gradN(2,2), 0,          gradN(2,3);
     gradN(2,1), gradN(1,1), gradN(2,2), gradN(1,2), gradN(2,3), gradN(1,3)];

BDB = B'*(Material.Stiffness)*B;

k = (Material.Thickness)*Area*BDB;

k = 0.5*(k+k');

m = ((Material.Density)*(Material.Thickness)*Area/12.0) * ...
    [[2, 0, 1, 0, 1, 0];
     [0, 2, 0, 1, 0, 1];
     [1, 0, 2, 0, 1, 0];
     [0, 1, 0, 2, 0, 1];
     [1, 0, 1, 0, 2, 0];
     [0, 1, 0, 1, 0, 2]];

m = 0.5*(m+m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
