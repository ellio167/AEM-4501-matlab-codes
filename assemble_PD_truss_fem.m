function [M, K] = assemble_PD_truss_fem(PD)
%
% Function to assemble FEM matrices of 3D truss problems
% using a general mesh specified with a PD data structure
%
% Synopsis:
%     [M, K]  =   assemble_PD_truss_fem(PD)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx3 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx2 matrix of node numbers (for each element)
%        NM         =   Number of Materials Sets (MatsSets)
%        MatsSets   =   NMx1 matrix of Matlab structures containing (at least)
%                       the following fields
%           E          =   Young's Modulus
%           A          =   X-sec. area
%           rho        =   Mass density per unit length
%        ElmMats    =   NEx1 matrix of Material set (for each element)
%        EqnNumbering = Node to global DOF numbering as @(Node, NodeDOF)()
%
% Output: (NTot = N * NodeDOFs)
%     M          =   (NTot)x(NTot) lumped global mass matrix
%     K          =   (NTot)x(NTot) global stiffness matrix
%
%
% By: Ryan S. Elliott -- Jan. 2015
%

NodeDOFs = 3;

NTot = NodeDOFs * PD.N;

M = zeros(NTot, NTot);
K = zeros(NTot, NTot);

for i = 1:PD.NE
  % Set local element i node numbers
  Node1 = PD.ElmConnect(i, 1);
  Node2 = PD.ElmConnect(i, 2);
  % Set local element i node positions
  Pos1 = PD.NodePos(Node1,:)';
  Pos2 = PD.NodePos(Node2,:)';
  % Set local element i values for MatsSet
  MatSet = PD.ElmMats(i);
  Material = PD.MatsSets(MatSet);
  % compute element stiffness and mass matrix
  [m, k] = truss_local_matrices(Pos1, Pos2, Material);

  % Set global connectivity for element
  G = zeros(NodeDOFs,2);
  for i = 1:NodeDOFs
    G(i,1) = PD.EqnNumbering(Node1, i);
    G(i,2) = PD.EqnNumbering(Node2, i);
  end;

  % Define Range variable for element
  Range = [G(:,1); G(:,2)];

  % add element i contribution to global mass matrix
  M(Range, Range) = M(Range, Range) + m;

  % add element i contribution to global stiffness matrix
  K(Range, Range) = K(Range, Range) + k;
end;
