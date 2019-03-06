function [Disp, Reac, Force, Strain] = gmesh_bar_static(N, NodePos, ...
                                                        NE, ElmConnect, ...
                                                        E, A, ...
                                                        LEndNode, REndNode, ...
                                                        Fz, qz)
%
% Function to solve static FEM 1D tension/torsion problems
% using a general mesh.
%
% Synopsis:
%     [Disp, Reac, Force, Strain]  =  gmesh_bar_static(N, NodePos, ...
%                                                      NE, ElmConnect, ...
%                                                      E, A, ...
%                                                      LEndNode, REndNode, ...
%                                                      Fz, qz)
%
% Input:
%     N          =   Number of nodes in mesh (numbered 1:N)
%     NodePos    =   Nx1 matrix of nodal positions
%     NE         =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect =   NEx2 matrix of node numbers (for each element)
%     E          =   NEx1 matrix Young's (tension) or Shear Modulus (torsion)
%     A          =   NEx1 matrix X-sec. area (tension) or polar moment (torsion)
%     LEndNode   =   Node number for Left hand end of bar
%     REndNode   =   Node number for Right hand end of bar
%     Fz         =   End load force (tension) or torque (torsion)
%     qz         =   Function handle to distributed force per unit length
%                    or torque per unit length qz(z) as @(z)()
%
% Output:
%     Disp       =   Nx2 matrix with nodal positions and displacements
%     Reac       =   Nx1 matrix with nodal reaction forces
%     Force      =   NEx1 vector with element forces
%     Strain     =   NEx1 vector with element strain
%
%
% By: Ryan S. Elliott -- Jan. 2015

Disp = zeros(N,2);
Force = zeros(NE,1);
Strain = zeros(NE,1);

% Compute global values and ignore M
[M, K, GFqz] = assemble_gmesh_bar_fem(N, NodePos, ...
                                      NE, ElmConnect, ...
                                      E, A, ones(N,1), ...
                                      qz);

% Apply global force boundary conditions
%
% Add P to F(REndNode)
GFz = zeros(N,1);
GFz(REndNode) = Fz;

F = GFz + GFqz;

% Apply global displacenet boundary conditions
%
Ftilde = F;
Ktilde = K;

% Set Ftilde(LEndNode) to value of node 1 DISPLACEMENT*Ktilde(LEndNode,LEndNode)
Ftilde(LEndNode) = 0.0 * Ktilde(LEndNode,LEndNode);
% Set Ktilde(LEndNode,:) to provide the equation:
%    Ktilde(LEndNode,LEndNode)*U = Ftilde(LEndNode)*Ktilde(LEndNode,LEndNode)
tmp = Ktilde(LEndNode, LEndNode);
Ktilde(LEndNode,1:N) = 0.0;
Ktilde(LEndNode, LEndNode) = tmp;

% set nodal locations
Disp(:,1) = NodePos;

% find displacements
Disp(:,2) = linsolve(Ktilde,Ftilde);

% find reactions
Reac = -F + K*Disp(:,2);

% compute strains and stresses
for i = 1:NE
  Node1 = ElmConnect(i,1);
  Node2 = ElmConnect(i,2);
  Pos1 = NodePos(Node1);
  Pos2 = NodePos(Node2);
  len = Pos2 - Pos1;
  Strain(i) = (Disp(Node2,2) - Disp(Node1,2))/len;
  Force(i) = Strain(i)*E(i)*A(i);
end;
