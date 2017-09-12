function [M, K, Fqz] = assemble_gmesh_bar_fem(N, NodePos, ...
                                              NE, ElmConnect, ...
                                              E, A, rho, ...
                                              qz)
%
% Function to assemble FEM matrices of 1D tension/torsion problem
% using a general mesh.
%
% Synopsis:
%     [M, K, Fqz]  =   assemble_gmesh_bar_fem(N, NodePos, ...
%                                             NE, ElmConnect, ...
%                                             E, A, rho, ...
%                                             qz)
%
% Input:
%     N          =   Number of nodes in mesh (numbered 1:N)
%     NodePos    =   Nx1 matrix of nodal positions
%     NE         =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect =   NEx2 matrix of node numbers (for each element)
%     E          =   NEx1 matrix Young's (tension) or Shear Modulus (torsion)
%     A          =   NEx1 matrix X-sec. area (tension) or polar moment (torsion)
%     rho        =   NEx1 matrix Mass density per unit length (tension)
%                    or centrodial moment of inertia (torsion) per unit length
%     qz         =   Function handle to distributed force per unit length
%                    or torque per unit length qz(z) as @(z)()
%
% Output:
%     M          =   NxN lumped global mass matrix
%     K          =   NxN global stiffness matrix
%     Fqz        =   Nx1 global force vector
%
%
% By: Ryan S. Elliott -- Jan. 2015
%

M = zeros(N,N);
K = zeros(N,N);
Fqz = zeros(N,1);

for i = 1:NE
  % Set local element i node numbers
  Node1 = ElmConnect(i, 1);
  Node2 = ElmConnect(i, 2);
  % Set local element i node positions
  Pos1 = NodePos(Node1);
  Pos2 = NodePos(Node2);
  % Set local element i values for E, A, rho
  e = E(i);
  a = A(i);
  r = rho(i);
  % compute element stiffness matrix
  k = local_stiffness(Pos1, Pos2, e, a);
  % compute element mass matrix
  m = local_mass(Pos1, Pos2, r);
  % compute elment i local force vector using 5 point Gaussian Quadrature
  ff = quadrature(Pos1, Pos2, qz);

  % Set global connectivity for element
  G1 = Node1;
  G2 = Node2;

  % Define Range variable for element
  Range = [G1, G2];

  % add element i contribution to global mass matrix
  M(Range, Range) = M(Range, Range) + m;

  % add element i contribution to global stiffness matrix
  K(Range, Range) = K(Range, Range) + k;

  % add element i contribution to global force vector
  Fqz(Range) = Fqz(Range) + ff;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = local_stiffness(Pos1, Pos2, E, A)
%   function to compute local element stiffness matrix

len = Pos2 - Pos1;

k = (E*A/len) * [[1.0, -1.0];[-1.0, 1.0]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = local_mass(Pos1, Pos2, rho)
%   function to compute local element mass matrix

len = Pos2 - Pos1;

m = (len*rho/2.0) * [[1.0, 0.0];[0.0, 1.0]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ff = quadrature(Pos1, Pos2, f)
%   function to integrate using 5 point GQ the force vector on element

% shape functions
N1 = @(s)((1.0-s)/2.0);
N2 = @(s)((s+1.0)/2.0);

% local nodal locations
lz(1) = Pos1;
lz(2) = Pos2;

% Mapping from master element coordinate (s) to global coordinate (z)
z = @(s)((lz(1)+lz(2))/2.0 + s*(lz(2)-lz(1))/2.0);
% Jacobian of mapping
jac = (lz(2)-lz(1))/2.0;

% GQ points
s(1) = -1.0;
s(2) = -sqrt(3.0/7.0);
s(3) = 0.0;
s(4) = sqrt(3.0/7.0);
s(5) = 1.0;
% GQ weights
w(1) = 1.0/10.0;
w(2) = 49.0/90.0;
w(3) = 32.0/45.0;
w(4) = 49.0/90.0;
w(5) = 1.0/10.0;
% compute function values at quadrature points
gq(1) = z(s(1));
gq(2) = z(s(2));
gq(3) = z(s(3));
gq(4) = z(s(4));
gq(5) = z(s(5));

ff = zeros(2,1);
for i = 1:5
  ff(1) = ff(1) + jac*w(i)*f(gq(i))*N1(s(i));
  ff(2) = ff(2) + jac*w(i)*f(gq(i))*N2(s(i));
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
