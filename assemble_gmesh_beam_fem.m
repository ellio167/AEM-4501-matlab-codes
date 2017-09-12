function [M, K, Fqy] = assemble_gmesh_beam_fem(N, NodePos, ...
                                               NE, ElmConnect, ...
                                               E, Ix, rho, ...
                                               qy, EqnNumbering)
%
% Function to assemble FEM matrices of plane beam problems
% using a general mesh.
%
% Synopsis:
%     [M, K, Fqy]  =   assemble_gmesh_beam_fem(N, NodePos, ...
%                                              NE, ElmConnect, ...
%                                              E, Ix, rho, ...
%                                              qy, EqnNumbering)
%
% Input:
%     N            =   Number of nodes in mesh (numbered 1:N)
%     NodePos      =   Nx1 matrix of nodal positions
%     NE           =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect   =   NEx2 matrix of node numbers (for each element)
%     E            =   NEx1 matrix Young's Modulus
%     Ix           =   NEx1 matrix X-sec. moment of inertia
%     rho          =   NEx1 matrix Mass density per unit length
%     qy           =   Function handle to distributed shear force
%                      per unit length qy(z) as @(z)()
%     EqnNumbering =   Node to global DOF numbering as @(Node, NodeDOF)()
%
% Output: (NTot = N * NodeDOFs)
%     M          =   (NTot)x(NTot) lumped global mass matrix
%     K          =   (NTot)x(NTot) global stiffness matrix
%     Fqy        =   (NTot)x1 global force vector
%
%
% By: Ryan S. Elliott -- Jan. 2015
%

NodeDOFs = 2;

NTot = NodeDOFs * N;

M = zeros(NTot, NTot);
K = zeros(NTot, NTot);
Fqy = zeros(NTot, 1);

for i = 1:NE
  % Set local element i node numbers
  Node1 = ElmConnect(i, 1);
  Node2 = ElmConnect(i, 2);
  % Set local element i node positions
  Pos1 = NodePos(Node1);
  Pos2 = NodePos(Node2);
  % Set local element i values for E, Ix, rho
  e = E(i);
  inertia = Ix(i);
  r = rho(i);
  % compute element stiffness matrix
  k = local_stiffness(Pos1, Pos2, e, inertia);
  % compute element mass matrix
  m = local_mass(Pos1, Pos2, r);
  % compute elment i local force vector using 5 point Gaussian Quadrature
  ff = quadrature(Pos1, Pos2, qy);

  % Set global connectivity for element
  G1 = EqnNumbering(Node1, 1);
  G2 = EqnNumbering(Node1, 2);
  G3 = EqnNumbering(Node2, 1);
  G4 = EqnNumbering(Node2, 2);

  % Define Range variable for element
  Range = [G1, G2, G3, G4];

  % add element i contribution to global mass matrix
  M(Range, Range) = M(Range, Range) + m;

  % add element i contribution to global stiffness matrix
  K(Range, Range) = K(Range, Range) + k;

  % add element i contribution to global force vector
  Fqy(Range) = Fqy(Range) + ff;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = local_stiffness(Pos1, Pos2, E, I)
%   function to compute local element stiffness matrix

len = Pos2 - Pos1;

k = (E*I/len^3) * ...
  [[     12.0,   6.0*len,    -12.0,   6.0*len];
   [  6.0*len, 4.0*len^2, -6.0*len, 2.0*len^2];
   [    -12.0,  -6.0*len,     12.0,  -6.0*len];
   [  6.0*len, 2.0*len^2, -6.0*len, 4.0*len^2]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = local_mass(Pos1, Pos2, rho)
%   function to compute local element mass matrix

len = Pos2 - Pos1;

m = (len*rho/420.0) * ...
  [[    156.0,   22.0*len,      54.0,  -13.0*len],
   [ 22.0*len,  4.0*len^2,  13.0*len, -3.0*len^2],
   [     54.0,   13.0*len,     156.0,  -22.0*len],
   [-13.0*len, -3.0*len^2, -22.0*len,  4.0*len^2]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ff = quadrature(Pos1, Pos2, f)
%   function to integrate using 5 point GQ the force vector on element

len = Pos2 - Pos1;

% shape functions
N1 = @(s)(0.25*(1.0-s)^2*(2.0+s));
N2 = @(s)(0.125*len*(1.0-s)^2*(1.0+s));
N3 = @(s)(0.25*(1.0+s)^2*(2.0-s));
N4 = @(s)(0.125*len*(1.0+s)^2*(s-1.0));

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

ff = zeros(4,1);
for i = 1:5
  ff(1) = ff(1) + jac*w(i)*f(gq(i))*N1(s(i));
  ff(2) = ff(2) + jac*w(i)*f(gq(i))*N2(s(i));
  ff(3) = ff(3) + jac*w(i)*f(gq(i))*N3(s(i));
  ff(4) = ff(4) + jac*w(i)*f(gq(i))*N4(s(i));
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
