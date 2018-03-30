function [M, K, Fq] = assemble_PD_frame_fem(PD)
%
% Function to assemble FEM matrices of 1D frame problems
% using a general mesh specified with a PD data structure
%
% Synopsis:
%     [M, K, Fq]  =   assemble_PD_beam_fem(PD)
%
% Input:
%     PD         =   Matlab structure with the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx1 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx2 matrix of node numbers (for each element)
%        NM         =   Number of Materials Sets (MatsSets)
%        MatsSets   =   NMx1 matrix of Matlab structures containing (at least)
%                       the following fields
%           E          =   Young's Modulus
%           G          =   Shear Modulus
%           A          =   X-sec. area
%           J          =   X-sec. torsional constant
%           Ix         =   X-sec. moment of inertia about x axis
%           Iy         =   X-sec. moment of inertia about y axis
%           Ixy        =   X-sec. cross moment of inertia
%           MassMoment =   X-sec. mass moment of inertia per unit length
%           rho        =   Mass density per unit length
%        ElmMats    =   NEx1 matrix of Material set (for each element)
%        qx         =   Function handle to distributed shear force in the x-dir
%                       per unit length qx(z) as @(z)()
%        qy         =   Function handle to distributed shear force in the y-dir
%                       per unit length qy(z) as @(z)()
%        qz         =   Function handle to distributed axial force in the z-dir
%                       per unit length qz(z) as @(z)()
%        EqnNumbering =   Node to global DOF numbering as @(Node, NodeDOF)()
%
% Output: (NTot = N * NodeDOFs)
%     M          =   (NTot)x(NTot) lumped global mass matrix
%     K          =   (NTot)x(NTot) global stiffness matrix
%     Fq         =   (NTot)x1 global force vector
%
%
% By: Ryan S. Elliott -- Jan. 2015
%

NodeDOFs = 6;

NTot = NodeDOFs * PD.N;

M = zeros(NTot, NTot);
K = zeros(NTot, NTot);
Fq = zeros(NTot, 1);

for i = 1:PD.NE
  % Set local element i node numbers
  Node1 = PD.ElmConnect(i, 1);
  Node2 = PD.ElmConnect(i, 2);
  % Set local element i node positions
  Pos1 = PD.NodePos(Node1);
  Pos2 = PD.NodePos(Node2);
  % Set local element i values for MatsSet
  MatSet = PD.ElmMats(i);
  Material = PD.MatsSets(MatSet);
  % compute element stiffness matrix
  k = local_stiffness(Pos1, Pos2, Material);
  % compute element mass matrix
  m = local_mass(Pos1, Pos2, Material);
  % compute elment i local force vector using 5 point Gaussian Quadrature
  ff = quadrature(Pos1, Pos2, PD.qx, PD.qy, PD.qz);

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

  % add element i contribution to global force vector
  Fq(Range) = Fq(Range) + ff;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = local_stiffness(Pos1, Pos2, Material)
%   function to compute local element stiffness matrix
%

% Sanity check on Material
if (Material.Ixy != 0)
  error('Ixy is non-zero.  Must use principal axes.');
end;


NodeDOFs = 6; % [ux, uy, uz, Tx, Ty, Tz] (T is rotation)

k = zeros(2*NodeDOFs, 2*NodeDOFs);

% Bar elemnet
BarRange = [3, 9];
% Torsion element
TorsionRange = [6, 12];
% Beam bending in x-dir
BeamXRange = [1, 4, 7, 10];
% Beam bending in y-dir
BeamYRange = [2, 5, 8, 11];

k(BarRange,BarRange) = k(BarRange,BarRange) ...
 + local_bar_stiffness(Pos1, Pos2, Material.E, Material.A);

k(TorsionRange, TorsionRange) = k(TorsionRange, TorsionRange) ...
 + local_bar_stiffness(Pos1, Pos2, Material.G, Material.J);

k(BeamXRange, BeamXRange) = k(BeamXRange, BeamXRange) ...
 + local_beam_stiffness(Pos1, Pos2, Material.E, Material.Iy);

k(BeamYRange, BeamYRange) = k(BeamYRange, BeamYRange) ...
 + local_beam_stiffness(Pos1, Pos2, Material.E, Material.Ix);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = local_bar_stiffness(Pos1, Pos2, E, A)
%   function to compute local element stiffness matrix for a bar

len = Pos2 - Pos1;

k = (E*A/len) * [[1.0, -1.0];[-1.0, 1.0]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = local_beam_stiffness(Pos1, Pos2, E, I)
%   function to compute local element stiffness matrix for a beam bending in
%   the xz plane

len = Pos2 - Pos1;

k = (E*I/len^3) * ...
  [[     12.0,   6.0*len,    -12.0,   6.0*len];
   [  6.0*len, 4.0*len^2, -6.0*len, 2.0*len^2];
   [    -12.0,  -6.0*len,     12.0,  -6.0*len];
   [  6.0*len, 2.0*len^2, -6.0*len, 4.0*len^2]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = local_mass(Pos1, Pos2, Material)
%   function to compute local element mass matrix

NodeDOFs = 6; % [ux, uy, uz, Tx, Ty, Tz] (T is rotation)

m = zeros(2*NodeDOFs, 2*NodeDOFs);

% Bar elemnet
BarRange = [3, 9];
% Torsion element
TorsionRange = [6, 12];
% Beam bending in x-dir
BeamXRange = [1, 4, 7, 10];
% Beam bending in y-dir
BeamYRange = [2, 5, 8, 11];

m(BarRange,BarRange) = m(BarRange,BarRange) ...
 + local_bar_mass(Pos1, Pos2, Material.rho);

m(TorsionRange, TorsionRange) = m(TorsionRange, TorsionRange) ...
 + local_bar_mass(Pos1, Pos2, Material.MassMoment);

m(BeamXRange, BeamXRange) = m(BeamXRange, BeamXRange) ...
 + local_beam_mass(Pos1, Pos2, Material.rho);

m(BeamYRange, BeamYRange) = m(BeamYRange, BeamYRange) ...
 + local_beam_mass(Pos1, Pos2, Material.rho);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = local_bar_mass(Pos1, Pos2, rho)
%   function to compute local element mass matrix for a bar

len = Pos2 - Pos1;

m = (len*rho/2.0) * [[1.0, 0.0];[0.0, 1.0]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = local_beam_mass(Pos1, Pos2, rho)
%   function to compute local element mass matrix for a beam

len = Pos2 - Pos1;

m = (len*rho/420.0) * ...
  [[    156.0,   22.0*len,      54.0,  -13.0*len],
   [ 22.0*len,  4.0*len^2,  13.0*len, -3.0*len^2],
   [     54.0,   13.0*len,     156.0,  -22.0*len],
   [-13.0*len, -3.0*len^2, -22.0*len,  4.0*len^2]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ff = quadrature(Pos1, Pos2, qx, qy, qz)
%   function to integrate using 5 point GQ the force vector on a frame element

NodeDOFs = 6; % [ux, uy, uz, Tx, Ty, Tz] (T is rotation)

ff = zeros(2*NodeDOFs, 1);

% Bar elemnet
BarRange = [3, 9];
% Torsion element
TorsionRange = [6, 12];
% Beam bending in x-dir
BeamXRange = [1, 4, 7, 10];
% Beam bending in y-dir
BeamYRange = [2, 5, 8, 11];

ff(BarRange) = ff(BarRange) + bar_quadrature(Pos1, Pos2, qz);

ff(BeamYRange) = ff(BeamYRange) + beam_quadrature(Pos1, Pos2, qy);

ff(BeamXRange) = ff(BeamXRange) + beam_quadrature(Pos1, Pos2, qx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ff = bar_quadrature(Pos1, Pos2, f)
%   function to integrate using 5 point GQ the force vector on a bar element

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
function ff = beam_quadrature(Pos1, Pos2, f)
%   function to integrate using 5 point GQ the force vector on beam element

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
