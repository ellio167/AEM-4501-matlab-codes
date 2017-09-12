function [Disp, Reac, Force, Strain] = bar_static(N, L, E, A, Fz, qz)
%
% Function to solve static FEM 1D tension/torsion problems
% using a uniform mesh.
%
% Synopsis:
%     [Disp, Reac, Force, Strain]  =  bar_static(N, L, E, A, Fz, qz)
%
% Input:
%     N          =   Number of Nodes (number of elements is N-1)
%     L          =   Length of bar
%     E          =   Young's (tension) or Shear Modulus (torsion)
%     A          =   X-section Area (tension) or polar moment (torsion)
%     Fz         =   End load force (tension) or torque (torsion)
%     qz         =   Function handle to distributed force per unit length qz(z)
%                    or torque per unit length Mx(z) as @(z)()
%
% Output:
%     Disp       =   Nx2 matrix with nodal positions and displacements
%     Reac       =   Nx1 matrix with nodal reaction forces
%     Force      =   (N-1)x1 vector with element forces
%     Strain     =   (N-1)x1 vector with element strains
%
%
% By: Ryan S. Elliott -- Jan. 2015

Disp = zeros(N,2);
Force = zeros(N-1,1);
Strain = zeros(N-1,1);

% Compute global values and ignore M
[M, K, GFqz] = assemble_bar_fem(N, L, E, A, 1.0, qz);

% Apply global force boundary conditions
%
GFz = zeros(N,1);
GFz(N) = Fz;

F = GFz + GFqz;

% Apply global displacenet boundary conditions
%
Ftilde = F;
Ktilde = K;

% Set Ftilde(1) to value of node 1 DISPLACEMENT*Ktilde(1,1)
Ftilde(1) = 0.0 * Ktilde(1,1);
% Set Ktilde(1,:) to provide the equation Ktilde(1,1)*U1 = Ftilde(1)*Ktilde(1,1)
Ktilde(1,2:N) = 0.0;

% set nodal locations
for i = 1:N
  Disp(i,1) = (L*(i-1))/(N-1);
end;

% find displacements
Disp(:,2) = linsolve(Ktilde,Ftilde);

% find reactions
Reac = -F + K*Disp(:,2);

% compute strains and forces
len = L/(N-1);
for i = 1:(N-1)
  Strain(i) = (Disp(i+1,2) - Disp(i,2))/len;
  Force(i) = Strain(i)*E*A;
end;
