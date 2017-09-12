function [Freq, Modes] = gmesh_bar_modes(N, NodePos, ...
                                         NE, ElmConnect, ...
                                         E, A, rho, ...
                                         LEndNode, REndNode, Nmodes)
%
% Function to FEM solve for free vib. modes of 1D tension/torsion problems
% using a general mesh.
%
% Synopsis:
%     [Freq, Modes] = gmesh_bar_modes(N, NodePos, ...
%                                     NE, ElmConnect, ...
%                                     E, A, rho, ...
%                                     LEndNode, REndNode, Nmodes)
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
%     LEndNode   =   Node number for Left hand end of bar
%     REndNode   =   Node number for Right hand end of bar
%     Nmodes     =   Number of modes to solve for
%
% Output:
%     Freq       =   (Nmodes)x1 vector of vibration frequencies
%     Modes      =   (N)x(Nmodes) matrix of vibration modes
%
%
% By: Ryan S. Elliott -- Jan. 2015

Freq = zeros(Nmodes,1);
Modes = zeros(N,Nmodes);

% Compute global values (and ignore F; which will be zero)
[M, K, GFqz] = assemble_gmesh_bar_fem(N, NodePos, ...
                                      NE, ElmConnect, ...
                                      E, A, rho, ...
                                      @(z)(0));

% Apply glabal displacement boundary conditions (left-fixed; right-free)
% work with deflated matrices to eliminate homogeneous BC at left end
%
Range = [1:LEndNode-1, LEndNode+1:N];
Mhat = M(Range, Range);
Khat = K(Range, Range);

[V,D] = eigs(Khat, Mhat, Nmodes,'sm');

Freq = sqrt(diag(D));
Modes(Range,:) = V;
