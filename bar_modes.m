function [Freq, Modes] = bar_modes(N, L, E, A, rho, Nmodes)
%
% Function to FEM solve for free vib. modes of 1D tension/torsion problems
% using a uniform mesh.
%
% Synopsis:
%     [Freq, Modes]  =  bar_modes(N, L, E, A, rho, Nmodes)
%
% Input:
%     N          =   Number of Nodes (number of elements is N-1)
%     L          =   Length of bar
%     E          =   Young's (tension) or Shear Modulus (torsion)
%     A          =   X-section Area (tension) or polar moment (torsion)
%     rho        =   Mass density per unit length (tension)
%                    or centrodial moment of inertia (torsion) per unit length
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
[M, K, GFqz] = assemble_bar_fem(N, L, E, A, rho, @(z)(0.0));

% Apply glabal displacement boundary conditions
% work with deflated matrices to eliminate homogeneous BC at left end
%
Mhat = M(2:N,2:N);
Khat = K(2:N,2:N);

[V,D] = eigs(Khat, Mhat, Nmodes,'sm');

Freq = sqrt(diag(D));
Modes = [zeros(1,Nmodes); V];
