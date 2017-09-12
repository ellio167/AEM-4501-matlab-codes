function [U, Z, Moment, D2uyDZ2] = gmesh_beam_static(N, NodePos, ...
                                                     NE, ElmConnect, ...
                                                     E, Ix, ...
                                                     LEndNode, REndNode, ...
                                                     Vy, Mx, qy)
%
% Function to solve static FEM cantilevered plane beam problems
% using a general mesh.
%
% Synopsis:
%     [U, Z, Moment, D2uyDz2] =  gmesh_beam_static(N, NodePos, ...
%                                                  NE, ElmConnect, ...
%                                                  E, Ix, ...
%                                                  LEndNode, REndNode, ...
%                                                  Vy, Mx, qy)
%
% Input:
%     N          =   Number of nodes in mesh (numbered 1:N)
%     NodePos    =   Nx1 matrix of nodal positions
%     NE         =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect =   NEx2 matrix of node numbers (for each element)
%     E          =   NEx1 matrix Young's Modulus
%     Ix         =   NEx1 matrix X-sec. moment of inertia
%     LEndNode   =   Node number for Left hand end of beam
%     REndNode   =   Node number for Right hand end of beam
%     Vy         =   Shear force applied to Right hand end of beam
%     Mx         =   Bending moment applied to Right hand end of beam
%     qy         =   Function handle to distributed shear force
%                    per unit length qy(z) as @(z)()
%
% Output:
%     U          =   2Nx1 solution vector from basic FEM solution
%     Z          =   (PlotPtsPerElm*NE+1)x2 matrix with nodal (z) position
%                    and transverse disp.
%     Moment     =   NEx1 vector with element midpoint Moments
%     D2uyDZ2    =   NEx1 vector with element midpoint 2nd deriv. of uy(z)
%
%
% By: Ryan S. Elliott -- Jan. 2015

PlotPtsPerElm = 10;
Z = zeros(PlotPtsPerElm*NE+1,2);
Moment = zeros(NE,1);
D2uyDZ2 = zeros(NE,1);

% Compute global values and ignore M
[M, K, Fqy] = assemble_gmesh_beam_fem(N, NodePos, ...
                                      NE, ElmConnect, ...
                                      E, Ix, ones(N,1), ...
                                      qy, @EqnNumbering);

% Apply global shear and moment boundary conditions
%
%
% determine equation number for REndNode disp.
dEQN = EqnNumbering(REndNode, 1);
% determine equation number for REndNode rot.
rEQN = EqnNumbering(REndNode, 2);

Fy = zeros(size(Fqy));
% Add V to Fy(dEQN)
Fy(dEQN) = Fy(dEQN) + Vy;
% Subtract Mom to F(rEQN) [subtract since Mx is conjugate to -d2uydz2]
Fy(rEQN) = Fy(rEQN) - Mx;

F = Fy + Fqy;

% Apply global displacenet boundary conditions
%
Ftilde = F;
Ktilde = K;

% determine equation number for LEndNode disp.
dEQN = EqnNumbering(LEndNode, 1);
% determine equation number for LEndNode rot.
rEQN = EqnNumbering(LEndNode, 2);

% Set Ftilde(dEQN) to value of node 1 DISPLACEMENT*Ktilde(dEQN, dEQN)
Ftilde(dEQN) = 0.0 * Ktilde(dEQN, dEQN);
%
% Set Ktilde(dEQN, :) to provide the equation:
%    Ktilde(dEQN, dEQN)*U = Ftilde(dEQN)*Ktilde(dEQN, dEQN)
tmp = Ktilde(dEQN, dEQN);
Ktilde(dEQN, :) = 0.0;
Ktilde(dEQN, dEQN) = tmp;
%
% Set Ftilde(rEQN) to value of node 1 DISPLACEMENT*Ktilde(rEQN, rEQN)
Ftilde(rEQN) = 0.0 * Ktilde(rEQN, rEQN);
%
% Set Ktilde(rEQN, :) to provide the equation:
%    Ktilde(rEQN, rEQN)*U = Ftilde(rEQN)*Ktilde(rEQN, rEQN)
tmp = Ktilde(rEQN, rEQN);
Ktilde(rEQN, :) = 0.0;
Ktilde(rEQN, rEQN) = tmp;
%


% find displacements and rotations
U = linsolve(Ktilde,Ftilde);

% FEM Interpolate results for output and ploting
[Z, Moment, D2uyDZ2] = beam_interpolate_results( ...
                         N, NodePos, ...
                         NE, ElmConnect, ...
                         E, Ix, ...
                         LEndNode, REndNode, ...
                         PlotPtsPerElm, U, @EqnNumbering);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = EqnNumbering(Node, NodeDOF)
%   Function to compute global equation number (d) from Node and nodal dof num.

NodeDOFs = 2;
d = NodeDOFs * (Node-1) + NodeDOF;
