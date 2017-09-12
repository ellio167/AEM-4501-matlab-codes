function [Freq, RawModes, M, K, InterpModes, Zpos] = gmesh_beam_modes(...
                                                      N, NodePos, ...
                                                      NE, ElmConnect, ...
                                                      E, Ix, rho, ...
                                                      LEndNode, ...
                                                      REndNode, ...
                                                      Nmodes)
%
% Function to solve for free vib. modes of cantilevered plane beam problems
%
% Synopsis:
%     [Freq, RawModes, M, K, InterpModes, Zpos]  =  gmesh_beam_modes(...
%                                                    N, NodePos, ...
%                                                    NE, ElmConnect, ...
%                                                    E, Ix, rho, ...
%                                                    LEndNode, REndNode, ...
%                                                    Nmodes)
%
% Input:
%     N          =   Number of nodes in mesh (numbered 1:N)
%     NodePos    =   Nx1 matrix of nodal positions
%     NE         =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect =   NEx2 matrix of node numbers (for each element)
%     E          =   NEx1 matrix Young's Modulus
%     Ix         =   NEx1 matrix X-sec. moment of inertia
%     rho        =   NEx1 matrix Mass density per unit length
%     LEndNode   =   Node number for Left hand end of beam
%     REndNode   =   Node number for Right hand end of beam
%     Nmodes     =   Number of modes to solve for
%
% Output:
%     Freq       =   (Nmodes)x1 vector of vibration frequencies
%     RawModes   =   2NxNmodes matrix of vibration mode disp. and rot.
%     M          =   2Nx2N Mass Matrix
%     K          =   2Nx2N Stiffness Matrix
%     InterpModes=   (PlotPtsPerElm*NE+1)x(Nmodes) matrix of vibration modes
%     Zpos       =   (PlotPtsPerElm*NE+1)x1 matrix of nodal (z) pos.
%
%
% By: Ryan S. Elliott -- Jan. 2015, Feb. 2015

PlotPtsPerElm = 10;
Freq = zeros(Nmodes,1);
Modes = zeros(PlotPtsPerElm*NE+1,Nmodes);

% Compute global values and ignore M
[M, K, GFy] = assemble_gmesh_beam_fem(N, NodePos, ...
                                      NE, ElmConnect, ...
                                      E, Ix, rho, ...
                                      @(z)(0), @EqnNumbering);

% Apply global displacenet boundary conditions
%
% determine equation number for LEndNode disp.
dEQN = EqnNumbering(LEndNode, 1);
% determine equation number for LEndNode rot.
rEQN = EqnNumbering(LEndNode, 2);

% work with deflated matrices to eliminate homogeneouse BCs at left end
if (dEQN < rEQN)
  fstEQN = dEQN;
  sndEQN = rEQN;
else
  fstEQN = rEQN;
  sndEQN = dEQN;
end
Range = [1:fstEQN-1, fstEQN+1:sndEQN-1, sndEQN+1:2*N];
Mhat = M(Range, Range);
Khat = K(Range, Range);


[V,D] = eigs(Khat, Mhat, Nmodes, 'sm');

Freq = sqrt(diag(D));

Mode_i = zeros(size(GFy));
for i = 1:Nmodes
  % mode i
  Mode_i(Range) = V(:,i);
  RawModes(:,i) = Mode_i;
  % FEM Interpolate results for output and ploting
  [Z, Moment, D2uyDZ2] = beam_interpolate_results(...
                           N, NodePos, ...
                           NE, ElmConnect, ...
                           E, Ix, ...
                           LEndNode, REndNode, ...
                           PlotPtsPerElm, Mode_i, @EqnNumbering);
  InterpModes(:,i) = Z(:,2);
end;
Zpos = Z(:,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = EqnNumbering(Node, NodeDOF)
%   function to compute global equation number (d) from Node and nodal dof num.

NodeDOFs = 2;
d = NodeDOFs * (Node-1) + NodeDOF;
