function PD = frame_interpolate_results(PD)
%
% Function to create interpolated frame deformation curve
%
% Synopsis:
%     PD = frame_interpolate_results(PD)
%
% Input:
%     PD               =   Matlab structure output from PD_frame_static.
%                          Only those fields necessary are listed
%        N             =   Number of nodes in mesh (numbered 1:N)
%        NodePos       =   Nx1 matrix of nodal positions
%        NE            =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect    =   NEx2 matrix of node numbers (for each element)
%        EqnNumbering  =   Node to global DOF numbering as @(Node, NodeDOF)()
%        U             =   (NTot)x1 Global solution vector
%        PlotPtsPerElm =   Number of points to plot within each element
%
% Output:
%     PD               = Input data structure with the following fields added
%        Z             = (PlotPtsPerElm*NE+1)x7 matrix with nodal (z) position
%                        and ux, uy, uz, theta_x, theta_y, theta_z displacements
%                        and rotations.
%        Momment       = NEx3 matrix with element midpoint Moments: Mx, My, Mz
%        Force         = NEx3 matrix with element midpoint Forces: Fx, Fy, Fz
%

% @@@ THIS NEEDS WORK WHERE THE @@@'s ARE

NodeDOFs = 6;

% Assume elements of mesh are in order from left to right
DS = 2.0/PD.PlotPtsPerElm;

k = 1;
for i = 1:PD.NE
  % local nodal locations
  Node1 = PD.ElmConnect(i,1);
  Node2 = PD.ElmConnect(i,2);
  G = zeros(NodeDOFs,2);
  for i = 1:NodeDOFs
    G(i,1) = PD.EqnNumbering(Node1, i);
    G(i,2) = PD.EqnNumbering(Node2, i);
  end;
  
  lz(1) = PD.NodePos(Node1);
  lz(2) = PD.NodePos(Node2);
  
  len = lz(2) - lz(1);
  % shape functions
  BarN1 = @(s)((1.0-s)/2.0);
  BarN2 = @(s)((1.0+s)/2.0);
  BeamN1 = @(s)(0.25*(1.0-s)^2*(2.0+s));
  BeamN2 = @(s)(0.125*len*(1.0-s)^2*(1.0+s));
  BeamN3 = @(s)(0.25*(1.0+s)^2*(2.0-s));
  BeamN4 = @(s)(0.125*len*(1.0+s)^2*(s-1.0));
  % 1st derivative of shape functions
  BarN1ds = @(s)(-0.5);
  BarN2ds = @(s)(0.5);
% @@@  BeamN1ds = @(s)();
% @@@  BeamN2ds = @(s)();
% @@@  BeamN3ds = @(s)();
% @@@  BeamN4ds = @(s)();
  % 2nd derivative of shape functions
  BarN1d2 = @(s)(-0.5);
  BarN2d2 = @(s)(0.5);
  BeamN1d2 = @(s)(1.5*s);
  BeamN2d2 = @(s)(0.25*len*(3.0*s-1.0));
  BeamN3d2 = @(s)(-0.5*s);
  BeamN4d2 = @(s)(0.25*len*(3.0*s+1.0));
  
  
  % Mapping from master element coordinate (s) to global coordinate (z)
  z = @(s)((lz(1)+lz(2))/2.0 + s*(lz(2)-lz(1))/2.0);
  % Jacobian of mapping
  jac = (lz(2)-lz(1))/2.0;

  for s = -1:DS:1-DS
    Z(k,1) = z(s);
    Z(k,2) = U(G(1,1))*BeamN1(s) + U(G(5,1))*BeamN2(s) ...
           + U(G(1,2))*BeamN3(s) + U(G(5,2))*BeamN4(s);
    Z(k,3) = U(G(2,1))*BeamN1(s) + U(G(4,1))*BeamN2(s) ...
           + U(G(2,2))*BeamN3(s) + U(G(4,2))*BeamN4(s);
    Z(k,4) = U(G(3,1)*BarN1(s) + U(G(3,2))*BarN2(s);
% @@@    Z(k,5) =
% @@@    Z(k,6) =
% @@@    Z(k,7) =
    
    k = k+1;
  end;

% @@@   D2uyDZ2(i) = (U(G1)*BeamN1d2(0) + U(G2)*BeamN2d2(0) ...
% @@@                + U(G3)*BeamN3d2(0) + U(G4)*BeamN4d2(0)) / (jac^2);
% @@@ 
% @@@   Moment(i) = - E(i)*Ix(i)*D2uyDZ2(i);
end;
% Last point
s=1;
% @@@ Z(k,1) = z(s);
% @@@ ...
% @@@