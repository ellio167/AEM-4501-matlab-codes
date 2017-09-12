function [Z, Moment, D2uyDZ2] = beam_interpolate_results(N, NodePos, ...
                                                         NE, ElmConnect, ...
                                                         E, Ix, ...
                                                         LEndNode, REndNode, ...
                                                         PlotPtsPerElm, U, ...
                                                         EqnNumbering)
%
% Function to create interpolated beam uy(z) curve and Moment and 2nd deriv
%
% Synopsis:
%     [Z, Moment, D2uyDZ2] = beam_interpolate_results(N, NodePos, ...
%                                                     NE, ElmConnect, ...
%                                                     E, Ix, ...
%                                                     LEndNode, REndNode, ...
%                                                     PlotPtsPerElm, U, ...
%                                                     EqnNumbering)
%
% Input:
%     N            =   Number of nodes in mesh (numbered 1:N)
%     NodePos      =   Nx1 matrix of nodal positions
%     NE           =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect   =   NEx2 matrix of node numbers (for each element)
%     E            =   NEx1 matrix Young's Modulus
%     Ix           =   NEx1 matrix X-sec. moment of inertia
%     rho          =   NEx1 matrix Mass density per unit length
%     LEndNode     =   Node number for Left hand end of bar
%     REndNode     =   Node number for Right hand end of bar
%     U            =   2*Nx1 matrix of global solution results
%     EqnNumbering =   Node to global DOF numbering as @(Node, NodeDOF)()
%
% Output:
%     Z            =   (PlotPtsPerElm*NE+1)x2 matrix with nodal (z) position
%                      and transverse disp.
%     Moment       =   NEx1 vector with element midpoint Moments
%     D2uyDZ2      =   NEx1 vector with element midpoint 2nd deriv. of uy(z)
%

% Assume elements of mesh are in order from left to right
DS = 2.0/PlotPtsPerElm;

k = 1;
for i = 1:NE
  % local nodal locations
  Node1 = ElmConnect(i,1);
  Node2 = ElmConnect(i,2);
  G1 = EqnNumbering(Node1,1);
  G2 = EqnNumbering(Node1,2);
  G3 = EqnNumbering(Node2,1);
  G4 = EqnNumbering(Node2,2);
  
  lz(1) = NodePos(Node1);
  lz(2) = NodePos(Node2);
  
  len = lz(2) - lz(1);
  
  % shape functions
  N1 = @(s)(0.25*(1.0-s)^2*(2.0+s));
  N2 = @(s)(0.125*len*(1.0-s)^2*(1.0+s));
  N3 = @(s)(0.25*(1.0+s)^2*(2.0-s));
  N4 = @(s)(0.125*len*(1.0+s)^2*(s-1.0));
  
  % Mapping from master element coordinate (s) to global coordinate (z)
  z = @(s)((lz(1)+lz(2))/2.0 + s*(lz(2)-lz(1))/2.0);
  % Jacobian of mapping
  jac = (lz(2)-lz(1))/2.0;

  for s = -1:DS:1-DS
    Z(k,1) = z(s);
    Z(k,2) = U(G1)*N1(s) + U(G2)*N2(s) + U(G3)*N3(s) + U(G4)*N4(s);
    k = k+1;
  end;

  % 2nd derivative of shape functions
  N1d2 = @(s)(1.5*s);
  N2d2 = @(s)(0.25*len*(3.0*s-1.0));
  N3d2 = @(s)(-0.5*s);
  N4d2 = @(s)(0.25*len*(3.0*s+1.0));
  D2uyDZ2(i) = (U(G1)*N1d2(0) + U(G2)*N2d2(0) ...
               + U(G3)*N3d2(0) + U(G4)*N4d2(0)) / (jac^2);

  Moment(i) = - E(i)*Ix(i)*D2uyDZ2(i);
end;
% Last point
Z(k,1) = z(1);
Z(k,2) = U(G1)*N1(1) + U(G2)*N2(1) + U(G3)*N3(1) + U(G4)*N4(1);
