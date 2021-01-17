function PD = PD_truss_static(PD)
%
% Function to solve static FEM 3D truss problems
% using a general mesh specified with a PD data structure
%
% Synopsis:
%     PD =  PD_truss_static(PD)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx3 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx2 matrix of node numbers (for each element)
%        NM         =   Number of Materials Sets (MatsSets)
%        MatsSets   =   NMx1 matrix of Matlab structures containing (at least)
%                       the following fields
%           E          =   Young's Modulus
%           A          =   X-sec. area
%        ElmMats    =   NEx1 matrix of Material set (for each element)
%        BCType     =   Nx3 matrix of 0's (force-type) and 1's (disp-type)
%        BCVal      =   Nx3 matrix of BC signed magnitudes
%
% Output:
%     PD            =   The original PD data structure supplied as input with
%                       additional field as listed below
%        U          =   Nx3 Global displacement solution vector
%        R          =   Nx3 Reaction force vector
%        ElmForce   =   NEx1 Force in each bar
%        ElmStress  =   NEx1 Stress in each bar
%
%
% By: Ryan S. Elliott -- Feb. 2015

PD.U = zeros(PD.N,3);
PD.R = zeros(PD.N,3);
PD.ElmForce = zeros(PD.NE,1);
PD.ElmStress = zeros(PD.NE,1);

% create copy of input
PDcpy = PD;

% check for duplicate bars
[~,unique_bars,~] = unique(sort(PDcpy.ElmConnect,2), 'rows');
if (size(unique_bars) ~= PDcpy.NE)
  fprintf(1, 'NOTE: Duplicate bars are present!\n');
  fprintf(1, '  Duplicate bar numbers are:\n');
  duplicate_bars = setdiff([1:PDcpy.NE],unique_bars);
  fprintf(1, '\t%i\n', sort(duplicate_bars));
  clear duplicate_bars;
end
clear unique_bars;

% overwrite density with arbitrary value
for i=1:PDcpy.NM
  PDcpy.MatsSets(i).rho = 1.0;
end

% Compute global values and ignore M
PDcpy.EqnNumbering = @(Node, NodeDOF)(3*(Node-1) + NodeDOF);
[M, K] = assemble_PD_truss_fem(PDcpy);

% Apply global force boundary conditions
%
F = zeros(3*PDcpy.N);
for i=1:PDcpy.N
  for j=1:3
    if (PDcpy.BCType(i,j) == 0)
      F(PDcpy.EqnNumbering(i,j)) = PDcpy.BCVal(i,j);
    end
  end
end

% Apply global displacenet boundary conditions
%
Ftilde = F;
Ktilde = K;

for i=1:PDcpy.N
  for j=1:3
    if (PDcpy.BCType(i,j) == 1)
      EQN = PDcpy.EqnNumbering(i,j);
      tmp = Ktilde(EQN,EQN);
      if (tmp == 0)
        tmp = max(diag(Ktilde));
      end;
      Ftilde(EQN) = PDcpy.BCVal(i,j)*tmp;
      Ktilde(EQN,:) = 0.0;
      Ktilde(EQN,EQN) = tmp;
    end
  end
end

% find displacements
U = linsolve(Ktilde,Ftilde);
% find reactions
R = -F + K*U;

% update original PD with results
for i=1:PDcpy.N
  for j=1:3
    EQN=PDcpy.EqnNumbering(i,j);
    PD.U(i,j) = U(EQN);
    PD.R(i,j) = R(EQN);
  end
end

% compute force and stress
for i=1:PDcpy.NE
  [m, k, T] = truss_local_matrices(PDcpy.NodePos(PDcpy.ElmConnect(i,1),:)', ...
                                   PDcpy.NodePos(PDcpy.ElmConnect(i,2),:)', ...
                                   PDcpy.MatsSets(PDcpy.ElmMats(i)));
  f = T'*k*[PD.U(PDcpy.ElmConnect(i,1),:)'; PD.U(PDcpy.ElmConnect(i,2),:)'];
  PD.ElmForce(i) = f(6);
  PD.ElmStress(i) = PD.ElmForce(i)/PDcpy.MatsSets(PDcpy.ElmMats(i)).A;
end
