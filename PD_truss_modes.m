function PD = PD_truss_modes(PD, Nmodes)
%
% Function to solve free vib. modes of 3D truss problems
% using a general mesh specified with a PD data structure
%
% Synopsis:
%     PD =  PD_truss_modes(PD, Nmodes)
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
%           rho        =   Mass density per unit length
%        ElmMats    =   NEx1 matrix of Material set (for each element)
%        BCType     =   Nx3 matrix of 0's (force-type) and 1's (disp-type)
%     Nmodes     =   Number of modes to solve for
%
% Output:
%     PD            =   The original PD data structure supplied as input with
%                       additional field as listed below
%        Nmodes     =   Number of modes
%        Freq       =   (Nmodes)x1 vector of vibration frequencies
%        Modes      =   Nx3xNmodes array of vibration mode displacements
%
%
% By: Ryan S. Elliott -- Feb. 2015

PD.Nmodes = Nmodes;
PD.Freq = zeros(Nmodes,1);
PD.Modes = zeros(PD.N,3,Nmodes);

% check for duplicate bars
[~,unique_bars,~] = unique(sort(PD.ElmConnect,2), 'rows');
if (size(unique_bars) ~= PD.NE)
  fprintf(1, 'NOTE: Duplicate bars are present!\n');
  fprintf(1, '  Duplicate bar numbers are:\n');
  duplicate_bars = setdiff([1:PD.NE],unique_bars);
  fprintf(1, '\t%i\n', sort(duplicate_bars));
  clear duplicate_bars;
end
clear unique_bars;

% Compute global values and ignore M
PD.EqnNumbering = @(Node, NodeDOF)(3*(Node-1) + NodeDOF);
[M, K] = assemble_PD_truss_fem(PD);

% Apply global (zero) displacement boundary conditions
%
Range = [];
r = 1;
for i=1:PD.N
  for j=1:3
    if (PD.BCType(i,j) == 0)
      EQN = PD.EqnNumbering(i,j);
      Range(r) = EQN;
      r = r+1;
    end
  end
end
%
Mhat = M(Range, Range);
Khat = K(Range, Range);

% find modes
[V,D] = eigs(Khat, Mhat, Nmodes, -0.1);

PD.Freq = sqrt(diag(D));
for k=1:Nmodes
  Mode(1:3*PD.N) = 0.0;
  Mode(Range) = V(:,k);
  for i=1:PD.N
    for j=1:3
      EQN=PD.EqnNumbering(i,j);
      PD.Modes(i,j,k) = Mode(EQN);
    end
  end
end

PD=rmfield(PD, 'EqnNumbering');
