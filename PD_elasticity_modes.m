function PD = PD_elasticity_modes(PD, plot_flag, Nmodes)
% PD_elasticity_modes - Function solves 2D elasticity for free vibration
%                       problems
%
% Synopsis:
%     PD =  PD_elasticity_modes(PD, plot_flag)
%           - Solves for all modes of structure
%     PD =  PD_elasticity_modes(PD, plot_flag, Nmodes)
%           - Solves for Nmodes of structure
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        DistFunc   =   Distance Function handle for 2D x-section domain
%                       as @(r)(), with r an Nx2 array of [x,y] values.
%                       The function must be zero on the boundary, negative
%                       inside the domain and positive outside the domain.
%                       For example, an elliptical domain (a=4, b=2):
%                                 @(r)(r(:,1).^2/(4^2) + r(:,2).^2/(2^2) - 1.0)
%        InitEdgeLen=   Initial mesh edge length
%        BBox       =   2x2 matrix [xmin,ymin; xmax,ymax] Domain bounding box
%        Material   =   Matlab structure containing (at least) the following
%                       fields
%           Stiffness  =  3x3 symmetrix 2D stiffness matrix
%           Density    =  Mass density per unit volume
%           Thickness  =  Thickness in third-direction (1.0 for plane-strain)
%     Nmodes        =   Number of modes to solve for
%     plot_flag  =   1 to generate plots, 0 to skip plots (default)
%
% Output:
%     PD            =   The original PD data structure supplied as input with
%                       additional field as listed below
%        N          =   The number of nodes in the generated mesh
%        NE         =   The number of elements in the generated mesh
%        NodePos    =   Nx2 array of the node positions
%        ElmConnect =   NEx3 array of the element connectivity
%        Nmodes     =   Number of modes
%        FreqSq     =   (Nmodes)x1 vector of squared vibration frequencies
%        Modes      =   Nx3xNmodes array of vibration mode displacements

%
% By: Amartya Banerjee, Ryan S. Elliott, Lincoln L. Priebe -- Apr. 2015, Apr. 2018

%
% The mesh generation is performed by the routines from the DistMesh suite of
% programs (developed by Persson and Strang)
%
% Please visit http://persson.berkeley.edu/distmesh/ for more details on
% the DistMesh suite of programs.
%
% The program PD_torsion.m is partially based on the FEM solver described by
% John Coady https://www.particleincell.com/2012/matlab-fem/
%

% set default value for plot_flag if not provided
switch nargin
  case 1
    plot_flag = 0;
end

% create mesh:
[PD.NodePos, PD.ElmConnect] = distmesh2d(PD.DistFunc,...
                                         @huniform,...
                                         PD.InitEdgeLen,...
                                         PD.BBox,...
                                         [],...
                                         plot_flag);

PD.N = size(PD.NodePos,1);
PD.NE = size(PD.ElmConnect,1);



% Detect the boundary points
BoundaryRange=unique(boundedges(PD.NodePos, PD.ElmConnect));

% Compute global matrices
PD.EqnNumbering = @(Node, NodeDOF)(2*(Node-1) + NodeDOF);
[M, K] = assemble_PD_elasticity_fem(PD);


%% SECTION 2 - Make changes to K to specify first 3 modes as translation in x,y and rot.

T = zeros((PD.N*2),1);
Tx=T;Ty=T;
Tx(1:2:length(T)) = 1; Tx=Tx/sqrt(Tx'*M*Tx);  % normalize through the mass matrix
Ty(2:2:length(T)) = 1; Ty=Ty/sqrt(Ty'*M*Ty);  % normalize through the mass matrix
theta = .01;
t=[0 theta; -theta 0];
R=[];    %initialize R matrix
for i=1:(PD.N)
    u=[PD.NodePos(i,1);PD.NodePos(i,2)];   %x,y coordinate of each node
    r=t*u;
    R=cat(1,R,r);
end
R=R-Tx*dot(R,M*Tx)-Ty*dot(R,M*Ty); %force R to be M-orthogonal to Tx,Ty
R=R/sqrt(R'*M*R);

Ax=(0.01)*(M*Tx)*(M*Tx)';
Ay=(0.02)*(M*Ty)*(M*Ty)';
Ar=(0.03)*(M*R)*(M*R)';

K = K + Ax + Ay + Ar;



% Unconstrained free vibration has no BCs to apply
%

% find modes
if nargin < 3
    [V, D] = eigs(K, M, PD.N*2, 'smallestreal');  %solve for all modes
    Nmodes=size(D,1);
else
    [V, D] = eigs(K, M, Nmodes, 'smallestreal');  %solve for Nmodes
end

PD.Nmodes = Nmodes;
PD.FreqSq = zeros(Nmodes, 1);
PD.Modes = zeros(PD.N, 2, Nmodes);

PD.FreqSq = diag(D);
for k=1:Nmodes
  Mode = V(:,k);
  for i=1:PD.N
    for j=1:2
      EQN=PD.EqnNumbering(i,j);
      PD.Modes(i,j,k) = Mode(EQN);
    end
  end
end

PD=rmfield(PD, 'EqnNumbering');
