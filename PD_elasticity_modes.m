function PD = PD_elasticity_modes(PD, Nmodes, plot_flag)
%
% Function to solve 2D elasticity free vibration problems
%
% Synopsis:
%     PD =  PD_elasticity_modes(PD, Nmodes, plot_flag)
%
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
% By: Amartya Banerjee, Ryan S. Elliott -- Apr. 2015, Apr. 2018

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

PD.Nmodes = Nmodes;
PD.FreqSq = zeros(Nmodes, 1);
PD.Modes = zeros(PD.N, 2, Nmodes);

% Detect the boundary points
BoundaryRange=unique(boundedges(PD.NodePos, PD.ElmConnect));

% Compute global matrices
PD.EqnNumbering = @(Node, NodeDOF)(2*(Node-1) + NodeDOF);
[M, K] = assemble_PD_elasticity_fem(PD);

% Unconstrained free vibration has no BCs to apply
%

% find modes
[V, D] = eigs(K, M, Nmodes, 'SA');

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
