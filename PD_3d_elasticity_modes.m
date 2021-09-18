function PD = PD_3d_elasticity_modes(PD, plot_flag, Nmodes)
% PD_3d_elasticity_modes - Function solves 3D elasticity for free vibration
%                          problems
%
% Synopsis:
%     PD =  PD_3d_elasticity_modes(PD, plot_flag)
%           - Solves for all modes of structure
%     PD =  PD_3d_elasticity_modes(PD, plot_flag, Nmodes)
%           - Solves for Nmodes of structure
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        DistFunc   =   Distance Function handle for 3D domain
%                       as @(r)(), with r an Nx3 array of [x,y,z] values.
%                       The function must be zero on the boundary, negative
%                       inside the domain and positive outside the domain.
%                       For example, an ellipsoidal domain (a=4, b=2, c=3):
%                          @(r)(r(:,1).^2/(4^2) + r(:,2).^2/(2^2) + r(:,3).^2/(3^2) - 1.0)
%        InitEdgeLen=   Initial mesh edge length
%        BBox       =   2x3 matrix [xmin,ymin,zmin; xmax,ymax,zmax]
%                       bounding box
%        Material   =   Matlab structure containing (at least) the following
%                       fields
%           Stiffness  =  6x6 symmetrix 3D stiffness matrix
%           Density    =  Mass density per unit volume
%     Nmodes        =   Number of modes to solve for
%     plot_flag  =   1 to generate plots, 0 to skip plots (default)
%
% Output:
%     PD            =   The original PD data structure supplied as input with
%                       additional field as listed below
%        N          =   The number of nodes in the generated mesh
%        NE         =   The number of elements in the generated mesh
%        NodePos    =   Nx3 array of the node positions
%        ElmConnect =   NEx4 array of the element connectivity
%        Nmodes     =   Number of modes
%        FreqSq     =   (Nmodes)x1 vector of squared vibration frequencies
%        Modes      =   Nx3xNmodes array of vibration mode displacements

%
% By: Ryan S. Elliott, Amartya Banerjee, Lincoln L. Priebe -- Apr. 2015, Apr. 2018, Apr. 2021

%
% The mesh generation is performed by the routines from the DistMesh suite of
% programs (developed by Persson and Strang)
%
% Please visit http://persson.berkeley.edu/distmesh/ for more details on
% the DistMesh suite of programs.
%
% The program PD_3d_elasticity_modes.m is partially based on the FEM solver
% described by John Coady https://www.particleincell.com/2012/matlab-fem/
%

% set default value for plot_flag if not provided
switch nargin
  case 1
    plot_flag = 0;
end

% create mesh:
[PD.NodePos, PD.ElmConnect] = distmeshnd(PD.DistFunc,...
                                         @huniform,...
                                         PD.InitEdgeLen,...
                                         PD.BBox,...
                                         []);
PD.N = size(PD.NodePos,1);
PD.NE = size(PD.ElmConnect,1);


% Detect the boundary points
BoundaryRange=unique(boundedges(PD.NodePos, PD.ElmConnect));

% Compute global matrices
PD.EqnNumbering = @(Node, NodeDOF)(3*(Node-1) + NodeDOF);
[M, K] = assemble_PD_3d_elasticity_fem(PD);

PD.K=K;
PD.M=M;

%% SECTION 2 - Make changes to K to specify first 6 modes as translation in x,y,z and rots.

T = zeros((PD.N*3),1);
Tx=T;Ty=T;Tz=T;
Tx(1:3:length(T)) = 1; Tx=Tx/sqrt(Tx'*M*Tx);  % normalize through the mass matrix
Ty(2:3:length(T)) = 1; Ty=Ty/sqrt(Ty'*M*Ty);  % normalize through the mass matrix
Tz(3:3:length(T)) = 1; Tz=Tz/sqrt(Tz'*M*Tz);  % normalize through the mass matrix
theta = .001;
rot_z=[0 theta 0; -theta 0 0; 0 0 0];
rot_x=[0 0 0; 0 0 theta; 0 -theta 0];
rot_y=[0 0 theta; 0 0 0; -theta 0 0];
%initialize R matries
R_z=[];
R_x=[];
R_y=[];
for i=1:(PD.N)
    r=[PD.NodePos(i,1);PD.NodePos(i,2);PD.NodePos(i,3)];   %x,y,z coordinate of each node
    u_z=rot_z*r;
    u_x=rot_x*r;
    u_y=rot_y*r;
    R_z=cat(1,R_z,u_z);
    R_x=cat(1,R_x,u_x);
    R_y=cat(1,R_y,u_y);
end
% force R_z to be M-orthogonal to Tx,Ty,Tz
R_z=R_z-Tx*dot(R_z,M*Tx)-Ty*dot(R_z,M*Ty)-Tz*dot(R_z,M*Tz);
R_z=R_z/sqrt(R_z'*M*R_z);
% force R_x to be M-orthogonal to Tx,Ty,Tz,R_z
R_x=R_x-Tx*dot(R_x,M*Tx)-Ty*dot(R_x,M*Ty)-Tz*dot(R_x,M*Tz)-R_z*dot(R_x,M*R_z);
R_x=R_x/sqrt(R_x'*M*R_x);
% force R_x to be M-orthogonal to Tx,Ty,Tz,R_z,R_x
R_y=R_y-Tx*dot(R_y,M*Tx)-Ty*dot(R_y,M*Ty)-Tz*dot(R_y,M*Tz)-R_z*dot(R_y,M*R_z)-R_y*dot(R_y,M*R_x);
R_y=R_y/sqrt(R_y'*M*R_y);

Ax=(0.01)*(M*Tx)*(M*Tx)';
Ay=(0.02)*(M*Ty)*(M*Ty)';
Az=(0.03)*(M*Tz)*(M*Tz)';
Arz=(0.04)*(M*R_z)*(M*R_z)';
Arx=(0.05)*(M*R_x)*(M*R_x)';
Ary=(0.06)*(M*R_y)*(M*R_y)';

K = K + Ax + Ay + Az + Arz + Arx + Ary;



% Unconstrained free vibration has no BCs to apply
%

% find modes
if nargin < 3
    [V, D] = eigs(K, M, PD.N*3, 'smallestreal');  %solve for all modes
    Nmodes=size(D,1);
else
    opts.tol = 1e-8;  % put in to fixup convergence, not sure why it should be needed
    opts.spdB = true;
    opts.disp = true;
    [V, D] = eigs((K+K')/2.0, (M+M')/2.0, Nmodes, 'smallestreal',opts);  %solve for Nmodes
end

PD.Nmodes = Nmodes;
PD.FreqSq = zeros(Nmodes, 1);
PD.Modes = zeros(PD.N, 3, Nmodes);

PD.FreqSq = diag(D);
for k=1:Nmodes
  Mode = V(:,k);
  for i=1:PD.N
    for j=1:3
      EQN=PD.EqnNumbering(i,j);
      PD.Modes(i,j,k) = Mode(EQN);
    end
  end
end

PD=rmfield(PD, 'EqnNumbering');
