function PD = PD_torsion(PD, plot_flag)
%
% Function to solve Prandtl torsion problems
%
% Synopsis:
%     PD =  PD_torsion(PD)
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
%        RHS        =   Right-hand side of torsion equation (-2*G*dTheta/dz)
%                       Note: Use a value of -1.0 if you are only interested in
%                             computing the J value.
%     plot_flag  =   1 to generate plots, 0 to skip plots (default)
%
% Output:
%     PD            =   The original PD data structure supplied as input with
%                       additional field as listed below
%        N          =   The number of nodes in the generated mesh
%        NE         =   The number of elements in the generated mesh
%        NodePos    =   Nx2 array of the node positions
%        ElmConnect =   NEx3 array of the element connectivity
%        J          =   The domain's torsional rigidity constant
%        phi        =   Nx1 array of Prandtl stress function values at the nodes
%        ElmCenterX =   NEx1 array of X coordinates for the element center
%        ElmCenterY =   NEx1 array of Y coordinates for the element center
%        ShearStress=   NEx2 array of [Sigma_xz, Sigma_yz] values at elm centers
%
%
% By: Amartya Banerjee, Ryan S. Elliott -- Apr. 2015

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
                                         plot_flag,...
                                         []);

PD.N = size(PD.NodePos,1);
PD.NE = size(PD.ElmConnect,1);
PD.ElmCenterX = zeros(PD.NE,1);
PD.ElmCenterY = zeros(PD.NE,1);
PD.ShearStress = zeros(PD.NE,2);

% Detect the boundary points
BoundaryRange=unique(boundedges(PD.NodePos, PD.ElmConnect));

% Compute global matrices
[K, F] = assemble_PD_torsion_fem(PD);

% Apply global phi value boundary conditions
%
Ftilde = F;
Ktilde = K;

Ktilde(BoundaryRange,:) = 0.0;
for i=BoundaryRange'
  Ftilde(i) = 0.0*K(i,i);
  Ktilde(i,i) = K(i,i);
end

% find phi
PD.phi = linsolve(Ktilde,Ftilde);

% Compute integral of phi: use the one point rule
% Compute various other quantities, too.
IntegralPhi = 0.0;
for i = 1:PD.NE
  % Set local element i node numbers
  Node1 = PD.ElmConnect(i, 1);
  Node2 = PD.ElmConnect(i, 2);
  Node3 = PD.ElmConnect(i, 3);
  % Set local element i node positions
  Pos1 = PD.NodePos(Node1,:);
  Pos2 = PD.NodePos(Node2,:);
  Pos3 = PD.NodePos(Node3,:);
  Pe = [ones(3,1),[Pos1; Pos2; Pos3]]; % 3 by 3 with rows=[1 xcorner ycorner]
  Area = abs(det(Pe))/2.0; % area of triangle e = half of parallelogram area

  % Compute integral
  IntegralPhi = IntegralPhi ...
                + ((Area/3.0) * (PD.phi(Node1)+PD.phi(Node2)+PD.phi(Node3)));

  % compute element center
  center = (Pos1+Pos2+Pos3)/3.0;
  PD.ElmCenterX(i) = center(1);
  PD.ElmCenterY(i) = center(2);

  % compute solution gradient and shear stresses
  grad = linsolve(Pe, [PD.phi(Node1); PD.phi(Node2); PD.phi(Node3)]);
  PD.ShearStress(i,:) = [grad(3), -grad(2)];
end

% Compute J
PD.J = -(4.0/PD.RHS)*IntegralPhi;

% Do some post processing here
if(plot_flag == 1)
    % Plot the FEM approximation phi(x,y)
    figure;
    trisurf(PD.ElmConnect,PD.NodePos(:,1),PD.NodePos(:,2),PD.phi,...
            'edgecolor','none','facecolor','interp');
    view(3);
    %axis equal;
    colorbar;
    %daspect([1 1 1]);
    title('Surface / contour plot of phi')

    % Plot vector field of shear stresses
    figure;
    quiver(PD.ElmCenterX, PD.ElmCenterY, ...
           PD.ShearStress(:,1), PD.ShearStress(:,2));
    view(2);
    axis equal;
    daspect([1 1 1]);
    title('Vector field plot of shear stresses');
end
