function PD = PD_torsion_poly(PD, plot_flag)
%
% Function to solve Prandtl torsion problems
%
% Synopsis:
%     PD =  PD_torsion(PD)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        VertexList =   (N+1)x2 matrix with counter-clockwise cycle of polygon
%                       vertex (x,y) locations.  First and last rows must be
%                       identical.
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
%        DistFunc   =   Distance Function handle for 2D x-section domain
%                       as @(r)(), with r an Nx2 array of [x,y] values.
%                       The function must be zero on the boundary, negative
%                       inside the domain and positive outside the domain.
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

% set default value for plot_flag if not provided
switch nargin
  case 1
    plot_flag = 0;
end

% set DistFunc
PD.DistFunc = @(p)(dpoly(p,PD.VertexList));

% call PD_torsion
PD = PD_torsion(PD, plot_flag);