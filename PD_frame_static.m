function PD = PD_fram_static(PD)
%
% Function to solve static FEM cantilevered frame problems
% using a general mesh specified with a PD data structure
%
% Synopsis:
%     PD =  PD_frame_static(PD)
%
% Input:
%     PD         =   Matlab structure with the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx1 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx2 matrix of node numbers (for each element)
%        NM         =   Number of Materials Sets (MatsSets)
%        MatsSets   =   NMx1 matrix of Matlab structures containing (at least)
%                       the following fields
%           E          =   Young's Modulus
%           G          =   Shear Modulus
%           A          =   X-sec. area
%           J          =   X-sec. torsional constant
%           Ix         =   X-sec. moment of inertia about x axis
%           Iy         =   X-sec. moment of inertia about y axis
%           Ixy        =   X-sec. cross moment of inertia
%           MassMoment =   X-sec. mass moment of inertia per unit length
%           rho        =   Mass density per unit length
%        ElmMats    =   NEx1 matrix of Material set (for each element)
%        Fx         =   Shear force component applied to Right hand end of frame
%        Fy         =   Shear force component applied to Right hand end of frame
%        Fz         =   Axial force component applied to Right hand end of frame
%        Mx         =   Bending moment component applied to Right hand end of
%                       frame
%        My         =   Bending moment component applied to Right hand end of
%                       frame
%        Mz         =   Torsion moment component applied to Right hand end of
%                       frame
%        qx         =   Function handle to distributed shear force in the x-dir
%                       per unit length qx(z) as @(z)()
%        qy         =   Function handle to distributed shear force in the y-dir
%                       per unit length qy(z) as @(z)()
%        qz         =   Function handle to distributed axial force in the z-dir
%                       per unit length qz(z) as @(z)()
%
% Output:
%     PD            =   The original PD data structure supplied as input with
%                       additional field as listed below
%        U          =   (NTot)x1 Global solution vector
%        Z          =   (PlotPtsPerElm*NE+1)x7 matrix with nodal (z) position
%                       and ux, uy, uz, theta_x, theta_y, theta_z displacements
%                       and rotations.
%        Momment    =   NEx3 matrix with element midpoint Moments: Mx, My, Mz
%        Force      =   NEx3 matrix with element midpoint Forces: Fx, Fy, Fz
%
%
% By: Ryan S. Elliott -- Feb. 2015

PlotPtsPerElm = 10;
Z = zeros(PlotPtsPerElm*PD.NE+1,7);
Moment = zeros(PD.NE,1);
Force = zeros(PD.NE,1);

% Compute global values and ignore M
[M, K, GFq] = assemble_PD_frame_fem(PD);

% Apply global force and moment boundary conditions
%
%
% determine equation numbers for REndNode disps.
dxEQN = PD.EqnNumbering(PD.REndNode, 1);
dyEQN = PD.EqnNumbering(PD.REndNode, 2);
dzEQN = PD.EqnNumbering(PD.REndNode, 3);
% determine equation numbers for REndNode rots.
rxEQN = PD.EqnNumbering(PD.REndNode, 4);
ryEQN = PD.EqnNumbering(PD.REndNode, 5);
rzEQN = PD.EqnNumbering(PD.REndNode, 6);

Fn = zeros(size(GFq));
% Add forces to Fn
Fn(dxEQN) = Fn(dxEQN) + PD.Fx;
Fn(dyEQN) = Fn(dyEQN) + PD.Fy;
Fn(dzEQN) = Fn(dzEQN) + PD.Fz;
% Add/subtract moments to Fn
Fn(rxEQN) = Fn(rxEQN) - PD.Mx;
Fn(ryEQN) = Fn(ryEQN) + PD.My;
Fn(rzEQN) = Fn(rzEQN) + PD.Mz;

F = Fn + GFq;

% Apply global displacenet boundary conditions
%
Ftilde = F;
Ktilde = K;

% determine equation numbers for LEndNode disps.
dxEQN = PD.EqnNumbering(PD.LEndNode, 1);
dyEQN = PD.EqnNumbering(PD.LEndNode, 2);
dzEQN = PD.EqnNumbering(PD.LEndNode, 3);
% determine equation numbers for LEndNode rots.
rxEQN = PD.EqnNumbering(PD.LEndNode, 4);
ryEQN = PD.EqnNumbering(PD.LEndNode, 5);
rzEQN = PD.EqnNumbering(PD.LEndNode, 6);

% Set Ftilde(daEQN) to value of node 1 DISPLACEMENT*Ktilde(daEQN, daEQN)
Ftilde(dxEQN) = 0.0 * Ktilde(dxEQN, dxEQN);
Ftilde(dyEQN) = 0.0 * Ktilde(dyEQN, dyEQN);
Ftilde(dzEQN) = 0.0 * Ktilde(dzEQN, dzEQN);
%
% Set Ktilde(daEQN, :) to provide the equation:
%    Ktilde(daEQN, daEQN)*U = Ftilde(daEQN)*Ktilde(daEQN, daEQN)
xtmp = Ktilde(dxEQN, dxEQN);
ytmp = Ktilde(dyEQN, dyEQN);
ztmp = Ktilde(dzEQN, dzEQN);
Ktilde(dxEQN, :) = 0.0;
Ktilde(dyEQN, :) = 0.0;
Ktilde(dzEQN, :) = 0.0;
Ktilde(dxEQN, dxEQN) = xtmp;
Ktilde(dyEQN, dyEQN) = ytmp;
Ktilde(dzEQN, dzEQN) = ztmp;
%
% Set Ftilde(raEQN) to value of node 1 DISPLACEMENT*Ktilde(raEQN, raEQN)
Ftilde(rxEQN) = 0.0 * Ktilde(rxEQN, rxEQN);
Ftilde(ryEQN) = 0.0 * Ktilde(ryEQN, ryEQN);
Ftilde(rzEQN) = 0.0 * Ktilde(rzEQN, rzEQN);
%
% Set Ktilde(raEQN, :) to provide the equation:
%    Ktilde(raEQN, raEQN)*U = Ftilde(raEQN)*Ktilde(raEQN, raEQN)
xtmp = Ktilde(rxEQN, rxEQN);
ytmp = Ktilde(ryEQN, ryEQN);
ztmp = Ktilde(rzEQN, rzEQN);
Ktilde(rxEQN, :) = 0.0;
Ktilde(ryEQN, :) = 0.0;
Ktilde(rzEQN, :) = 0.0;
Ktilde(rxEQN, rxEQN) = xtmp;
Ktilde(ryEQN, ryEQN) = ytmp;
Ktilde(rzEQN, rzEQN) = ztmp;
%


% find displacements and rotations
PD.U = linsolve(Ktilde,Ftilde);


% @@ NEEDS TO BE FINISHED
% FEM Interpolate results for output and ploting
%PD.PlotPtsPerElm = 10;
%PD = frame_interpolate_results(PD);
