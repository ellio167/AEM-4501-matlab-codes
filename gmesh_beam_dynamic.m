function sol = gmesh_beam_dynamic(InitialCond, Tfinal, ...
                                  N, NodePos, ...
                                  NE, ElmConnect, ...
                                  E, Ix, rho, ...
                                  LEndNode, REndNode, ...
                                  Vy, Mx, qy)
%
% Function to solve dynamic FEM cantilevered plan bema problems
% using a general mesh.
%
% Synopsis:
%     sol        =  gmesh_beam_dynamic(InitialCond, Tfinal, ...
%                                      N, NodePos, ...
%                                      NE, ElmConnect, ...
%                                      E, Ix, rho, ...
%                                      LEndNode, REndNode, ...
%                                      Vy, Mx, qy)
% Input:
%     InitialCond=   2Nx1 matrix with initial displacement (rows 1:2N)
%                    and velocities (rows 2N+1:4N)
%     Tfinal     =   Final time for integration
%     N          =   Number of nodes in mesh (numbered 1:N)
%     NodePos    =   Nx1 matrix of nodal positions
%     NE         =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect =   NEx2 matrix of node numbers (for each element)
%     E          =   NEx1 matrix Young's (tension) or Shear Modulus (torsion)
%     Ix         =   NEx1 matrix X-sec. moment of inertia
%     rho        =   NEx1 matrix Mass density per unit length
%     LEndNode   =   Node number for Left hand end of beam
%     REndNode   =   Node number for Right hand end of beam
%     Vy         =   Shear force applied to Right hand end of beam as @(t)()
%     Mx         =   Bending moment applied to Right hand end of beam as @(t)()
%     qy         =   Function handle to distributed shear force
%                    per unit length qy(t,z) as @(t,z)()
%
% Output:
%     sol        =   ode45 solution object with "pos" column and "interp"
%                    matrix added.  pos is the interpolaged positions associated
%                    with the interpolated displacements.  interp
%                    contains the interpolated solution corresponding to y
%
%
% By: Ryan S. Elliott -- Jan. 2015

NodeDOFs = 2;

NTot = NodeDOFs * N;

% Compute global (constant) Mass matrix
[M, K, GFqy] = assemble_gmesh_beam_fem(N, NodePos, ...
                                       NE, ElmConnect, ...
                                       E, Ix, rho, ...
                                       @(z)(qy(0,z)), @EqnNumbering);

% Set Mass matrix for ode solver
MM = [[eye(NTot), zeros(NTot,NTot)];[zeros(NTot,NTot), M]];
odeOPTs = odeset('Mass', MM);

sol = ode45(@(t,Z)(get_rhs(t, Z, N, NodePos, ...
                                    NE, ElmConnect, ...
                                    E, Ix, rho, ...
                                    LEndNode, REndNode, ...
                                    Vy, Mx, qy)), ...
            [0, Tfinal], InitialCond, odeOPTs);

% Post process results for plotting
PlotPtsPerElm = 10;

% Create pos
[Z, Moment, D2uyDZ2] = beam_interpolate_results(N, NodePos, NE, ElmConnect, ...
                                                E, Ix, LEndNode, REndNode, ...
                                                PlotPtsPerElm, InitialCond, ...
                                                @EqnNumbering);
sol.pos = Z(:,1);

% Create interp
for i = 1:size(sol.y,2)
  [Z, Moment, D2uyDZ2] = beam_interpolate_results(N, NodePos, NE, ElmConnect,...
                                                  E, Ix, LEndNode, REndNode, ...
                                                  PlotPtsPerElm, sol.y(:,i), ...
                                                  @EqnNumbering);
  sol.interp(i,:) = Z(:,2)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = get_rhs(t, Z, N, NodePos, ...
                                NE, ElmConnect, ...
                                E, Ix, rho, ...
                                LEndNode, REndNode, ...
                                Vy, Mx, qy)

NodeDOFs = 2;

NTot = NodeDOFs * N;

[M, K, GFqy] = assemble_gmesh_beam_fem(N, NodePos, ...
                                       NE, ElmConnect, ...
                                       E, Ix, rho, ...
                                       @(z)(qy(t,z)), @EqnNumbering);

% Apply global shear and moment boundary conditions
%
GFy = zeros(size(GFqy));

%
% determine equation number for REndNode disp.
dEQN = EqnNumbering(REndNode, 1);
% determine equation number for REndNode rot.
rEQN = EqnNumbering(REndNode, 2);

% Set V in GFy(dEQN)
GFy(dEQN) = Vy(t);
% Set Mom in GFy(rEQN)
GFy(rEQN) = - Mx(t);

F = GFy + GFqy;

FF = [zeros(NTot,1); F];

% Compute righ hand side of equation of motion
KK = [[zeros(NTot,NTot), eye(NTot)];[-K, zeros(NTot,NTot)]];

% Account for disp-type BC's by forcing the associated velocities
% and accelerations to zero
%
% determine equation number for LEndNode disp.
dEQN = EqnNumbering(LEndNode, 1);
% determine equation number for LEndNode rot.
rEQN = EqnNumbering(LEndNode, 2);
%
FF(dEQN) = 0.0;        % zero value
KK(dEQN,:) = 0.0;      % velocity
FF(dEQN+NTot) = 0.0;   % zero value
KK(dEQN+NTot,:) = 0.0; % acceleration

FF(rEQN) = 0.0;        % zero value
KK(rEQN,:) = 0.0;      % velocity
FF(rEQN+NTot) = 0.0;   % zero value
KK(rEQN+NTot,:) = 0.0; % acceleration

rhs = FF + KK*Z;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = EqnNumbering(Node, NodeDOF)
%   Function to compute global equation number (d) from Node and nodal dof num.

NodeDOFs = 2;
d = NodeDOFs * (Node-1) + NodeDOF;
