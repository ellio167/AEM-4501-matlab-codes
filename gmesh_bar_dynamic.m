function sol = gmesh_bar_dynamic(InitialCond, Tfinal, ...
                                 N, NodePos, ...
                                 NE, ElmConnect, ...
                                 E, A, rho, ...
                                 LEndNode, REndNode, ...
                                 Fz, qz)
%
% Function to solve dynamic FEM 1D tension/torsion problems
% using a general mesh.
%
% Synopsis:
%     sol        =  gmesh_bar_dynamic(InitialCond, Tfinal, ...
%                                     N, NodePos, ...
%                                     NE, ElmConnect, ...
%                                     E, A, rho, ...
%                                     LEndNode, REndNode, ...
%                                     EndNode, Fz, qz)
% Input:
%     InitialCond=   2Nx1 matrix with initial displacement (rows 1:N)
%                    and velocities (rows N+1:2N)
%     Tfinal     =   Final time for integration
%     N          =   Number of nodes in mesh (numbered 1:N)
%     NodePos    =   Nx1 matrix of nodal positions
%     NE         =   Number of elements in mesh (numbered 1:NE)
%     ElmConnect =   NEx2 matrix of node numbers (for each element)
%     E          =   NEx1 matrix Young's (tension) or Shear Modulus (torsion)
%     A          =   NEx1 matrix X-sec. area (tension) or polar moment (torsion)
%     rho        =   NEx1 matrix Mass density per unit length (tension)
%                    or centrodial moment of inertia (torsion) per unit length
%     LEndNode   =   Node number for Left hand end of bar
%     REndNode   =   Node number for Right hand end of bar
%     Fz         =   End load force (tension) or torque (torsion) as @(t)()
%     qz         =   Function handle to distributed force per unit length 
%                    or torque per unit length qz(t,z) as @(t,z)()
%
% Output:
%     sol        =   ode45 solution object with "pos" column added
%
%
% By: Ryan S. Elliott -- Jan. 2015

% Compute global (constant) Mass matrix
[M, K, GFqz] = assemble_gmesh_bar_fem(N, NodePos, ...
                                      NE, ElmConnect, ...
                                      E, A, rho, ...
                                      @(z)(qz(0,z)));

% Set Mass matrix for ode solver
MM = [[eye(N), zeros(N,N)];[zeros(N,N), M]];
odeOPTs = odeset('Mass', MM);

sol = ode45(@(t,Z)(get_rhs(t, Z, N, NodePos, ...
                                 NE, ElmConnect, ...
                                 E, A, rho, ...
                                 LEndNode, REndNode, ...
                                 Fz, qz)), ...
            [0, Tfinal], InitialCond, odeOPTs);
sol.pos = NodePos;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = get_rhs(t, Z, N, NodePos, ...
                             NE, ElmConnect, ...
                             E, A, rho, ...
                             LEndNode, REndNode, ...
                             Fz, qz)

[M, K, GFqz] = assemble_gmesh_bar_fem(N, NodePos, ...
                                      NE, ElmConnect, ...
                                      E, A, rho, ...
                                      @(z)(qz(t,z)));

% Apply global force boundary conditions
%
% Add P(t) to GFz(REndNode)
GFz = zeros(N);
GFz(REndNode) = GFz(REndNode) + Fz(t);

F = GFz + GFqz;

FF = [zeros(N,1); F];

% Compute righ hand side of equation of motion
KK = [[zeros(N,N), eye(N)];[-K, zeros(N,N)]];

% Account for disp-type BC's by forcing the associated velocities
% and accelerations to zero
%
FF(LEndNode) = 0.0;     % zero value
KK(LEndNode,:) = 0.0;   % velocity
FF(LEndNode+N) = 0.0;   % zero value
KK(LEndNode+N,:) = 0.0; % acceleration

rhs = FF + KK*Z;
