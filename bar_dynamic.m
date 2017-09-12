function sol = bar_dynamic(InitialCond, Tfinal, N, L, E, A, rho, damp, Fz, qz)
%
% Function to solve dynamic FEM 1D tension/torsion problems
% using a uniform mesh.
%
% Synopsis:
%     sol        =  bar_dynamic(InitialCond, Tfinal, N, L, E, A, rho, damp ...
%                               Fz, qz)
%
% Input:
%     InitialCond=   2Nx1 matrix with initial displacement (rows 1:N)
%                    and velocities (rows N+1:2N)
%     Tfinal     =   Final time for integration
%     N          =   Number of Nodes (number of elements is N-1)
%     L          =   Length of bar
%     E          =   Young's (tension) or Shear Modulus (torsion)
%     A          =   X-section Area (tension) or polar moment (torsion)
%     rho        =   Mass density per unit length (tension)
%     damp       =   Rayleigh Damping coefficient for each node
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
[M, K, GFqz] = assemble_bar_fem(N, L, E, A, rho, @(z)(qz(0,z)));

% Set Mass matrix for ode solver
MM = [[eye(N), zeros(N,N)];[zeros(N,N), M]];
odeOPTs = odeset('Mass', MM);

sol = ode45(@(t,Z)(get_rhs(t, Z, N, L, E, A, rho, damp, Fz, qz)), ...
            [0, Tfinal], InitialCond, odeOPTs);
sol.pos = [0:N-1]'*L/(N-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = get_rhs(t, Z, N, L, E, A, rho, damp, Fz, qz)

[M, K, GFqz] = assemble_bar_fem(N, L, E, A, rho, @(z)(qz(t,z)));

% Create diagonal Rayleigh damping matrix
C = damp*eye(N);

% Add P(t) to GFz to apply global force BCs
GFz = zeros(N,1);
GFz(N) = Fz(t);

F = GFz + GFqz;

FF = [zeros(N,1); F];

% Compute righ hand side of equation of motion
KK = [[zeros(N,N), eye(N)];[-K, -C]];

% Account for disp-type BC's by forcing the associated velocities
% and accelerations to zero
FF(1) = 0.0;     % zero value
KK(1,:) = 0.0;   % velocity
FF(1+N) = 0.0;   % zero value
KK(1+N,:) = 0.0; % acceleration

rhs = FF + KK*Z;
