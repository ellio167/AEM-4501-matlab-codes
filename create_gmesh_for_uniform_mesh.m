function [N,NodePos,NE,ElmConnect,E,A,rho,LEndNode,REndNode] = ...
   create_gmesh_for_uniform_mesh(Nodes, L, e, a, r, seed)
%
% Function to create general mesh data structures for a uniform
% mesh, with an option to randomize the node numbering.
%
% Synopsis:
%     [N,NodePos,NE,ElmConnect,E,A,rho,LEndNode,REndNode] = ...
%        create_gmesh_for_uniform_mesh(Nodes, L, e, a, r, seed)
%
% Input:
%     Nodes      =   Number of Nodes (number of elements is N-1)
%     L          =   Length of bar
%     e          =   Young's (tension) or Shear Modulus (torsion)
%     a          =   X-section Area (tension) or polar moment (torsion)
%     r          =   Mass density per unit length
%     seed       =   Seed for random number generator.  If positive, nodes will
%                    be numbered randomly along the bar.  Otherwise, nodes will
%                    be numbered sequentially from left to right.
% 
%

N = Nodes;
NE = N - 1;

if (seed > 0.0)
  rand('seed', seed);
  Map = randperm(N);
else
  Map = [1:N];
end;

LEndNode = Map(1);
REndNode = Map(N);

Pos = [0:L/NE:L]';
NodePos(Map,1) = Pos;

j = 1;
for i = 1:NE
  ElmConnect(i,:) = [Map(j), Map(j+1)];
  j = j+1;
end;

E = e*ones(NE,1);
A = a*ones(NE,1);
rho = r*ones(NE,1);
