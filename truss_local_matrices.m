function [m, k, T] = truss_local_matrices(Pos1, Pos2, Material)
%
% Function to compute the 6x6 mass and stiffness matricies for a
% 3D truss element with nodes at Pos1 and Pos2 and made of Material

kp = Pos2 - Pos1;
len = norm(kp);
kp = kp/len;

if ((kp(1) == 0) && (kp(2) == 0))
  jp = [0;1;0];
else
  jp = [kp(2); -kp(1); 0];
  jp = jp/norm(jp);
end

ip = cross(jp, kp);

i = [1;0;0];
j = [0;1;0];
k = [0;0;1];
t(1,1) = dot(ip,i);
t(2,1) = dot(ip,j);
t(3,1) = dot(ip,k);
t(1,2) = dot(jp,i);
t(2,2) = dot(jp,j);
t(3,2) = dot(jp,k);
t(1,3) = dot(kp,i);
t(2,3) = dot(kp,j);
t(3,3) = dot(kp,k);

T = [[t, zeros(3,3)];[zeros(3,3), t]];


k = zeros(6,6); % note reuse of k variable!
k(3,3) = 1;
k(3,6) = -1;
k(6,3) = -1;
k(6,6) = 1;

k = (Material.E*Material.A/len)*T*k*T';

m = eye(6);
m = (Material.rho*len/2)*m;
