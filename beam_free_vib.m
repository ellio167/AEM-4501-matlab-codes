close all;
clear all;

Nodes= 11;
L= 5;
E= 90e9;
I= 2.5e-6;
rho= 10;

% create uniform gmesh
[N,NP,NE,EC,E,Ix,R,LEN,REN] = create_gmesh_for_uniform_mesh(Nodes,L,E,I,rho,0);

% get smallest frequency as good inverse time scale
Freq=gmesh_beam_modes(N,NP,NE,EC,E,Ix,R,LEN,REN,2*N-3);
min_fundamental_period=2*pi/Freq(1);
max_fundamental_period=2*pi/Freq(2*N-3);
fundamental_period=max_fundamental_period;

% define loading
Fz=@(t)(0);
qz=@(t,z)(0);

Tfinal = 2*fundamental_period

[U,Z]=gmesh_beam_static(N,NP,NE,EC,E,Ix,LEN,REN, -5000, -15000, @(z)(0));
InitialCond=[U; zeros(2*N,1)];

if exist('beam_free_vib_saved.mat', 'file') == 2
  load('beam_free_vib_saved.mat')
else
  t1=cputime();
  sol=gmesh_beam_dynamic(InitialCond, Tfinal, ...
                         N,NP,NE,EC,E,Ix,R,LEN,REN, ...
                         @(t)(0), @(t)(0), @(t,z)(0));
  t2=cputime();
  Ttot=t2-t1

  save('beam_free_vib_saved.mat');
end

[Freq, RawModes, M, K, InterpModes, Zpos] = ...
   gmesh_beam_modes(N,NP,NE,EC,E,Ix,R,LEN,REN,2*N-3);
for i=1:2*N-3
  omega(i) = Freq(2*N-2-i);
  Phi(i,:) = RawModes(:,2*N-2-i)';
  InterpPhi(i,:) = InterpModes(:,2*N-2-i)';
  A(i) = dot(U, M*Phi(i,:)');
end

uc(:,:,1) = A(1)*cos(omega(1)*sol.x')*Phi(1,:);
ic(:,:,1) = A(1)*cos(omega(1)*sol.x')*InterpPhi(1,:);
for i=2:2*N-3
  uc(:,:,i) = uc(:,:,i-1) + A(i)*cos(omega(i)*sol.x')*Phi(i,:);
  ic(:,:,i) = ic(:,:,i-1) + A(i)*cos(omega(i)*sol.x')*InterpPhi(i,:);
end

figure;
plot(sol.x,sol.y(2*N-3,:),'r-', sol.x,uc(:,2*N-3,2*N-3),'b-');

figure;
hold on;
order = 5;
t = floor(0.1*size(sol.x,2));
plot(sol.pos,sol.interp(t,:),'r-',sol.pos,ic(t,:,order),'b-');
t = floor(0.2*size(sol.x,2));
plot(sol.pos,sol.interp(t,:),'r-',sol.pos,ic(t,:,order),'b-');
t = floor(0.4*size(sol.x,2));
plot(sol.pos,sol.interp(t,:),'r-',sol.pos,ic(t,:,order),'b-');
t = floor(0.7*size(sol.x,2));
plot(sol.pos,sol.interp(t,:),'r-',sol.pos,ic(t,:,order),'b-');
t = floor(0.9*size(sol.x,2));
plot(sol.pos,sol.interp(t,:),'r-',sol.pos,ic(t,:,order),'b-');
