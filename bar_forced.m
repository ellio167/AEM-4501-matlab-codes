L= 5;
N= 20;
E= 70e9;
A= 1e-4;
rho= 2.7*(1000/(100^3));
damp= 1e2;

% get smallest frequency as good inverse time scale
fundamental_freq=bar_modes(N, L, E, A, rho, 1);
fundamental_period=2*pi/fundamental_freq;

% define loading
forcing_freq = 0.8*fundamental_freq;
forcing_period = 2*pi/forcing_freq;
Fz=@(t)((5e4)*cos(forcing_freq*t - pi/2.0));
qz=@(t,z)(0);

Tfinal = 10*forcing_period;

InitialCond=zeros(2*N,1);

sol=bar_dynamic(InitialCond, Tfinal, N, L, E, A, rho, damp, Fz, qz);

% number of time steps in solution
Nt = size(sol.x,2);

% Find limits for plotting
ylimit = max(max(abs(sol.y(1:N,:))));

% repeat numtimes
numtimes = 0;
% run for 10 seconds
delay = 10.0/Nt;

flnm = 'bar_forced_1.gif';
figure(1);
for n=1:Nt
  plot(sol.pos, sol.y(1:N,n), 'o-');
  hold on;
  quiver(sol.pos(N),sol.y(N,n), 0, 0.3*ylimit*Fz(sol.x(n))/(5e4));
  hold off;
  axis normal;
  axis([0,1.05*L, -1.1*ylimit, 1.1*ylimit]);
  xlabel('pos (m)');
  ylabel('disp (m)');
  frame=getframe(1);
  im = frame2im(frame);
  [A,map] = rgb2ind(im,256);
  if n == 1
    imwrite(A,map,flnm,'gif','LoopCount',numtimes,'DelayTime',delay);
  else
    imwrite(A,map,flnm,'gif','WriteMode','append','DelayTime',delay);
  end
end;
close(1);

flnm = 'bar_forced_2.gif';
Scale = 0.1*L/ylimit;
figure(2);
for n=1:Nt
  plot(sol.pos + Scale*sol.y(1:N,n), zeros(N,1), 'o-');
  hold on;
  quiver(sol.pos(N) + Scale*sol.y(N,n), 0, 0.1*L*Fz(sol.x(n))/(5e4), 0);
  hold off;
  axis normal;
  axis([0,1.3*L, -L/10.0, L/10.0]);
  frame=getframe(2);
  im = frame2im(frame);
  [A,map] = rgb2ind(im,256);
  if n == 1
    imwrite(A,map,flnm,'gif','LoopCount',numtimes,'DelayTime',delay);
  else
    imwrite(A,map,flnm,'gif','WriteMode','append','DelayTime',delay);
  end
end;
close(2);
