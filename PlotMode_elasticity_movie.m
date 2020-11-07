function PlotMode_elasticity_movie(PD,scl,ModeNum)
% PlotMode_elasticity_movie - Function creats a movie for a solved PD structure's
%                             mode shape for free vibration problems
%
% Synopsis:
%     PlotMode_elasticity_movie(PD,scl,ModeNum)
%           - NOTICE: This function relies on simpplot.m to run
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   The number of nodes in mesh
%        NE         =   The number of elements in mesh
%        NodePos    =   Nx2 array of the node positions
%        ElmConnect =   NEx3 array of the element connectivity
%        Nmodes     =   Number of modes
%        FreqSq     =   (Nmodes)x1 vector of squared vibration frequencies
%        Modes      =   Nx3xNmodes array of vibration mode displacements
%
%     scl           =   Scale factor of amplitude of modal vibration
%     ModeNum       =   Mode number selected to plot
%
% Output:
%     None
%
% By: Lincoln L. Priebe -- Apr. 2018
%

x=PD.NodePos(:,1);
y=PD.NodePos(:,2);
ux=PD.Modes(:,1,ModeNum); %max deflection in x
uy=PD.Modes(:,2,ModeNum); %max deflection in y

t=PD.ElmConnect;

freq=sqrt(PD.FreqSq(ModeNum));
omega=2*pi*freq;

%select number of cycles to plot
cycles = 5;pts=cycles*20;
if isreal(freq)
    fprintf('--> Frequency selected is %.2d Hz\n',freq)
    dt=1/freq;
    T=cycles*dt;
    time=linspace(0,T,pts);
else
    T=log(10)/imag(omega);
    time=linspace(0,T,pts);
end
fprintf('--> total time of movie is %.2d sec\n\n',T)

pNormal=[x,y];
pDisp=[x+scl*ux,y+scl*uy]; %max dispacement(with scale)

% Find initial plot axis limits
    xMin=min(min(pNormal(:,1)),min(pDisp(:,1)));
    xMax=max(max(pNormal(:,1)),max(pDisp(:,1)));
    xAbsMax=max(abs(xMin),abs(xMax));
    yMin=min(min(pNormal(:,2)),min(pDisp(:,2)));
    yMax=max(max(pNormal(:,2)),max(pDisp(:,2)));
    yAbsMax=max(abs(yMin),abs(yMax));

figure(1)
for i=1:length(time)

    Amplitude=scl*cos(omega*time(i));
    xPosOfT=x+ux*Amplitude;
    yPosOfT=y+uy*Amplitude;


%see if axis should get bigger with time during the first cycle
        if i < pts/cycles
            currentxMax=abs(max(abs(min(xPosOfT)),max(xPosOfT)));
            currentyMax=abs(max(abs(min(yPosOfT)),max(yPosOfT)));
            if currentxMax>xAbsMax
                xAbsMax=currentxMax;
            end
            if currentyMax>yAbsMax
                yAbsMax=currentyMax;
            end
        end

    pDisp=[xPosOfT,yPosOfT];
    simpplot(pNormal,t,'',[1,1,1]);  %continue to plot normal structure
    hold on;
    simpplot(pDisp,t,'',[1,.3,.3]);  %plot changing displaced structure
    axis([-xAbsMax xAbsMax -yAbsMax yAbsMax])
    pause(0);                        %plot as fast as possible
    hold off;
end
