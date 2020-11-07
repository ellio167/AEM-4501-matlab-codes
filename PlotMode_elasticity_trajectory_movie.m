function PlotMode_elasticity_trajectory_movie(PD,TimeScale,scl)
% PlotMode_elasticity_trajectory_movie - Function creats a movie for a solved PD structure's
%                                        dynamic response to initial displacements
%
% TODOs:
%  * Change to take initial displacement field as input
%    - compute modal decomposition to get the "Avec" amplitudes
%  * Probably get rid of "scl" and just let it be handled by the initial displacement values.
%
% Synopsis:
%     PlotMode_elasticity_trajectory_movie(PD,TimeScale,scl)
%           - NOTICE: This function relies on simpplot.m to run
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   The number of nodes in mesh
%        NE         =   The number of elements in mesh
%        NodePos    =   Nx2 array of the node positions
%        ElmConnect =   NEx3 array of the element connectivity
%        Nmodes     =   N*2 , number of modes of entire structure
%        FreqSq     =   Nmodesx1 vector of squared vibration frequencies
%        Modes      =   Nx3xNmodes array of vibration mode displacements
%        Avec       =   Nmodesx1 vector of modal amplitudes
%     TimeScale     =   Vibration frequency of which plotting is based on
%                       = 'high' for highest frequency (slowest movie)
%                       = 'low' for lowest frequency (lowest movie)
%                       = 'MaxA' for most dominant frequency
%     scl           =   Scale factor of amplitudes of modal vibration
%
% Output:
%     None
%
% By: Lincoln L. Priebe -- Apr. 2018
%



x=PD.NodePos(:,1);
y=PD.NodePos(:,2);
ux=PD.Modes(:,1,(1:1:PD.Nmodes)); %max deflection in x
uy=PD.Modes(:,2,(1:1:PD.Nmodes)); %max deflection in y

t=PD.ElmConnect;


omega=sqrt(PD.FreqSq(1:1:PD.Nmodes));

% Select cycles to plot / choose time scale of movie
cycles = 5;pts=cycles*20;
if isreal(sqrt(PD.FreqSq(1)))
    if strcmp(TimeScale,'high')
       freqPlot=omega(PD.Nmodes)/2/pi;  %movie based on lowest freq
       fprintf('\n--> Highest frequency is %.2d Hz\n',freqPlot)
    elseif strcmp(TimeScale,'low')
       freqPlot=omega(1)/2/pi;          %movie based on highest freq
       fprintf('\n--> Lowest frequency is %.2d Hz\n',freqPlot)
    elseif strcmp(TimeScale,'MaxA')
       [~, MaxAIndex]=max(abs(PD.Avec));
       freqPlot=omega(MaxAIndex)/2/pi;  %movie based on most dominant freq
       fprintf('\n--> Most Dominant frequency is %.2d Hz\n',freqPlot)
    else
       error('Either enter "high", "low", or "MaxA" for Timescale');
    end

    dt=1/freqPlot;
    T=cycles*dt;
    time=linspace(0,T,pts);
    pDisp=[x+scl*ux,y+scl*uy];
    fprintf('--> total time of movie is %.2d sec\n\n',T)
else
    dt=log(10)/imag(omega(1));     %used omega based on lowest omega
    T=10*dt;
    time=linspace(0,T,10*pts);
    pDisp=[x+100*ux,y+100*uy];
    fprintf('--> total time of movie is %.2d sec\n\n',T)
end

pNormal=[x,y];
% Use to find initial plot axis limits
    xMin=min(min(pNormal(:,1)),min(pDisp(:,1)));
    xMax=max(max(pNormal(:,1)),max(pDisp(:,1)));
    xAbsMax=max(abs(xMin),abs(xMax));
    yMin=min(min(pNormal(:,2)),min(pDisp(:,2)));
    yMax=max(max(pNormal(:,2)),max(pDisp(:,2)));
    yAbsMax=max(abs(yMin),abs(yMax));

% Initialize
Amplitude=zeros(PD.Nmodes,1);
A=Amplitude;
xPos=zeros(PD.N,PD.Nmodes);
yPos=xPos;
xPosOfT=zeros(PD.N,1);
yPosOfT=zeros(PD.N,1);

% Create Movie
figure(1)
for i=1:length(time)
    for j = 1:PD.Nmodes
        A(j)=PD.Avec(j);
        Amplitude(j)=A(j)*scl*cos(omega(j)*time(i));
        xPos(:,j)=x+ux(:,:,j)*Amplitude(j);          %each collumn is xPos at jth mode number
        yPos(:,j)=y+uy(:,:,j)*Amplitude(j);          %each collumn is yPos at jth mode number
    end
    for k = 1:PD.N                                   %add up nodal displacements from every mode number
            xPosOfT(k)=sum(xPos(k,:));
            yPosOfT(k)=sum(yPos(k,:));
    end

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
            %pause()     %allows for visualization of initial displacements
        end

    pDisp=[xPosOfT,yPosOfT];
    simpplot(pNormal+15,t,'',[1,1,1]);  %continue to plot normal structure
    hold on;
    simpplot(pDisp,t,'',[1,.3,.3]);     %plot changing displaced structure
    axis([-xAbsMax xAbsMax -yAbsMax yAbsMax])
    pause(0);                           %plot as fast as possible
    hold off;
end
