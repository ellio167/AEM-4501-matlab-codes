function PlotMode_elasticity(PD,scl,ModeNum)
% PlotMode_elasticity - Function plots a solved PD structure's mode shape for free vibration problems
%
% Synopsis:
%     PlotMode_elasticity(PD,scl,ModeNum)
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
%     scl           =   Scale factor for modal displacements
%     ModeNum       =   Mode number selected to plot
%
% Output:
%     None
%
% By: Lincoln L. Priebe -- Apr. 2018
%

figh = figure;

x=PD.NodePos(:,1);
y=PD.NodePos(:,2);
ux=scl*PD.Modes(:,1,ModeNum);
uy=scl*PD.Modes(:,2,ModeNum);

pNormal=[x,y];
pDisp=[x+ux,y+uy];
t=PD.ElmConnect;

% create bounds of plot
xMin=min(min(pNormal(:,1)),min(pDisp(:,1)));
xMax=max(max(pNormal(:,1)),max(pDisp(:,1)));
xAbsMax=max(abs(xMin),abs(xMax));
yMin=min(min(pNormal(:,2)),min(pDisp(:,2)));
yMax=max(max(pNormal(:,2)),max(pDisp(:,2)));
yAbsMax=max(abs(yMin),abs(yMax));

simpplot(pNormal,t,'',[1,1,1]);    %plot normal structure
hold on;
simpplot(pDisp,t,'',[1,.3,.3]);    %plot displaced structure

axis([-xAbsMax xAbsMax -yAbsMax yAbsMax])

title(strcat('Mode Number-',int2str(ModeNum)));
legend('Original',strcat('Mode scaled by-',num2str(scl)),'location','bestoutside');
