function PlotTrussMode(PD, scl, ModeNumber, plotlab)
%
% PlotTrussMode   function to plot a truss normal mode graphically
%
% Synopsis:
%     PlotTrussMode(PD, scl, ModeNumber, plotlab)
%
% Input:
%     PD         =   Matlab structure with (at least) the following fields
%        N          =   Number of nodes in mesh (numbered 1:N)
%        NodePos    =   Nx3 matrix of nodal positions
%        NE         =   Number of elements in mesh (numbered 1:NE)
%        ElmConnect =   NEx2 matrix of node numbers (for each element)
%        NM         =   Number of Materials Sets (MatsSets)
%        MatsSets   =   NMx1 matrix of Matlab structures containing (at least)
%                       the following fields
%           E          =   Young's Modulus
%           A          =   X-sec. area
%           rho        =   Mass density per unit length
%        ElmMats    =   NEx1 matrix of Material set (for each element)
%        BCType     =   Nx3 matrix of 0's (disp-type) and 1's (force-type)
%        Nmodes     =   Number of modes in the Modes field
%        Freq       =   (Nmodes)x1 vector of vibration frequencies
%        Modes      =   Nx3xNmodes array of vibration mode displacements
%
%     scl        =   number by which to scale the displacements
%     ModeNumber =   number of mode (in 1:Nmodes) to be ploted
%     plotlab    =   'y' or 'n',  If 'y', plot node and element labels
%
% Output:
%     none
%

% By: Ryan S. Elliott -- Sept. 2009; Feb. 2015

figh = figure;
%axis equal;

Matclrmp = lines;
Stressclrmp = cool;

NumNodes = PD.N;
NumElms  = PD.NE;

for i = 1:NumElms
  x = [PD.NodePos(PD.ElmConnect(i,1),1); PD.NodePos(PD.ElmConnect(i,2),1)];
  y = [PD.NodePos(PD.ElmConnect(i,1),2); PD.NodePos(PD.ElmConnect(i,2),2)];
  z = [PD.NodePos(PD.ElmConnect(i,1),3); PD.NodePos(PD.ElmConnect(i,2),3)];
  ux = scl*[PD.Modes(PD.ElmConnect(i,1),1,ModeNumber);
            PD.Modes(PD.ElmConnect(i,2),1,ModeNumber)];
  uy = scl*[PD.Modes(PD.ElmConnect(i,1),2,ModeNumber);
            PD.Modes(PD.ElmConnect(i,2),2,ModeNumber)];
  uz = scl*[PD.Modes(PD.ElmConnect(i,1),3,ModeNumber);
            PD.Modes(PD.ElmConnect(i,2),3,ModeNumber)];

  wdth = 1;
  clr = [0,0,0];

  ln(i) = line(x,y,z, ...
               'LineWidth', wdth, ...
               'Color', clr ...
              );
  defln(i) = line(x+ux,y+uy,z+uz, ...
                  'LineWidth', 3, ...
                  'Color', clr ...
                 );
end;

xRange = [min(PD.NodePos(:,1)+scl*PD.Modes(:,1,ModeNumber)), ...
          max(PD.NodePos(:,1)+scl*PD.Modes(:,1,ModeNumber))];
yRange = [min(PD.NodePos(:,2)+scl*PD.Modes(:,2,ModeNumber)), ...
          max(PD.NodePos(:,2)+scl*PD.Modes(:,2,ModeNumber))];
zRange = [min(PD.NodePos(:,3)+scl*PD.Modes(:,3,ModeNumber)), ...
          max(PD.NodePos(:,3)+scl*PD.Modes(:,3,ModeNumber))];

xticlen = 0.05*(xRange(2)-xRange(1));
yticlen = 0.05*(yRange(2)-yRange(1));
zticlen = 0.05*(zRange(2)-zRange(1));

for i = 1:NumElms
  x = [PD.NodePos(PD.ElmConnect(i,1),1); PD.NodePos(PD.ElmConnect(i,2),1)];
  y = [PD.NodePos(PD.ElmConnect(i,1),2); PD.NodePos(PD.ElmConnect(i,2),2)];
  z = [PD.NodePos(PD.ElmConnect(i,1),3); PD.NodePos(PD.ElmConnect(i,2),3)];

  if (plotlab == 'y')
    lnnum(i) = text(sum(x)/2 + .25*xticlen, sum(y)/2 + .25*yticlen, sum(z)/2 + .25*zticlen, ...
                 num2str(i), ...
                 'HorizontalAlignment', 'left', ...
                 'FontWeight', 'bold' ...
                );
  end;
end;

fc = zeros(NumNodes,3);

for i = 1:NumNodes
  x = PD.NodePos(i,1);
  y = PD.NodePos(i,2);
  z = PD.NodePos(i,3);

  if (PD.BCType(i,1) == 0)
    clr = 'm';
  else
    clr = 'c';
  end;
  if (~((PD.BCType(i,1) == 0) && (PD.BCVal(i,1) == 0.0)))
   fc(i,1) = line([x-xticlen,x+xticlen],[y,y],[z,z],'Color',clr,'LineWidth',2);
  end;

  if (PD.BCType(i,2) == 0)
    clr = 'm';
  else
    clr = 'c';
  end;
  if (~((PD.BCType(i,2) == 0) && (PD.BCVal(i,2) == 0.0)))
   fc(i,2) = line([x,x],[y-yticlen,y+yticlen],[z,z],'Color',clr,'LineWidth',2);
  end;

  if (PD.BCType(i,3) == 0)
    clr = 'm';
  else
    clr = 'c';
  end;
  if (~((PD.BCType(i,3) == 0) && (PD.BCVal(i,3) == 0.0)))
   fc(i,3) = line([x,x],[y,y],[z-zticlen,z+zticlen],'Color',clr,'LineWidth',2);
  end;

  nd(i) = line(x,y,z, ...
               'Marker', '.', ...
               'MarkerSize', 20 ...
              );

 if (plotlab == 'y')
   ndnum(i) = text(x+.5*xticlen,y+.5*yticlen,z+.5*zticlen, ...
                 num2str(i), ...
                 'HorizontalAlignment', 'left', ...
                 'Color', 'blue', ...
                 'FontWeight', 'bold', ...
                 'FontSize', 15 ...
                    );
  end;
end;

Range = [min([xRange(1), yRange(1), zRange(1)]),...
         max([xRange(2), yRange(2), zRange(2)])];
xlabel('X');
xlim(Range+(Range(2)-Range(1))*[-0.1,0.1]);
ylabel('Y');
ylim(Range+(Range(2)-Range(1))*[-0.1,0.1]);
zlabel('Z');
zlim(Range+(Range(2)-Range(1))*[-0.1,0.1]);
