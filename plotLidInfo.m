function plotLidInfo(Geometry,Grid,Lid,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotLidInfo.m
% Plot lid data
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% June 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure(indFigure); clf;
orient landscape;
colormap(jet);

% % lid depth plot
subplot(2,2,1); hold on;
[C,h]=contourf(Grid.x,Grid.y,Lid.Depth,[0:10:100]);
clabel(C,h,[0:10:100]);
axis equal; 
axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]); 
box on;
title('Permeability Barrier Depth [km]');
set(gca,'Clim',[0,100]); colorbar;

% % lid temperature plot
subplot(2,2,2); hold on;
[C,h]=contourf(Grid.x,Grid.y,Lid.T,[1200:10:1400]);
clabel(C,h,[1200:50:1400]);
axis equal; 
axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]); 
box on;
title('Permeability Barrier Temperature [degC]');
set(gca,'Clim',[1200,1400]); colorbar;

% % excess temperature plot
subplot(2,2,3); hold on;
[C,h]=contour(Grid.x,Grid.y,Lid.ExcessT,[0:10:250],'k');
clabel(C,h,[00:50:250]);
[C,h]=contourf(Grid.x,Grid.y,Lid.T-Lid.T_Solidus,[0:10:250]);
clabel(C,h,[0:50:250]);
axis equal; 
axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]); 
box on;
title('ExcessT and T-T\_Solidus at Permeability Barrier [degC]');
set(gca,'Clim',[0,250]); colorbar;

% % slope plot
subplot(2,2,4); hold on;
[C,h]=contourf(Grid.x,Grid.y,log10(Lid.Slope),[-3:0.5:1]);
clabel(C,h,[-3:1:1]);
axis equal; 
axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]); 
box on;
title('log10(Slope)');
set(gca,'Clim',[-3,1]); 
colorbar;