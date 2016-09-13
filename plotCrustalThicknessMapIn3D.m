function plotCrustalThicknessMapIn3D(SeafloorCrust,PlateBoundary_x,PlateBoundary_y,ModelBoundary_x,ModelBoundary_y,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotCrustalThicknessMapIn3D.m
% Plot 3D crustal thickness map 
% Hailong Bai
% October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

figure(indFigure); clf; hold on;

surface(SeafloorCrust.TileCenter_x,SeafloorCrust.TileCenter_y,SeafloorCrust.Thickness); shading flat;
RectangleHeight=[0,15];
for iRectangle=1:numel(PlateBoundary_x)-1;
    PlateBoundaryRectangle(iRectangle)=fill3(PlateBoundary_x(iRectangle+[0,0,1,1,0]),PlateBoundary_y(iRectangle+[0,0,1,1,0]),...
        RectangleHeight([1,2,2,1,1]),'m');   
end
alpha(PlateBoundaryRectangle,0.5);
box on;

set(gca,'DataAspectRatio',[1,1,0.5]);
set(gca,'xlim',ModelBoundary_x);
set(gca,'ylim',ModelBoundary_y);
set(gca,'zlim',RectangleHeight);
set(gca,'ztick',[]);

% Position=get(gca,'position');
% set(gca,'position',[Position(1)-0.07,Position(2)+0.11,Position(3),Position(4)]);
set(gca,'projection','perspective');
% camva(30);
CameraTarget=camtarget;
InitialCameraPosition=campos;
CameraPosition=CameraTarget+(InitialCameraPosition-CameraTarget)/5;
set(gca,'cameraPosition',CameraPosition);

Light1=lightangle(150,20);
Light2=lightangle(-60,50);
view(75,25);