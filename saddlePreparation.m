function Saddle=saddlePreparation(Geometry,Grid,Res,Lid,SaddleInitialPoint,SaddleSearchRadius,SaddleFlatCutoff,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saddlePreparation.m
% Determine location of saddle lines in lid
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% June 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |......
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates [km]
%   Res
%       |.Saddle            : Distance between points along saddle line [km]
%   Lid
%       |.Depth             : Depth of lid [km]
%       |.DepthGrad         : Slope of lid in x or y direction
%       |.SlopeGrad         : Slope gradient in x or y direction
%       |......
%   SaddleInitialPoint      : Starting point for saddle point search [km]
%   SaddleSearchRadius      : Radius for saddle line search [km]
%   SaddleFlatCutoff        : Criterion for saddle flatness
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Saddle(ind)
%       |.Point             : Position of saddle point [km]
%       |.Down              : Downgoing saddle line coordinates [km]
%       |.DisDown           : Distance along downgoing saddle line [km]
%       |.Up                : Upgoing saddle line coordinates [km]
%       |.DisUp             : Distance along upgoing saddle line [km]
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   indSaddle               : Index of saddle
%   SaddlePointSearchLine   : Saddle point search path coordinates [km]
%   SaddlePointSearchDis    : Distance along saddle point search path [km]
%   SaddlePoint_x, y, z     : Saddle point coordinates [km]
%   SaddleLineSearchDirection : Direction of saddle line search [rad]
%   SaddleLineSearchCircle
%       |_x, _y, _z         : Saddle line search circle coordinates [km]
%   SearchDepth             : Depth of starting point of saddle line [km]
%   indSearch               : Index of chosen saddle line seeds
%   SaddleLine
%       |_Down1, _Down2     : Downgoing saddle line coordinates [km]
%       |_Up1, _Up2         : Upgoing saddle line coordinates [km]
%   SaddleLineDis
%       |_Down1, _Down2     : Distance along downgoing saddle lines [km]
%       |_Up1, _Up2         : Distance along upgoing saddle lines [km]
%   NewSaddleLine
%       |_Up1, _Up2         : Sorted upgoing saddle lines coordinates [km]
%   NewSaddleLineDis
%       |_Up1, _Up2         : Distance along sorted upgoing saddle lines [km]
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   alongSlopeSlopeGrad
%   alongSlopeGrad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

% % plot map of lid depth and slope
figure(indFigure); clf; hold on;
colormap(summer);
contourf(Grid.x,Grid.y,Lid.Depth,'edgecolor','none'); colorbar;
contour(Grid.x,Grid.y,log10(Lid.Slope),[-3:0.1:0.5]);
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Permeability Barrier Depth and Saddle Lines');

% % calculate saddle line
indSaddle=0; % initialize saddle index
for i=1:size(SaddleInitialPoint,1); % work on each initial saddle point
    disp(sprintf('>>> Working on saddle %d/%d',i,size(SaddleInitialPoint,1)));
    
    % % find saddle point from the initial point
    [SaddlePointSearchDis,SaddlePointSearchLine]=ode45(@(l,P) alongSlopeSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,SaddleFlatCutoff,1),... 
        [0,20],[SaddleInitialPoint(i,1),SaddleInitialPoint(i,2)],optimset('display','off'));
    plot(SaddlePointSearchLine(:,1),SaddlePointSearchLine(:,2),'k');
    plot(SaddlePointSearchLine(end,1),SaddlePointSearchLine(end,2),'om');

    % % saddle point coordinates
    SaddlePoint_x=SaddlePointSearchLine(end,1);
    SaddlePoint_y=SaddlePointSearchLine(end,2);
    SaddlePoint_z=interp2(Grid.x,Grid.y,Lid.Depth,SaddlePoint_x,SaddlePoint_y);
    disp(sprintf('>>>    Flat area at (%g, %g, %g)',SaddlePoint_x,SaddlePoint_y,SaddlePoint_z));

    % % find saddle lines along a circle around the saddle point
    SaddleLineSearchDirection=linspace(0,pi,100);
    SaddleLineSearchCircle_x=SaddlePoint_x+cos(SaddleLineSearchDirection)*SaddleSearchRadius;
    SaddleLineSearchCircle_y=SaddlePoint_y+sin(SaddleLineSearchDirection)*SaddleSearchRadius;
    SaddleLineSearchCircle_z=interp2(Grid.x,Grid.y,Lid.Depth,SaddleLineSearchCircle_x,SaddleLineSearchCircle_y);

    % % downgoing saddle lines
    disp('>>>    Searching down');
    [SearchDepth,indSearch]=max(SaddleLineSearchCircle_z); % starting from the deepest points
    [SaddleLineDis_Down1,SaddleLine_Down1]=ode23(@(l,P) alongSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,2),...
        [0:Res.Saddle:1000],[SaddleLineSearchCircle_x(indSearch),SaddleLineSearchCircle_y(indSearch)],optimset('display','off'));
    [SaddleLineDis_Down2,SaddleLine_Down2]=ode23(@(l,P) alongSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,2),...
        [0:Res.Saddle:1000],[SaddlePoint_x,SaddlePoint_y]*2-[SaddleLineSearchCircle_x(indSearch),SaddleLineSearchCircle_y(indSearch)],optimset('display','off'));
    plot(SaddleLine_Down1(:,1),SaddleLine_Down1(:,2),'r');
    plot(SaddleLine_Down2(:,1),SaddleLine_Down2(:,2),'b');
    
    % % upgoing saddle lines
    disp('>>>    Searching up')
    [SearchDepth,indSearch]=min(SaddleLineSearchCircle_z); % starting from the shallowest points
    [SaddleLineDis_Up1,SaddleLine_Up1]=ode23(@(l,P) alongSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,1),...
        [0:Res.Saddle:1000],[SaddleLineSearchCircle_x(indSearch),SaddleLineSearchCircle_y(indSearch)],odeset('reltol',1e-8)); %odeset('reltol',1e-8)
    [SaddleLineDis_Up2,SaddleLine_Up2]=ode23(@(l,P) alongSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,1),...
        [0:Res.Saddle:1000],[SaddlePoint_x,SaddlePoint_y]*2-[SaddleLineSearchCircle_x(indSearch),SaddleLineSearchCircle_y(indSearch)],odeset('reltol',1e-8));
    plot(SaddleLine_Up1(:,1),SaddleLine_Up1(:,2),'r');    
    plot(SaddleLine_Up2(:,1),SaddleLine_Up2(:,2),'b');
    
    % % prepare to number saddle lines based on the y coordinates
    if mean(SaddleLine_Up1(:,2))>mean(SaddleLine_Up2(:,2));
        NewSaddleLine_Up1=SaddleLine_Up1;
        NewSaddleLineDis_Up1=SaddleLineDis_Up1;
        NewSaddleLine_Up2=SaddleLine_Up2;
        NewSaddleLineDis_Up2=SaddleLineDis_Up2;
    else
        NewSaddleLine_Up1=SaddleLine_Up2;
        NewSaddleLineDis_Up1=SaddleLineDis_Up2;
        NewSaddleLine_Up2=SaddleLine_Up1;
        NewSaddleLineDis_Up2=SaddleLineDis_Up1;
    end    
    
    % % store and output saddle data
    indSaddle=indSaddle+1;
    Saddle(indSaddle).Point=[SaddlePoint_x,SaddlePoint_y,SaddlePoint_z];
    Saddle(indSaddle).Down=SaddleLine_Down1;
    Saddle(indSaddle).DisDown=SaddleLineDis_Down1+Res.Saddle;
    text(SaddleLine_Down1(ceil(size(SaddleLine_Down1,1)/2),1)-5,SaddleLine_Down1(ceil(size(SaddleLine_Down1,1)/2),2)-5,sprintf('%gD',indSaddle));
    Saddle(indSaddle).Up=NewSaddleLine_Up1;
    Saddle(indSaddle).DisUp=NewSaddleLineDis_Up1+Res.Saddle;
    text(NewSaddleLine_Up1(ceil(size(NewSaddleLine_Up1,1)/2),1)-5,NewSaddleLine_Up1(ceil(size(NewSaddleLine_Up1,1)/2),2)-5,sprintf('%gU',indSaddle));

    indSaddle=indSaddle+1;
    Saddle(indSaddle).Point=[SaddlePoint_x,SaddlePoint_y,SaddlePoint_z];
    Saddle(indSaddle).Down=SaddleLine_Down2;
    Saddle(indSaddle).DisDown=SaddleLineDis_Down2+Res.Saddle;
    text(SaddleLine_Down2(ceil(size(SaddleLine_Down2,1)/2),1)+5,SaddleLine_Down2(ceil(size(SaddleLine_Down2,1)/2),2)+5,sprintf('%gD',indSaddle));
    Saddle(indSaddle).Up=NewSaddleLine_Up2;
    Saddle(indSaddle).DisUp=NewSaddleLineDis_Up2+Res.Saddle;
    text(NewSaddleLine_Up2(ceil(size(NewSaddleLine_Up2,1)/2),1)+5,NewSaddleLine_Up2(ceil(size(NewSaddleLine_Up2,1)/2),2)+5,sprintf('%gU',indSaddle));

    disp(sprintf('>>> Saddle %d/%d completed',i,size(SaddleInitialPoint,1)));
end
return