function Shoulder=meltSwath(Model,Geometry,Grid,Res,Lid,Shoulder,LidDepthLimit,Switch_UseCOMSOL,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltSwath.m
% Connect adjacent melt trajectories into swaths
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Model                   : 3D model results file
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |......
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates [km]
%   Res
%       |.nMeltSwath        : Melt swath sampling size
%       |......
%   Lid
%       |.Depth             : Depth of lid [km]
%       |.T                 : Temperature of lid [degC]
%       |.Slope             : Slope of lid
%       |......
%   Shoulder(ind)
%       |.Trajectory(ind)
%           |.x, .y         : Coordinates of melt trajectory [km]
%           |.Distance      : Distance along melt trajectory [km]
%           |.Depth         : Depth of melt trajectory [km]
%           |.DepthGrad     : Slope of melt trajectory in x or y direction
%           |.Slope         : Slope along melt trajectory
%           |.T             : Temperature along melt trajectory [degC]
%           |.T_Solidus     : Solidus temperature along melt trajectory [degC]
%           |.LineStart     : Seed coordinate [km]
%   LidDepthLimit           : Lower limit for permeability barrier [km]
%   Switch_UseCOMSOL        : Switch for COMSOL model usage, 1: using COMSOL, 2: not using COMSOL
%   indFigure               : Figure index
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Shoulder(ind)
%       |.Swath(ind)
%           |.Polygon(ind)
%               |.x,.y      : Coordinates of polygon outline [km]
%           |.PolygonArea   : Area of each polygon [km^2]
%           |.SwathArea     : Area of swath [km^2]
%           |.PolygonCenter
%               |_x, _y, _z : Coordinates of polygon centers [km]
%               |_T         : Temperature at polygon centers [degC]
%           |.PolygonSlope  : Slope at polygon centers
%           |.Column
%               |_x, _y, _z : Coordinate at vertical sampling column below polygon [km]
%               |_T         : Temperature at sampling column [degC]
%               |_VerticalVelocity : Vertical velocity at sampling column [m/s]
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   nShoulder               : Number of melt shoulders
%   nLine                   : Number of melt trajectories in each shoulder
%   nPolygon                : Number of polygons in each melt swath
%   iShoulder               : Shoulder index
%   iSwath                  : Swath index
%   iTrajectory             : Trajectory index
%   iPolygon                : Polygon index
%   interpolant_vz          : Interpolant for vertical velocity
%   interpolant_T           : Interpolant for temperature
%   ExtendToSegmentDistance : When line ends get in range, extend line to plate boundary segment [km] 
%   Swath                   : Swath structure for swath information storage
%   CurrentLine             : Current line information
%   PreviousLine            : Previous line information
%   UpdatedLine             : Resampled line information
%   CurrentLineDis          : Distance along current line [km]
%   PreviousLineDis         : Distance along previous line [km]
%   UpdatedLineDis          : Distance along resampled line [km]
%   CurrentOnSegmentInfo    : Information of current line on-segment point
%   PreviousOnSegmentInfo   : Information of previous line on-segment point
%   PolygonArea             : Area of polygon [km^2]
%   PolygonCenter           : Coordinates and temperature information at polygon center
%   PolygonSlope            : Slope of polygon
%   Column_z                : Sampling column depth
%   Column_T                : Sampling column temperature
%   Column_VerticalVelocity : Sampling column vertical velocity
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   lineSampleDepth
%   clockwiseAreaCalc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

ExtendToSegmentDistance=30;
nShoulder=numel(Shoulder);
if ~Switch_UseCOMSOL
    interpolant_vz=scatteredInterpolant(Model{1},Model{2},Model{3},Model{6});
    interpolant_T=scatteredInterpolant(Model{1},Model{2},Model{3},Model{7});
end
for iShoulder=1:nShoulder; % shoulder index
    disp(sprintf('>>>    Defining swaths for shoulder %d',iShoulder))
    iSwath=0; % swath index
    iTrajectory=1; % trajectory index
    Swath=[];
    [CurrentLine,CurrentOnSegmentInfo]=lineSampleDepth(Geometry,Grid,Res,Lid,Shoulder(iShoulder).Trajectory(iTrajectory),...
        ExtendToSegmentDistance,LidDepthLimit);
    CurrentLineDis=CurrentLine(:,4);
    
    % % to form swath and polygons, number of samples on trajectories should be same
    nLine=numel(Shoulder(iShoulder).Trajectory);
    for iTrajectory=2:nLine;
        iSwath=iSwath+1;
        % store previous line info
        PreviousLineDis=CurrentLineDis;
        PreviousLine=CurrentLine;
        PreviousOnSegmentInfo=CurrentOnSegmentInfo;
        if ~isempty(Shoulder(iShoulder).Trajectory(iTrajectory).Distance); % If empty, use old Trajectory...
            [CurrentLine,CurrentOnSegmentInfo]=lineSampleDepth(Geometry,Grid,Res,Lid,Shoulder(iShoulder).Trajectory(iTrajectory),...
                ExtendToSegmentDistance,LidDepthLimit);
        end
        CurrentLineDis=CurrentLine(:,4);
        if numel(CurrentLineDis)<numel(PreviousLineDis); % previous trajectory has more sample points, resample current trajectory
            UpdatedLineDis=linspace(CurrentLineDis(1),CurrentLineDis(end),numel(PreviousLineDis));
            UpdatedLine=interp1(CurrentLineDis,CurrentLine,UpdatedLineDis);
            nPolygon=numel(UpdatedLineDis)-1; % Number of polygons
            for iPolygon=1:nPolygon; % form melt polygons between previous and current trajectory
                Swath(iSwath).Polygon(iPolygon).x=[PreviousLine(iPolygon+[0,1],1);UpdatedLine(iPolygon+[1,0],1)];
                Swath(iSwath).Polygon(iPolygon).y=[PreviousLine(iPolygon+[0,1],2);UpdatedLine(iPolygon+[1,0],2)];
            end
        else
            if numel(CurrentLineDis)>numel(PreviousLineDis); % current trajectory has more points, resample previous trajectory 
                UpdatedLineDis=linspace(PreviousLineDis(1),PreviousLineDis(end),numel(CurrentLineDis));
                PreviousLine=interp1(PreviousLineDis,PreviousLine,UpdatedLineDis);
                PreviousLineDis=UpdatedLineDis;
            end
            nPolygon=numel(CurrentLineDis)-1;
            for iPolygon=1:nPolygon; % form melt polygons between previous and current trajectory
                Swath(iSwath).Polygon(iPolygon).x=[PreviousLine(iPolygon+[0,1],1);CurrentLine(iPolygon+[1,0],1)];
                Swath(iSwath).Polygon(iPolygon).y=[PreviousLine(iPolygon+[0,1],2);CurrentLine(iPolygon+[1,0],2)];
            end
        end
        if PreviousOnSegmentInfo(3)~=CurrentOnSegmentInfo(3) % previous and current trajectories end up in different plate boundary segment
            % % add one more polygon to include plate boundary corner
            if CurrentOnSegmentInfo(3)>PreviousOnSegmentInfo(3);
                Swath(iSwath).Polygon(nPolygon+1).x=[PreviousLine(nPolygon+1,1);...
                    Geometry.PlateBoundary.x([PreviousOnSegmentInfo(3)+1:CurrentOnSegmentInfo(3)])';CurrentLine(nPolygon+1,1)];
                Swath(iSwath).Polygon(nPolygon+1).y=[PreviousLine(nPolygon+1,2);...
                    Geometry.PlateBoundary.y([PreviousOnSegmentInfo(3)+1:CurrentOnSegmentInfo(3)])';CurrentLine(nPolygon+1,2)];
            else
                Swath(iSwath).Polygon(nPolygon+1).x=[PreviousLine(nPolygon+1,1);...
                    Geometry.PlateBoundary.x([PreviousOnSegmentInfo(3):-1:CurrentOnSegmentInfo(3)+1])';CurrentLine(nPolygon+1,1)];
                Swath(iSwath).Polygon(nPolygon+1).y=[PreviousLine(nPolygon+1,2);...
                    Geometry.PlateBoundary.y([PreviousOnSegmentInfo(3):-1:CurrentOnSegmentInfo(3)+1])';CurrentLine(nPolygon+1,2)];
            end
            nPolygon=nPolygon+1;
        end
        
        % % get additional information about polygon
        PolygonArea=[]; Column_z=[]; PolygonCenter_x=[]; PolygonCenter_y=[]; PolygonCenter_z=[]; PolygonCenter_T=[];
        for iPolygon=1:nPolygon;
            PolygonArea(iPolygon)=clockwiseAreaCalc(Swath(iSwath).Polygon(iPolygon).x,Swath(iSwath).Polygon(iPolygon).y);
            PolygonCenter_x(iPolygon)=mean(Swath(iSwath).Polygon(iPolygon).x);
            PolygonCenter_y(iPolygon)=mean(Swath(iSwath).Polygon(iPolygon).y);
            PolygonCenter_z(iPolygon)=interp2(Grid.x,Grid.y,Lid.Depth,PolygonCenter_x(iPolygon),PolygonCenter_y(iPolygon));
            PolygonCenter_T(iPolygon)=interp2(Grid.x,Grid.y,Lid.T,PolygonCenter_x(iPolygon),PolygonCenter_y(iPolygon));
            PolygonSlope(iPolygon)=interp2(Grid.x,Grid.y,Lid.Slope,PolygonCenter_x(iPolygon),PolygonCenter_y(iPolygon));
            Column_z(iPolygon,[1:Res.nMeltSwath])=linspace(PolygonCenter_z(iPolygon),max(Geometry.ModelBoundary.z),Res.nMeltSwath);
        end
        Swath(iSwath).PolygonArea=PolygonArea'; % area of polygon
        Swath(iSwath).SwathArea=sum(PolygonArea); % area of swath
        Swath(iSwath).PolygonCenter_x=PolygonCenter_x'; % x coordinate of polygon center
        Swath(iSwath).PolygonCenter_y=PolygonCenter_y'; % y coordinate of polygon center
        Swath(iSwath).PolygonCenter_z=PolygonCenter_z'; % z coordinate of polygon center
        Swath(iSwath).PolygonCenter_T=PolygonCenter_T'; % temperature at polygon center
        Swath(iSwath).PolygonSlope=PolygonSlope'; % slope at polygon center
        Swath(iSwath).Column_x=repmat(PolygonCenter_x',[1,Res.nMeltSwath]); % x coordinate of vertical sampling column below polygon
        Swath(iSwath).Column_y=repmat(PolygonCenter_y',[1,Res.nMeltSwath]); % y coordinate of vertical sampling column below polygon
        Swath(iSwath).Column_z=Column_z; % z coordinate of vertical sampling column below polygon
        if Switch_UseCOMSOL
            [Column_T,Column_VerticalVelocity]=mphinterp(Model,{'T','w'},'coord',...
                [Swath(iSwath).Column_x(:),Swath(iSwath).Column_y(:),-Swath(iSwath).Column_z(:)]','unit',{'degC','m/s'},'dataset','dset1','solnum','end'); 
        else
            Column_T=interpolant_T(Swath(iSwath).Column_x(:),Swath(iSwath).Column_y(:),-Swath(iSwath).Column_z(:))-273.15;
            Column_VerticalVelocity=interpolant_vz(Swath(iSwath).Column_x(:),Swath(iSwath).Column_y(:),-Swath(iSwath).Column_z(:));
        end
        Swath(iSwath).Column_T=reshape(Column_T,size(Column_z)); % temperature of sampling column
        Swath(iSwath).Column_VerticalVelocity=reshape(Column_VerticalVelocity,size(Column_z)); % vertical velocity of sampling column
    end
    Shoulder(iShoulder).Swath=Swath;
end

%% Melt Swath Plot

figure(indFigure); clf; hold on;
contourf(Grid.x,Grid.y,-Lid.Depth,[0:-2:-100],'edgecolor','none'); colorbar;
set(gca,'Clim',[-100,0]);
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
for iShoulder=1:numel(Shoulder);
    for iSwath=1:numel(Shoulder(iShoulder).Swath);
        for iPolygon=1:numel(Shoulder(iShoulder).Swath(iSwath).Polygon);
              patch(Shoulder(iShoulder).Swath(iSwath).Polygon(iPolygon).x,Shoulder(iShoulder).Swath(iSwath).Polygon(iPolygon).y,...
                'y','facecolor','none','edgecolor','white','linewidth',0.7);
        end
    end
end
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Melt Swath and Polygons');

return