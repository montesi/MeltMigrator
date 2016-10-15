function Shoulder=meltTrajectory(Geometry,Grid,Res,Lid,Saddle,Seed,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltTrajectory.m
% Compute melt trajectories from seeds 
% Lines are grouped into shoulders separated by saddles
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |......
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates [km]
%   Res
%       |.MeltTrajectory    : Distance between adjacent points on melt trajectory [km]
%       |......
%   Lid
%       |.Depth             : Depth of lid [km]
%       |......
%   Saddle(ind)
%       |.Point             : Position of saddle point [km]
%       |.Down              : Downgoing saddle line coordinates [km]
%       |.DisDown           : Distance along downgoing saddle line [km]
%       |.Up                : Upgoing saddle line coordinates [km]
%       |.DisUp             : Distance along upgoing saddle line [km]
%   Seed(ind)
%       |.x, .y             : Coordinates of melt line seeds [km]
%   indFigure               : Figure index
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
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
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   iSaddle                 : Saddle index
%   iShoulder               : Shoulder index
%   iSeed                   : Seed index
%   iLine                   : Line index
%   iTrajectory             : Trajectory index
%   iShoulderBoundary       : Shoulder boundary index
%   iModifier               : Index modifier for adjacent saddle line selection
%   iDepthChange            : Index for points on melt line with depth change larger than criteria
%   iMinDepth               : Index for minimum depth along melt line [km]
%   iSelect                 : Index for selected points on melt line 
%   iTracker                : Tracker for progress in waitbar
%   nTracker                : Tracker for waitbar display
%   ShoulderBoundaryUp      : Coordinates of upgoing shoulder boundary [km]
%   ShoulderBoundaryDown    : Coordinates of downgoing shoulder boundary [km]
%   ShoulderBoundaryDisUp   : Distance along upgoing shoulder boundary [km]
%   ShoulderBoundaryDisDown : Distance along downgoing shoulder boundary [km]
%   ShoulderBoundary(ind)
%       |.Position          : Coordinates of shoulder boundary [km]
%       |.Distance          : Distance along shoulder boundary [km]
%       |.Segment           : Shoulder boundary segment built for intersection determination [km]
%   DistanceUp              : Upgoing trajectory search range [km]
%   DistanceDown            : Downgoing trajectory search range [km]
%   options                 : Options for ode solver used in trajectory search
%   LineUp                  : Coordinates of upgoing line [km]
%   LineDown                : Coordinates of downgoing line [km]
%   LineDisUp               : Distance along upgoing line [km]
%   LineDisDown             : Distance along downgoing line [km]
%   Line                    : Coordinates of melt line [km]
%   LineDis                 : Distance along melt line [km]
%   LineDepth               : Depth of melt line [km]
%   SeedSegment             : Seed segment built for intersection determination [km]
%   indIntersection         : Line intersection indicator
%   Intersection            : Line intersection structure
%   MinimumDepthChange      : Minimum depth change considered, tweek to clean up line tail [km]
%   MinDepth                : Minimum depth along melt line [km]
%   Trajectory              : Coordinates of melt trajectory [km]
%   WaitBar                 : Waitbar header
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   lineTruncate
%   lineSegmentIntersect
%   alongSlopeGrad
%   lineCleanUp
%   lineStore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Shoulder Generation

Shoulder=[]; % initialize trajectories
iShoulder=0; % shoulder index
DistanceUp=200; % uphill distance for melt line search
DistanceDown=300; % downhill distance for melt line search
options=odeset('reltol',1e-8); % optimset('display','off')
MinimumDepthChange=5e-1; % minimum depth change considered, tweek to clean up line tail noise

% % combine adjacent saddle lines to form boundaries for shoulders
for iSaddle=1:numel(Saddle);
    % % start with a downgoing saddle line
	[ShoulderBoundaryDisDown,ShoulderBoundaryDown]=lineTruncate(Saddle(iSaddle).DisDown,Saddle(iSaddle).Down);
	for iModifier=0:1; % modifier helps to select the next adjacent upgoing saddle line (clockwise)
		if rem(iSaddle,2)==1;
			[ShoulderBoundaryDisUp,ShoulderBoundaryUp]=lineTruncate(Saddle(iSaddle+iModifier).DisUp,Saddle(iSaddle+iModifier).Up);
		else
			[ShoulderBoundaryDisUp,ShoulderBoundaryUp]=lineTruncate(Saddle(iSaddle-iModifier).DisUp,Saddle(iSaddle-iModifier).Up);
		end
		iShoulder=iShoulder+1;
        % % each shoulder boundary is made of one upgoing saddle line and one downgoing saddle line
        ShoulderBoundary(iShoulder).Position=[flipud(ShoulderBoundaryDown);Saddle(iSaddle).Point(1:2);ShoulderBoundaryUp([2:end],:)];
        ShoulderBoundary(iShoulder).Distance=[flipud(-ShoulderBoundaryDisDown);0;ShoulderBoundaryDisUp([2:end])];
        % % generate seudo boundary segment for intersection determination
        ShoulderBoundary(iShoulder).Segment=[ShoulderBoundary(iShoulder).Position(1:end-1,:),ShoulderBoundary(iShoulder).Position(2:end,:)]; 
    end
end

%% Melt Trajectory Generation

iShoulder=0;
iTracker=0;
nTracker=0;
for iSeed=1:numel(Seed);
    nTracker=nTracker+numel(Seed(iSeed).x);
end
WaitBar=waitbar(0,'Pleas wait...');

for iSeed=1:numel(Seed);
    iShoulder=iShoulder+1; % shoulder index
    % % difference between 'Line' and 'Trajectory' here is numbering:
    % % 'Line' comes directly from seeds, its numbering goes with seed point number; 
    % % 'Trajectory' is numbered within one shoulder, and is set to 0 every time shoulder number changes
    iLine=0; % line index
    iTrajectory=0; % trajectory index

    SeedSegment=[Seed(iSeed).x(1:end-1,:),Seed(iSeed).y(1:end-1,:),...
        Seed(iSeed).x(2:end,:),Seed(iSeed).y(2:end,:)]; 
    clear Intersection;
    indIntersection=zeros(numel(ShoulderBoundary),numel(Seed(iSeed).x)-1);
    for iShoulderBoundary=1:numel(ShoulderBoundary);
        Intersection(iShoulderBoundary)=lineSegmentIntersect(ShoulderBoundary(iShoulderBoundary).Segment,SeedSegment);
        indIntersection(iShoulderBoundary,:)=max(Intersection(iShoulderBoundary).intAdjacencyMatrix);
    end
    
    while iLine<numel(Seed(iSeed).x);
        iTracker=iTracker+1;
        iLine=iLine+1;
%         disp(sprintf('>>>    Working on seed group %d, line %d/%d',iSeed,iLine,numel(Seed(iSeed).x)));
        [LineDisUp,LineUp]=ode23(@(l,P) alongSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,1),...
            [0:Res.MeltTrajectory:DistanceUp],[Seed(iSeed).x(iLine),Seed(iSeed).y(iLine)],options);
        [LineDisDown,LineDown]=ode23(@(l,P) alongSlopeGrad(P,Grid,Geometry.ModelBoundary,Lid,2),...
            [0:Res.MeltTrajectory:DistanceDown],[Seed(iSeed).x(iLine),Seed(iSeed).y(iLine)],options);
        
        % % treatment for trajectory data
        % % remove overlapping points
        [LineDisUp,LineUp]=lineTruncate(LineDisUp,LineUp);
        [LineDisDown,LineDown]=lineTruncate(LineDisDown,LineDown);
        % % combine downgoing and upgoing parts to assemble a full trajectory
        Line=[flipud(LineDown);LineUp([2:end],:)];
        LineDis=[flipud(-LineDisDown);LineDisUp([2:end])];
        % % truncate trajectory to its minimum depth
        LineDepth=interp2(Grid.x,Grid.y,Lid.Depth,Line(:,1),Line(:,2));
        iDepthChange=find(abs(diff(LineDepth))>MinimumDepthChange);
        [MinDepth,iMinDepth]=min(LineDepth(iDepthChange+1));
        iSelect=[iDepthChange(1:iMinDepth);iDepthChange(iMinDepth)+1]; 
        iSelect=iSelect(find(~isnan(LineDepth(iSelect))));
        if(isempty(iSelect));
            iSelect=1;
        end
        Line=Line(iSelect,:);
        LineDis=LineDis(iSelect,:);
        
        Trajectory=lineStore(Line,LineDis,Grid,Lid,[Seed(iSeed).x(iLine),Seed(iSeed).y(iLine)]);
        iTrajectory=iTrajectory+1;
        Shoulder(iShoulder).Trajectory(iTrajectory)=Trajectory;
        waitbar(iTracker/nTracker,WaitBar,sprintf('Store line %g-%g/%g as shoulder %g, line %g',...
            iSeed,iLine,numel(Seed(iSeed).x),iShoulder,iTrajectory));

        if iLine<numel(Seed(iSeed).x);
            if max(indIntersection(:,iLine)); % if seed group intersects with shoulder boundary, move on to new shoulder
                iShoulder=iShoulder+1;
                iTrajectory=0;
            end
        end
    end
end
delete(WaitBar);

%% Melt Trajectory Plot

figure(indFigure); clf; hold on;
contourf(Grid.x,Grid.y,-Lid.Depth,[0:-2:-100],'edgecolor','none'); colorbar;
set(gca,'Clim',[-100,0]);
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
for iShoulder=1:numel(Shoulder);
    for iTrajectory=1:numel(Shoulder(iShoulder).Trajectory)
        plot(Shoulder(iShoulder).Trajectory(iTrajectory).x,Shoulder(iShoulder).Trajectory(iTrajectory).y,'-w','linewidth',1);
    end
end
for iSeed=1:numel(Seed);
    plot(Seed(iSeed).x,Seed(iSeed).y,'m.');
end
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Melt Trajectories');

return
