function [SeafloorHistory]=crustalHistory(Geometry,Crust,SpreadingDirection,TrimDistance,SmoothingWidth,AccretWidth_Ridge,AccretWidth_Transform,CellSize,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crustalHistory.m
% Integrate crustal accretion profile along spreading direction and generate a map of crustal thickness
% Hailong Bai
% October 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.PlateBoundaryLength : Length of plate boundary segments [km]
%       |.PlateBoundaryAngleDeg : Angle of plate boundary segments (relative to x-direction) [deg]
%       |.PlateBoundaryType : Indicator of plate boundary type. 1: ridge, 2: transform, 3: oblique segment
%       |......
%   Crust
%       |.Thickness         : Distribution of erupted crust along plate boundary [km] 
%       |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%       |......
%   SpreadingDirection      : Direction of spreading [deg] 0: x-direction, 90: y-direction 
%   TrimDistance            : Distance to trim close to boundary [km]
%   SmoothingWidth          : Crust-level melt redistribution length (horizontal dike propagation distance, or smoothing width) [km]
%   AccretWidth_Ridge       : Accretion zone width along ridge segment [km]
%   AccretWidth_Transform   : Accretion zone width along transform faults [km] 
%   CellSize                : Size of tiles for crustal thickness calculation [km] 
%   indFigure               : Figure index
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   SeafloorHistory
%       |.Thickness  : Accummulated crustal thickness in the tile [km]
%       |.TileCenter_x, _y  : Tile center coordinates [km]
%       |.AverageSegmentCrust : Averaged accummulated crustal thickness along plate boundaries [km]
%       |.AverageWholdDomain : Averaged crustal thickness in the whole domain [km]
%       |.AverageRidgeDomain : Averaged crustal thickness in the ridge domain [km]
%       |.AverageTransformDomain : Averaged crustal thickness in the transform domain [km]
%       |.AverageObliqueDomain : Averaged crustal thickness in the oblique domain [km]
%       |.AccretWidth_Ridge : Accretion zone width along ridge segment [km]
%       |.AccretWidth_Transform : Accretion zone width along transform faults [km]
%       |.SmoothingWidth    : Smoothing width [km]
%       |.CellSize          : Size of tiles for crustal thickness calculation [km] 
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   iStart                  : Index of starting point of line segment of interest 
%   iEnd                    : Index of ending point of line segment of interest 
%   iSegment                : Index of plate boundary segment
%   iRidge                  : Index of tiles affected by ridge segment
%   iTransform              : Index of tiles affected by transform faults
%   iOblique                : Index of tiles affected by oblique segment
%   iTargetSegment          : Index of plate boundary segment of interest 
%   iRow                    : Index of row
%   iCol                    : Index of column
%   iTracker                : Tracker for progress in waitbar
%   nTracker                : Tracker for waitbar display
%   Tile
%       |.x, .y             : Coordinates of tile corners [km]
%       |.Center_x, _y      : Coordinate of tile center [km]
%       |.Area              : Area of tile [km^2]
%       |.CrustalThickness  : Thickness of crust in tile [km]
%       |.OriginPoint       : Coordinates of point on plate boundary affecting tile [km] 
%   EffectiveSmoothingWidth : Effective smoothing width (smoothing width cannot be larger than shortest plate boundary segment) [km] 
%   CrustalThickness        : Crustal thickness profile along plate boundaries [km] 
%   FiltedCrustalThickness  : Filted crustal thickness profile along plate boundaries [km] 
%   minDis                  : Distance between starting and ending point of a plate boundary segment [km] 
%   CompensateLength        : Distance added on the side of plate boundaries [km] 
%   PlateBoundary_x, _y     : Reoriented plate boundary coordinates [km]
%   ModelBoundary_x, _y     : Reoriented model boundary coordinates [km]
%   SamplingSize_x, _y      : Sampling size of reoriented geometry 
%   Vector_x, _y            : Sampling vector in x- and y-direction [km]
%   Grid_x                  : Sampling grid of reoriented geometry [km]
%   TilePolygon             : Tile polygon for distance determination 
%   SegmentPolygon          : Segment polygon for distance determination 
%   Distance                : Distance between tile and corresponding segment [km] 
%   AccretWidth             : Width of accretion zone aroung plate boundary segment of interest [km] 
%   AccretZoneBoundary      : Accretion zone boundary coordinates [km]
%   ObliqueAccretBoundary   : Oblique accretion zone boundary coordinates [km] 
%   Matrix                  : Matrix for calculation of length of tile projection on plate boundary 
%   Projection              : Length of tile projection on plate boundary [km] 
%   SegmentNumber           : Index of plate boundary segment tile corners associated with
%   AlongSegmentCoordinate  : Along plate boundary coordinates of tile corner projection [km]
%   DistanceToSegment       : Distance between tile corners and plate boundary [km]
%   Assign                  : Along plate boundary coordinates of tile [km] 
%   LeftPlate               : Accretion index of left-hand plate
%   RightPlate              : Accretion index of right-hand plate
%   RidgeDomainIndex        : Index of tiles affected by ridge segment
%   TransformDomainIndex    : Index of tiles affected by transform fault
%   ObliqueDomainIndex      : Index of tiles affected by oblique segment
%   AverageWidth            : Width of averaging along plate boundaries [km] 
%   AverageSegmentCrust     : Averaged crustal thickness along plate boundaries [km] 
%   SegmentSamplingSize     : Size of along plate boundary sampling
%   SegmentSampling         : Plate boundary sampling vector [km]
%   CurrentAverageSegmentCrust : Averaged crustal thickness along current plate boundary segment [km] 
%   WaitBar                 : Waitbar header
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   assignSegment
%   fastSmooth
%   clockwiseAreaCalc
%   minDistBetweenTtwoPolygons
%   plotCrustalThicknessMapIn3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Crustal Accretion Profile Along Plate Boundaries

% % replace original crustal thickness with filtered one
EffectiveSmoothingWidth=max(max(find((Crust.AlongPlateBoundaryCoordinate-Crust.AlongPlateBoundaryCoordinate(1))<SmoothingWidth)),1);
CrustalThickness=[];
iStart=1;
for iSegment=1:numel(Geometry.PlateBoundaryLength);
    [minDis,iEnd]=min(abs(Crust.AlongPlateBoundaryCoordinate-sum(Geometry.PlateBoundaryLength(1:iSegment))));
    FiltedCrustalThickness=fastSmooth(Crust.Thickness(iStart:iEnd),EffectiveSmoothingWidth,3,1);
    CrustalThickness=[CrustalThickness,FiltedCrustalThickness];
    iStart=1+iEnd;
end
Crust.Thickness=CrustalThickness;

%% Reoriented Plate Boundaires with Spreading in x-Direction for Easier Accretion Calculation 

CompensateLength=200; % distance from edges of plate boundary segments to model boundaries 
PlateBoundary_x=[0]; % initialize reoriented plate boundaries
PlateBoundary_y=[0];
% % reorient plate boundary segments
for iSegment=2:numel(Geometry.PlateBoundaryLength)+1;
    PlateBoundary_x(iSegment)=PlateBoundary_x(iSegment-1)+Geometry.PlateBoundaryLength(iSegment-1).*...
        cosd(Geometry.PlateBoundaryAngleDeg(iSegment-1)-SpreadingDirection);
    PlateBoundary_y(iSegment)=PlateBoundary_y(iSegment-1)+Geometry.PlateBoundaryLength(iSegment-1).*...
        sind(Geometry.PlateBoundaryAngleDeg(iSegment-1)-SpreadingDirection);
end
PlateBoundary_x=PlateBoundary_x+abs(min(PlateBoundary_x))*(min(PlateBoundary_x)<0)+CompensateLength;
PlateBoundary_y=PlateBoundary_y+abs(min(PlateBoundary_y))*(min(PlateBoundary_y)<0)+TrimDistance;
ModelBoundary_x=[0,max(PlateBoundary_x)+CompensateLength];
ModelBoundary_y=[min(PlateBoundary_y),max(PlateBoundary_y)];
% % grid reoriented geometry
SamplingSize_x=ceil(diff(ModelBoundary_x)/CellSize)+1;
SamplingSize_y=ceil(diff(ModelBoundary_y)/CellSize)+1;
Vector_x=unique([linspace(ModelBoundary_x(1),ModelBoundary_x(2),SamplingSize_x),PlateBoundary_x]); % include boundary turning points in the nodes
Vector_y=unique([linspace(ModelBoundary_y(1),ModelBoundary_y(2),SamplingSize_y),PlateBoundary_y]); 
[Grid_x,Grid_y]=meshgrid(Vector_x,Vector_y); % grid for accretion calculation

%% Accretion in x-Direction

iTracker=0; % waitbar progress tracker
nTracker=(size(Grid_x,1)-1)*(size(Grid_x,2)-1); % waitbar progress tracker
WaitBar=waitbar(0,'Please wait...');

warning off all;
% % generate tiles and accret crust in each tile along spreading direction
for iRow=1:size(Grid_x,1)-1; % work on every row
    for iCol=1:size(Grid_x,2)-1; % work on every column
        Tile(iRow,iCol).x=[Grid_x(iRow,iCol),Grid_x(iRow+1,iCol),Grid_x(iRow+1,iCol+1),Grid_x(iRow,iCol+1)]; % x coordinates of the tile (counterclockwise)
        Tile(iRow,iCol).y=[Grid_y(iRow,iCol),Grid_y(iRow+1,iCol),Grid_y(iRow+1,iCol+1),Grid_y(iRow,iCol+1)]; % y coordinates of the tile (counterclockwise)
        Tile(iRow,iCol).Center_x=mean(Tile(iRow,iCol).x); % x coordinate of the center of the tile
        Tile(iRow,iCol).Center_y=mean(Tile(iRow,iCol).y); % y coordinate of the center of the tile
        Tile(iRow,iCol).Area=abs(clockwiseAreaCalc(Tile(iRow,iCol).x,Tile(iRow,iCol).y)); % area of the tile
        Tile(iRow,iCol).CrustalThickness=0; % default value for crust thickness in the tile
        % % find the corresponding point on the segment
        for iSegment=1:numel(PlateBoundary_x)-1;
            if (PlateBoundary_y(iSegment)>=Tile(iRow,iCol).Center_y)&&(PlateBoundary_y(iSegment+1)<=Tile(iRow,iCol).Center_y);
                iTargetSegment=iSegment;
            end
        end
        Tile(iRow,iCol).OriginPoint=[(Tile(iRow,iCol).Center_y-PlateBoundary_y(iTargetSegment+1))*(PlateBoundary_x(iTargetSegment+1)-...
            PlateBoundary_x(iTargetSegment))/(PlateBoundary_y(iTargetSegment+1)-PlateBoundary_y(iTargetSegment))+...
            PlateBoundary_x(iTargetSegment+1),Tile(iRow,iCol).Center_y]; % coordinate of the corresponding point on boundaries
        
        TilePolygon.x=Tile(iRow,iCol).x; % tile polygon x coordinate, for distance determination
        TilePolygon.y=Tile(iRow,iCol).y; % tile polygon y coordinate 
        for iSegment=1:numel(PlateBoundary_x)-1; % determine spatial relation between the tile and every boundary segment
            % % segment polygon coordinates, for distance determination
            SegmentPolygon.x=[PlateBoundary_x(iSegment:iSegment+1),fliplr(PlateBoundary_x(iSegment:iSegment+1))]; % segment polygon x coordinate
            SegmentPolygon.y=[PlateBoundary_y(iSegment:iSegment+1),fliplr(PlateBoundary_y(iSegment:iSegment+1))]; % segment polygon y coordinate
            Distance=minDistBetweenTwoPolygons(SegmentPolygon,TilePolygon); % distance between tile polygon and segment polygon
            switch Geometry.PlateBoundaryType(iSegment); % set the accretion width based on boundary type
                case 1;
                    AccretWidth=AccretWidth_Ridge;
                case 2;
                    AccretWidth=AccretWidth_Transform;
                case 3;
                    AccretWidth=mean(AccretWidth_Transform,AccretWidth_Ridge);
            end
            if Distance<AccretWidth; % if the tile is withing the accretion box
                % % assign each point on the tile to the segment
                [SegmentNumber1,AlongSegmentCoordinate1,DistanceToSegment1]=assignSegment(Geometry,TilePolygon.x(1),TilePolygon.y(1));
                [SegmentNumber2,AlongSegmentCoordinate2,DistanceToSegment2]=assignSegment(Geometry,TilePolygon.x(2),TilePolygon.y(2));
                [SegmentNumber3,AlongSegmentCoordinate3,DistanceToSegment3]=assignSegment(Geometry,TilePolygon.x(3),TilePolygon.y(3));
                [SegmentNumber4,AlongSegmentCoordinate4,DistanceToSegment4]=assignSegment(Geometry,TilePolygon.x(4),TilePolygon.y(4));
                Assign=[AlongSegmentCoordinate1,AlongSegmentCoordinate2,AlongSegmentCoordinate3,AlongSegmentCoordinate4]; 
                % % starting and ending coordinate of the contributing portion
                % % determine the portion of the boundary affecting the tile
                iStart=ceil(min(Assign)); 
                iEnd=ceil(max(Assign));
                if iStart<=0;
                    iStart=1;
                end
                if iEnd>numel(Crust.Thickness);
                    iEnd=numel(Crust.Thickness);
                end
                
                switch Geometry.PlateBoundaryType(iSegment);
                    case 1; % if the tile is affected by a ridge segment
                        % % only consider tiles on side of ridge
                        if ((Tile(iRow,iCol).Center_y-PlateBoundary_y(iSegment))*(Tile(iRow,iCol).Center_y-PlateBoundary_y(iSegment+1)))<=0;
                            if Tile(iRow,iCol).Center_x<=Tile(iRow,iCol).OriginPoint(1,1); % if the tile is to the left of the ridge
                                AccretZoneBoundary=Tile(iRow,iCol).OriginPoint(1,1)-AccretWidth; % x coordinate of the left boundary of the accretion box
                                % % if the tile lies across the boundary
                                if (Tile(iRow,iCol).x(1)<AccretZoneBoundary)&&(Tile(iRow,iCol).x(3)>AccretZoneBoundary);
                                    % % crustal thickness in the tile with reduced size
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*... 
                                        abs(clockwiseAreaCalc([AccretZoneBoundary,AccretZoneBoundary,Tile(iRow,iCol).x(3:4)],Tile(iRow,iCol).y))/...
                                        (AccretWidth*abs(Tile(iRow,iCol).y(2)-Tile(iRow,iCol).y(1)));
                                else % if the tile is within the accretion box
                                    % % crustal thickness in the tile
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*...
                                        Tile(iRow,iCol).Area/(AccretWidth*abs(Tile(iRow,iCol).y(2)-Tile(iRow,iCol).y(1)));
                                end
                            else % if the tile is to the right of the ridge
                                AccretZoneBoundary=Tile(iRow,iCol).OriginPoint(1,1)+AccretWidth; % x coordinate of the right boundary of the accretion box
                                % % if the tile lies across the boundary
                                if (Tile(iRow,iCol).x(1)<AccretZoneBoundary)&&(Tile(iRow,iCol).x(3)>AccretZoneBoundary);
                                    % % crustal thickness in the tile with reduced size
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*... 
                                        abs(clockwiseAreaCalc([Tile(iRow,iCol).x(1:2),AccretZoneBoundary,AccretZoneBoundary],Tile(iRow,iCol).y))/...
                                        (AccretWidth*abs(Tile(iRow,iCol).y(2)-Tile(iRow,iCol).y(1)));
                                else % if the tile is within the accretion box
                                    % % crustal thickness in the tile
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*...
                                        Tile(iRow,iCol).Area/(AccretWidth*abs(Tile(iRow,iCol).y(2)-Tile(iRow,iCol).y(1)));
                                end
                            end
                        end
                    case 2; % if the tile is affected by a transform segment
                        % % only consider tiles on sides of transform
                        if ((Tile(iRow,iCol).Center_x-PlateBoundary_x(iSegment+1))*(Tile(iRow,iCol).Center_x-PlateBoundary_x(iSegment)))<=0; 
                            if Tile(iRow,iCol).Center_y<PlateBoundary_y(iSegment); % if the tile is below the transform
                                AccretZoneBoundary=PlateBoundary_y(iSegment)-AccretWidth; % y coordinate of the lower boundary of the accretion box
                                if (Tile(iRow,iCol).y(1)<AccretZoneBoundary)&&(Tile(iRow,iCol).y(3)>AccretZoneBoundary); % if the tile lies across the boundary
                                    % % crustal thickness in the tile with reduced size
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*... 
                                        abs(clockwiseAreaCalc(Tile(iRow,iCol).x,[AccretZoneBoundary,Tile(iRow,iCol).y(2:3),AccretZoneBoundary]))/...
                                        (AccretWidth*abs(Tile(iRow,iCol).x(4)-Tile(iRow,iCol).x(1)));
                                else % if the tile is within the accretion box
                                    % % crustal thickness in the tile
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*...
                                        Tile(iRow,iCol).Area/(AccretWidth*abs(Tile(iRow,iCol).x(4)-Tile(iRow,iCol).x(1))); 
                                end
                            elseif Tile(iRow,iCol).Center_y>PlateBoundary_y(iSegment); % if the tile is above the transform
                                AccretZoneBoundary=PlateBoundary_y(iSegment)+AccretWidth; % y coordinate of the upper boundary of the accretion box
                                if (Tile(iRow,iCol).y(1)<AccretZoneBoundary)&&(Tile(iRow,iCol).y(3)>AccretZoneBoundary); % if the tile lies across the boundary
                                    % % crustal thickness in the tile with reduced size
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*... 
                                        abs(clockwiseAreaCalc(Tile(iRow,iCol).x,[Tile(iRow,iCol).y(1),AccretZoneBoundary,AccretZoneBoundary,Tile(iRow,iCol).y(4)]))/...
                                        (AccretWidth*abs(Tile(iRow,iCol).x(4)-Tile(iRow,iCol).x(1)));
                                else % if the tile is within the accretion box
                                    % % crustal thickness in the tile
                                    Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*...
                                        Tile(iRow,iCol).Area/(AccretWidth*abs(Tile(iRow,iCol).x(4)-Tile(iRow,iCol).x(1))); 
                                end
                            end
                        end
                    case 3; % if the tile is affected by an oblique segment
                        % % coordinate of the boundaries of the oblique accretion box
                        ObliqueAccretBoundary1=[PlateBoundary_x(iSegment)+AccretWidth*sind(abs(Geometry.PlateBoundaryAngleDeg(iSegment))),...
                            PlateBoundary_y(iSegment)+AccretWidth*cosd(abs(Geometry.PlateBoundaryAngleDeg(iSegment))),...
                            PlateBoundary_x(iSegment)-AccretWidth*sind(abs(Geometry.PlateBoundaryAngleDeg(iSegment))),...
                            PlateBoundary_y(iSegment)-AccretWidth*cosd(abs(Geometry.PlateBoundaryAngleDeg(iSegment)))];
                        ObliqueAccretBoundary2=[PlateBoundary_x(iSegment+1)+AccretWidth*sind(abs(Geometry.PlateBoundaryAngleDeg(iSegment))),...
                            PlateBoundary_y(iSegment+1)+AccretWidth*cosd(abs(Geometry.PlateBoundaryAngleDeg(iSegment))),...
                            PlateBoundary_x(iSegment+1)-AccretWidth*sind(abs(Geometry.PlateBoundaryAngleDeg(iSegment))),...
                            PlateBoundary_y(iSegment+1)-AccretWidth*cosd(abs(Geometry.PlateBoundaryAngleDeg(iSegment)))];
                        % % only consider tiles on side of oblique segment
                        if (Tile(iRow,iCol).Center_x>=(Tile(iRow,iCol).Center_y-ObliqueAccretBoundary1(4))*(ObliqueAccretBoundary1(3)-ObliqueAccretBoundary1(1))...
                                /(ObliqueAccretBoundary1(4)-ObliqueAccretBoundary1(2))+ObliqueAccretBoundary1(3))...
                                &&(Tile(iRow,iCol).Center_x<=(Tile(iRow,iCol).Center_y-ObliqueAccretBoundary2(4))*(ObliqueAccretBoundary2(3)-ObliqueAccretBoundary2(1))...
                                /(ObliqueAccretBoundary2(4)-ObliqueAccretBoundary2(2))+ObliqueAccretBoundary2(3));
                            % % calculate the projection of the tile on the oblique segment
                            Matrix=[];
                            for iCorner=1:4;
                                Matrix1=[PlateBoundary_y(iSegment)-PlateBoundary_y(iSegment+1), -PlateBoundary_x(iSegment)+PlateBoundary_x(iSegment+1);...
                                    PlateBoundary_x(iSegment)-PlateBoundary_x(iSegment+1),PlateBoundary_y(iSegment)-PlateBoundary_y(iSegment+1)];
                                Matrix2=[PlateBoundary_x(iSegment+1)*(PlateBoundary_y(iSegment)-PlateBoundary_y(iSegment+1))-PlateBoundary_y(iSegment+1)*...
                                    (PlateBoundary_x(iSegment)-PlateBoundary_x(iSegment+1));...
                                    Tile(iRow,iCol).x(iCorner)*(PlateBoundary_x(iSegment)-PlateBoundary_x(iSegment+1))+Tile(iRow,iCol).y(iCorner)*...
                                    (PlateBoundary_y(iSegment)-PlateBoundary_y(iSegment+1))];
                                Matrix=[Matrix,Matrix1\Matrix2];
                            end
                            Projection1=sqrt((Matrix(1,1)-Matrix(1,3))^2+(Matrix(2,1)-Matrix(2,3))^2);
                            Projection2=sqrt((Matrix(1,2)-Matrix(1,4))^2+(Matrix(2,2)-Matrix(2,4))^2);
                            Projection=max(Projection1,Projection2); % the length of tile projection on oblique segment, use the longer one.
                            Tile(iRow,iCol).CrustalThickness=Tile(iRow,iCol).CrustalThickness+mean(Crust.Thickness(iStart:iEnd))*...
                                Tile(iRow,iCol).Area/(AccretWidth*Projection); % crustal thickness in the tile
                        end
                end              
            end
        end
        iTracker=iTracker+1;
        waitbar(iTracker/nTracker,WaitBar,sprintf('Working on Tile (%g,%g)',iRow,iCol));
    end
end
delete(WaitBar);
warning on all

% % accummulate crustal thickness
for iRow=1:size(Grid_x,1)-1;
    LeftPlate=[];
    RightPlate=[];
    for iCol=1:size(Grid_x,2)-1;
        if Tile(iRow,iCol).Center_x<=Tile(iRow,iCol).OriginPoint(1,1);
            LeftPlate=[LeftPlate;iCol]; % tiles to the left of the boundaries
        end
        if Tile(iRow,iCol).Center_x>Tile(iRow,iCol).OriginPoint(1,1);
            RightPlate=[RightPlate;iCol]; % tiles to the right of the boundaries
        end
    end
    for iLeft=1:LeftPlate(end-1); % summation to the left
        Tile(iRow,LeftPlate(end)-iLeft).CrustalThickness=Tile(iRow,LeftPlate(end)-iLeft+1).CrustalThickness+Tile(iRow,LeftPlate(end)-iLeft).CrustalThickness;
    end
    for iRight=RightPlate(2):RightPlate(end); % summation to the right
        Tile(iRow,iRight).CrustalThickness=Tile(iRow,iRight-1).CrustalThickness+Tile(iRow,iRight).CrustalThickness;
    end
    for iCol=1:size(Grid_x,2)-1; % store the information
        SeafloorHistory.Thickness(iRow,iCol)=Tile(iRow,iCol).CrustalThickness;
        SeafloorHistory.TileCenter_x(iRow,iCol)=Tile(iRow,iCol).Center_x;
        SeafloorHistory.TileCenter_y(iRow,iCol)=Tile(iRow,iCol).Center_y;
    end
end

%% Averaging Over Domains

RidgeDomainIndex=[];
TransformDomainIndex_x=[];
TransformDomainIndex_y=[];
ObliqueDomainIndex=[];
for iSegment=1:numel(Geometry.PlateBoundaryType);
    switch Geometry.PlateBoundaryType(iSegment);
        case 1
            iRidge=find(((SeafloorHistory.TileCenter_y(:,1)-PlateBoundary_y(iSegment)).*(SeafloorHistory.TileCenter_y(:,1)-PlateBoundary_y(iSegment+1)))<=0);
            RidgeDomainIndex=[RidgeDomainIndex;iRidge];
        case 2
            iTransform_x=find(((SeafloorHistory.TileCenter_y(:,1)-PlateBoundary_y(iSegment)+AccretWidth_Transform).*...
                (SeafloorHistory.TileCenter_y(:,1)-PlateBoundary_y(iSegment)-AccretWidth_Transform))<=0);
            iTransform_y=find(((SeafloorHistory.TileCenter_x(1,:)-PlateBoundary_x(iSegment+1)).*(SeafloorHistory.TileCenter_x(1,:)-PlateBoundary_x(iSegment)))<=0);
            TransformDomainIndex_x=[TransformDomainIndex_x;iTransform_x];
            TransformDomainIndex_y=[TransformDomainIndex_y;iTransform_y'];
        case 3
            iOblique=find(((SeafloorHistory.TileCenter_y(:,1)-PlateBoundary_y(iSegment)).*(SeafloorHistory.TileCenter_y(:,1)-PlateBoundary_y(iSegment+1)))<=0);
            ObliqueDomainIndex=[ObliqueDomainIndex;iOblique];
    end
end
AverageWholdDomain=mean(mean(SeafloorHistory.Thickness));
AverageRidgeDomain=mean(mean(SeafloorHistory.Thickness(RidgeDomainIndex,:)));
AverageTransformDomain=mean(mean(SeafloorHistory.Thickness(TransformDomainIndex_x,TransformDomainIndex_y)));
AverageObliqueDomain=mean(mean(SeafloorHistory.Thickness(ObliqueDomainIndex,:)));

AverageWidth=10;
AverageSegmentCrust=[];
iStart=1;
for iSegment=1:numel(Geometry.PlateBoundaryLength);
    SegmentSamplingSize=ceil(Geometry.PlateBoundaryLength(iSegment));
    SegmentSampling_x=linspace(PlateBoundary_x(iSegment),PlateBoundary_x(iSegment+1),SegmentSamplingSize);
    SegmentSampling_y=linspace(PlateBoundary_y(iSegment),PlateBoundary_y(iSegment+1),SegmentSamplingSize);
    for i=1:SegmentSamplingSize;
        iTargetSegment_x=find(((SeafloorHistory.TileCenter_y(:,1)-(SegmentSampling_y(i)-AverageWidth)).*...
            (SeafloorHistory.TileCenter_y(:,1)-(SegmentSampling_y(i)+AverageWidth)))<=0);
        iTargetSegment_y=find(((SeafloorHistory.TileCenter_x(1,:)-(SegmentSampling_x(i)-AverageWidth)).*...
            (SeafloorHistory.TileCenter_x(1,:)-(SegmentSampling_x(i)+AverageWidth)))<=0);
        CurrentAverageSegmentCrust=sum(sum(SeafloorHistory.Thickness(iTargetSegment_x,iTargetSegment_y)))...
            /(numel(SeafloorHistory.Thickness(iTargetSegment_x,iTargetSegment_y)));
        AverageSegmentCrust=[AverageSegmentCrust,CurrentAverageSegmentCrust];
    end
    SeafloorHistory.AverageSegmentCrust(iSegment)=mean(AverageSegmentCrust(iStart:end));
    iStart=1+ceil(SegmentSamplingSize);
end

SeafloorHistory.AverageSegmentCrust=AverageSegmentCrust;
SeafloorHistory.AverageWholdDomain=AverageWholdDomain;
SeafloorHistory.AverageRidgeDomain=AverageRidgeDomain;
SeafloorHistory.AverageTransformDomain=AverageTransformDomain;
SeafloorHistory.AverageObliqueDomain=AverageObliqueDomain;
SeafloorHistory.AccretWidth_Ridge=AccretWidth_Ridge;
SeafloorHistory.AccretWidth_Transform=AccretWidth_Transform;
SeafloorHistory.SmoothingWidth=SmoothingWidth;
SeafloorHistory.CellSize=CellSize;

%% Crustal Thickness Map Plot

figure(indFigure); clf; hold on;
pcolor(SeafloorHistory.TileCenter_x,SeafloorHistory.TileCenter_y,SeafloorHistory.Thickness); shading interp; colorbar;
plot(PlateBoundary_x,PlateBoundary_y,'w.-');
axis equal; axis tight;
caxis([0,15]);
title(sprintf('Average=%g, Ridge=%g, Transform=%g, Oblique=%g',AverageWholdDomain,AverageRidgeDomain,AverageTransformDomain,AverageObliqueDomain));

figure(indFigure+1); clf; hold on; box on;
plot([1:1:ceil(sum(Geometry.PlateBoundaryLength))],SeafloorHistory.AverageSegmentCrust,'linewidth',2);
for iSegment=1:numel(Geometry.PlateBoundaryLength);
    plot([1 ,1]*sum(Geometry.PlateBoundaryLength(1:iSegment)),[0,1]*12,'k','linewidth',1);
end
ylim([0,12]);
xlim([0,sum(Geometry.PlateBoundaryLength)]);
ylabel('Accumulated Crust (km)');
title('Averaged Crustal Thickness [km');

plotCrustalThicknessMapIn3D(SeafloorHistory,PlateBoundary_x,PlateBoundary_y,ModelBoundary_x,ModelBoundary_y,indFigure+2);