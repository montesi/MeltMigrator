function iExtraction=extractionDetermination(Geometry,LocalSwath,iPolygon,ExtractionWidth,ExtractionDepth,ExtractionSlope)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractionDetermination.m
% Determine if the patch of melt is extracted
% Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.PlateBoundaryType : Indicator of plate boundary type. 1: ridge, 2: transform, 3: oblique segment
%       |......
%   LocalSwath
%       |.Polygon(ind)
%           |.x,.y          : Coordinates of polygon outline [km]
%       |.PolygonCenter_z   : Polygon center depth [km]
%       |.PolygonSlope      : Slope at polygon centers
%       |......
%   iPolygon                : Index of current polygon
%   ExtractionWidth         : Melt extraction width [km]
%   ExtractionDepth         : Melt extraction depth [km]
%   ExtractionSlope         : Melt transport critical slope
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   iExtraction             : Indicator of extraction. 0: don't extract, transport to next polygon, 1: extract as crust, 2: extract as cryptic crust 
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   PlateBoundaryPolygon    : Plate boundary polygon for distance determination
%   RidgePolygon            : Ridge polygon for distance determination
%   LocalSwathPolygon       : Local polygon for distance determination
%   DistanceToSegment       : Distance between current polygon and plate boundary or ridge polygon 
%   iRidge                  : Index of ridge segments
%   iRidgePolygon           : Index of ridge polygons
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   minDistBetweenTwoPolygons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extraction Determination

% % construct plate boundary polygon and swath polygon for distance determination 
PlateBoundaryPolygon.x=[Geometry.PlateBoundary.x,fliplr(Geometry.PlateBoundary.x)];
PlateBoundaryPolygon.y=[Geometry.PlateBoundary.y,fliplr(Geometry.PlateBoundary.y)];
LocalSwathPolygon.x=LocalSwath.Polygon(iPolygon).x;
LocalSwathPolygon.y=LocalSwath.Polygon(iPolygon).y;
DistanceToSegment=minDistBetweenTwoPolygons(PlateBoundaryPolygon,LocalSwathPolygon);

if ~isnan(ExtractionWidth); % if extraction width is defined
    if DistanceToSegment<=ExtractionWidth; % if polygon is within extraction width
        if LocalSwath.PolygonCenter_z(iPolygon)<=ExtractionDepth; % if polygon is above extraction depth
            iExtraction=1; % extract
        else % if polygon is below extraction depth
            if LocalSwath.PolygonSlope(iPolygon)>=ExtractionSlope; % if polygon slope is larger than critical slope
                iExtraction=0; % transport to next polygon
            else % if polygon slope is smaller than critical slope
                iExtraction=2; % do not extract
            end
        end
    else % if polygon is not within extraction width
        if LocalSwath.PolygonSlope(iPolygon)>=ExtractionSlope; % if polygon slope is larger than critical slope
            iExtraction=0; % transport to next tile
        else % if polygon slope is smaller than critical slope
            iExtraction=2; % do not extract
        end
    end
else % if extraction width is not defined, melt ends up with ridge segment
    % define ridge segment polygons for distance determination
    iRidge=find(Geometry.PlateBoundaryType==1);
    for iRidgePolygon=1:numel(iRidge);
        RidgePolygon.x=[Geometry.PlateBoundary.x(iRidge(iRidgePolygon):iRidge(iRidgePolygon)+1),...
            fliplr(Geometry.PlateBoundary.x(iRidge(iRidgePolygon):iRidge(iRidgePolygon)+1))]; 
        RidgePolygon.y=[Geometry.PlateBoundary.y(iRidge(iRidgePolygon):iRidge(iRidgePolygon)+1),...
            fliplr(Geometry.PlateBoundary.y(iRidge(iRidgePolygon):iRidge(iRidgePolygon)+1))];
        DistanceToSegment(iRidgePolygon)=minDistBetweenTwoPolygons(RidgePolygon,LocalSwathPolygon);
    end
    if min(DistanceToSegment)==0; % extract melt only if swath ends on ridge segment 
        iExtraction=1;
    else
        if LocalSwath.PolygonSlope(iPolygon)>=ExtractionSlope;
            iExtraction=0;
        else
            iExtraction=2;
        end
    end
end