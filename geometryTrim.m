function Geometry=geometryTrim(Geometry,TrimDistance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometryTrim.m;
% Trim the edge of the model geometry to avoid boundary effects
% Find the intersection between plate boundary segments and new model boundaries
% Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry             
%       |.ModelBoundary     : Model geometry corner coordinates [km]
%       |.PlateBoundary     : Plate boundary segment coordinates [km]
%       |......
%   TrimDistance            : Distance to trim from edges [km]
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Geometry             
%       |.ModelBoundary     : Trimmed model geometry corner coordinates [km]
%       |.PlateBoundary     : Trimmed plate boundary segment coordinates [km]
%--------------------------------------------------------------------------
% INTERNAL ---------------------------------------------------------------- 
%   PlateSegment1           : Coordinate for first segment of plate boundary [km]
%   PlateSegment2           : Coordinate for last segment of plate boundary [km]
%   BoundarySegment         : Coordinates for new boundaries of model [km]
%   SegmentIntersection1    : Intersection relation between first plate boundary segment and model boundaries
%   SegmentIntersection2    : Intersection relation between last plate boundary segment and model boundaries
%   iInt1                   : Index for actual intersection between first plate boundary segment and model boundary
%   iInt2                   : Index for actual intersection between last plate boundary segment and model boundary
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS 
%   lineSegmentIntersect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

if TrimDistance~=0;   
    Geometry.ModelBoundary.x=[Geometry.ModelBoundary.x(1)+TrimDistance,Geometry.ModelBoundary.x(2)-TrimDistance]; 
    Geometry.ModelBoundary.y=[Geometry.ModelBoundary.y(1)+TrimDistance,Geometry.ModelBoundary.y(2)-TrimDistance];
    
    PlateSegment1=[Geometry.PlateBoundary.x(1),Geometry.PlateBoundary.y(1),Geometry.PlateBoundary.x(2),Geometry.PlateBoundary.y(2)];
    PlateSegment2=[Geometry.PlateBoundary.x(end-1),Geometry.PlateBoundary.y(end-1),Geometry.PlateBoundary.x(end),Geometry.PlateBoundary.y(end)];
    BoundarySegment=[Geometry.ModelBoundary.x(1),Geometry.ModelBoundary.y(1),Geometry.ModelBoundary.x(1),Geometry.ModelBoundary.y(2);...
        Geometry.ModelBoundary.x(1),Geometry.ModelBoundary.y(2),Geometry.ModelBoundary.x(2),Geometry.ModelBoundary.y(2);...
        Geometry.ModelBoundary.x(2),Geometry.ModelBoundary.y(2),Geometry.ModelBoundary.x(2),Geometry.ModelBoundary.y(1);...
        Geometry.ModelBoundary.x(2),Geometry.ModelBoundary.y(1),Geometry.ModelBoundary.x(1),Geometry.ModelBoundary.y(1)];
    SegmentIntersection1=lineSegmentIntersect(PlateSegment1,BoundarySegment);
    SegmentIntersection2=lineSegmentIntersect(PlateSegment2,BoundarySegment);
    iInt1=find(SegmentIntersection1.intAdjacencyMatrix==1);
    iInt2=find(SegmentIntersection2.intAdjacencyMatrix==1);
    if (length(iInt1)==1)&&(length(iInt2)==1);
        Geometry.PlateBoundary.x(1)=SegmentIntersection1.intMatrixX(iInt1);
        Geometry.PlateBoundary.y(1)=SegmentIntersection1.intMatrixY(iInt1);
        Geometry.PlateBoundary.x(end)=SegmentIntersection2.intMatrixX(iInt2);
        Geometry.PlateBoundary.y(end)=SegmentIntersection2.intMatrixY(iInt2);
    elseif (length(iInt1)==2)&&(length(Geometry.PlateBoundary.x)==2);
        Geometry.PlateBoundary.x(1)=SegmentIntersection1.intMatrixX(iInt1(1));
        Geometry.PlateBoundary.y(1)=SegmentIntersection1.intMatrixY(iInt1(1));
        Geometry.PlateBoundary.x(end)=SegmentIntersection1.intMatrixX(iInt1(2));
        Geometry.PlateBoundary.y(end)=SegmentIntersection1.intMatrixY(iInt1(2));
    else
        disp('!!! Discard distance is not proper');
    end
end

return