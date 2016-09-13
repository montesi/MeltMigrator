function [SegmentNumber,AlongSegmentCoordinate,DistanceToAllSegment]=assignSegment(Geometry,Position_x,Position_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assignSegment.m
% Assign nearest plate boundary segment with points on melt trajectory
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.PlateBoundaryLength : Length of each plate boundary segment [km]
%       |.PlateBoundaryAngleRad : Angle of each plate boundary segment [rad]
%       |......
%   Position
%       |_x, _y             : Coordinates of point [km]
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   SegmentNumber           : Index of plate boundary segment the line ends of 
%   AlongSegmentCoordinate  : Along plate boundary coordinates of point projection [km] 
%   DistanceToSegment       : Distance between point and plate boundary [km] 
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   iPosition               : Index for points
%   Diff_x                  : x-direction distance between point and plate boundary corners [km] 
%   Diff_y                  : y-direction distance between point and plate boundary corners [km]
%   LengthAlongSegment      : Coordinate of point projection along each segment [km]
%   DistanceToAllSegment    : Distance to all plate boundary segments [km]
%   IsIn                    : Determinator for point projection within plate boundary segments 
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

warning off MATLAB:divideByZero; 

for iPosition=1:length(Position_x);
    Diff_x=Position_x(iPosition)-Geometry.PlateBoundary.x(1:end-1);
    Diff_y=Position_y(iPosition)-Geometry.PlateBoundary.y(1:end-1);
    LengthAlongSegment=Diff_x.*cos(Geometry.PlateBoundaryAngleRad)+Diff_y.*sin(Geometry.PlateBoundaryAngleRad);
    DistanceToAllSegment=-Diff_x.*sin(Geometry.PlateBoundaryAngleRad)+Diff_y.*cos(Geometry.PlateBoundaryAngleRad);
    IsIn=(LengthAlongSegment>=0).*(LengthAlongSegment<=Geometry.PlateBoundaryLength);
    if max(IsIn)==0; % look for point projection outside plate boundary
        DistanceToAllSegment=((Geometry.PlateBoundary.x-Position_x(iPosition)).^2+(Geometry.PlateBoundary.y-Position_y(iPosition)).^2).^(1/2);
        [DistanceToSegment(iPosition),SegmentNumber(iPosition)]=min(abs(DistanceToAllSegment));
        AlongSegmentCoordinate(iPosition)=sum(Geometry.PlateBoundaryLength([1:SegmentNumber(iPosition)-1]));
        SegmentNumber(iPosition)=min(SegmentNumber(iPosition),numel(Geometry.PlateBoundaryLength_x));
    else
        [DistanceToSegment(iPosition),SegmentNumber(iPosition)]=min(abs(DistanceToAllSegment)./IsIn);
        AlongSegmentCoordinate(iPosition)=sum(Geometry.PlateBoundaryLength([1:(SegmentNumber(iPosition)-1)]))+...
        Diff_x(SegmentNumber(iPosition)).*cos(Geometry.PlateBoundaryAngleRad(SegmentNumber(iPosition)))+...
        Diff_y(SegmentNumber(iPosition)).*sin(Geometry.PlateBoundaryAngleRad(SegmentNumber(iPosition)));
    end
end

warning on MATLAB:divideByZero;

return