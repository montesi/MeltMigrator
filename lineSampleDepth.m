function [LineUpdated,OnSegmentInfo]=lineSampleDepth(Geometry,Grid,Res,Lid,Trajectory,ExtendToSegmentDistance,LidDepthLimit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lineSampleDepth.m
% Trucate and resample melt line
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.PlateBoundaryLength : Length of each plate boundary segment [km]
%       |......
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates [km]
%   Res
%       |.nMeltSwath        : Melt swath sampling size [km]
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
%   ExtendToSegmentDistance : When line ends get in range, extend line to plate boundary segment [km] 
%   LidDepthLimit           : Lower limit for permeability barrier [km]
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   LineUpdated
%       |(1,:) NewSample_x  : x-coordinate of resampled line [km]
%       |(2,:) NewSample_y  : y-coordinate of resampled line [km]
%       |(3,:) NewSample_z  : z-coordinate of resampled line [km]
%       |(4,:) NewSample_Distance : Distance along resampled line [km]
%       |(5,:) Sample_Slope : Slope of resampled line [km]
%   OnSegmentInfo
%       |(1) OnSegmentPoint_x : x-coordinate of resampled line end on segment [km]
%       |(2) OnSegmentPoint_y : y-coordinate of resampled line end on segment [km]
%       |(3) SegmentNumber  : Index of plate boundary segment the line ends of 
%       |(4) AlongSegmentCoordinate : Along plate boundary coordinates of polygon projection [km] 
%       |(5) DistanceToSegment : Distance between polygon and plate boundary [km] 
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   iWithinMelting          : Index for points on melt line with temperature higher than solidus
%   iBelowDepthLimit        : Index for points on melt line below lower limit of lid 
%   Sample_Distance         : Distance along sampled line [km]
%   Sample_x                : x-coordinate of sampled line [km]
%   Sample_y                : y-coordinate of sampled line [km]
%   Sample_z                : z-coordinate of sampled line [km]
%   DistanceAtCrossing      : Distance along sampled line when crossing lower lime of lid 
%   OnSegmentPoint_z        : z-coordinate of resampled line end on segment [km] 
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   assignSegment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Line Resampling

LidDepthLimit=min(LidDepthLimit,max(Geometry.ModelBoundary.z)); % set lower limit of lid

% consider only part of surface above solidus 
iWithinMelting=find(Trajectory.T>Trajectory.T_Solidus); % line is within melting region
if isempty(iWithinMelting); % line is not in melting region
    if numel(Trajectory)==1; % only one point on melt trajectory 
        Sample_Distance=Trajectory.Distance;
        Sample_x=Trajectory.x;
        Sample_y=Trajectory.y;
        Sample_z=Trajectory.Depth;
        Sample_Slope=Trajectory.Slope;
    else % more than one point on melt trajectory: sample mid-point of trajectory
        Sample_Distance=Trajectory.Distance(ceil(numel(Trajectory.T)./2));
        Sample_x=interp1(Trajectory.Distance,Trajectory.x,Sample_Distance);
        Sample_y=interp1(Trajectory.Distance,Trajectory.y,Sample_Distance);
        Sample_z=interp1(Trajectory.Distance,Trajectory.Depth,Sample_Distance);
        Sample_Slope=interp1(Trajectory.Distance,Trajectory.Slope,Sample_Distance);
    end
else % a part of line is within melting region
    iBelowDepthLimit=max(find(Trajectory.Depth(iWithinMelting)>LidDepthLimit)); % find points below lower limit
    if isempty(iBelowDepthLimit) % line is above lower limit
        Sample_Slope=Trajectory.Slope(iWithinMelting);
        Sample_Distance=Trajectory.Distance(iWithinMelting);
    else % a part of line is below lower limit
        if iWithinMelting(iBelowDepthLimit)==max(iWithinMelting) % all points below lower limit
            Sample_Slope=Trajectory.Slope(end);
            Sample_Distance=Trajectory.Distance(end);
        else % line crosses lower limit
            DistanceAtCrossing=Trajectory.Distance(iWithinMelting(iBelowDepthLimit))+...
                (Trajectory.Distance(iWithinMelting(iBelowDepthLimit)+1)-Trajectory.Distance(iWithinMelting(iBelowDepthLimit)))*...
                (LidDepthLimit-Trajectory.Depth(iWithinMelting(iBelowDepthLimit)))/...
                (Trajectory.Depth(iWithinMelting(iBelowDepthLimit)+1)-Trajectory.Depth(iWithinMelting(iBelowDepthLimit)));
            Sample_Slope=[interp1(Trajectory.Distance,Trajectory.Slope,DistanceAtCrossing);Trajectory.Slope(iWithinMelting(iBelowDepthLimit)+1:end)];
            Sample_Distance=[DistanceAtCrossing;Trajectory.Distance(iWithinMelting(iBelowDepthLimit)+1:end)];
        end
    end
    Sample_x=interp1(Trajectory.Distance,Trajectory.x,Sample_Distance);
    Sample_y=interp1(Trajectory.Distance,Trajectory.y,Sample_Distance);
    Sample_z=interp1(Trajectory.Distance,Trajectory.Depth,Sample_Distance);
end

%% Points Addition

% % extend current line to plate boundary segments
if ~isempty(Sample_x);
    [SegmentNumber,AlongSegmentCoordinate,DistanceToSegment]=assignSegment(Geometry,Sample_x(end),Sample_y(end));
    OnSegmentPoint_x=Geometry.PlateBoundary.x(SegmentNumber)+Geometry.PlateBoundaryLength_x(SegmentNumber)*...
        (AlongSegmentCoordinate-sum(Geometry.PlateBoundaryLength(1:SegmentNumber-1)))/Geometry.PlateBoundaryLength(SegmentNumber);
    OnSegmentPoint_y=Geometry.PlateBoundary.y(SegmentNumber)+Geometry.PlateBoundaryLength_y(SegmentNumber)*...
        (AlongSegmentCoordinate-sum(Geometry.PlateBoundaryLength(1:SegmentNumber-1)))/Geometry.PlateBoundaryLength(SegmentNumber);
    OnSegmentPoint_z=interp2(Grid.x,Grid.y,Lid.Depth,OnSegmentPoint_x,OnSegmentPoint_y);

    if  (DistanceToSegment<ExtendToSegmentDistance);
        Sample_Distance=[Sample_Distance;Sample_Distance(end)+DistanceToSegment];
        Sample_x=[Sample_x;OnSegmentPoint_x];
        Sample_y=[Sample_y;OnSegmentPoint_y];
        Sample_z=[Sample_z;OnSegmentPoint_z];
    end
end

% % resample lines to same size
if numel(Sample_Distance)>1;
    NewSample_Distance=linspace(Sample_Distance(1),Sample_Distance(end),Res.nMeltSwath);
    NewSample_x=interp1(Sample_Distance,Sample_x,NewSample_Distance);
    NewSample_y=interp1(Sample_Distance,Sample_y,NewSample_Distance);
    NewSample_z=interp1(Sample_Distance,Sample_z,NewSample_Distance);
    Sample_Slope=[interp2(Grid.x,Grid.y,Lid.Slope,NewSample_x,NewSample_y)];
else % duplicate point if line contains only one point
    NewSample_Distance=repmat(Sample_Distance(1),1,Res.nMeltSwath);
    NewSample_x=repmat(Sample_x(1),1,Res.nMeltSwath);
    NewSample_y=repmat(Sample_y(1),1,Res.nMeltSwath);
    NewSample_z=repmat(Sample_y(1),1,Res.nMeltSwath);
    Sample_Slope=[interp2(Grid.x,Grid.y,Lid.Slope,NewSample_x,NewSample_y)];
end

LineUpdated=[NewSample_x;NewSample_y;NewSample_z;NewSample_Distance;Sample_Slope]';
OnSegmentInfo=[OnSegmentPoint_x,OnSegmentPoint_y,SegmentNumber,AlongSegmentCoordinate,DistanceToSegment];

return