function [Line,LineDis]=lineCleanUp(Line,LineDis,Grid,Lid,Saddle,SeparationWidth,Switch_lineCleanUp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lineCleanUp.m
% Clean up melt trajectories, remove portion of line within a distance of permeability barrier crest to creat better-shaped melt swath
% LIMITATION: this script only works in simple ridge geometry setting with spreading in x-direction
% Hailong Bai
% June 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Line                    : Original melt line coordinates [km]
%   LineDis                 : Distance along original melt line [km]
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates [km]
%   Lid
%       |.Depth             : Depth of lid [km]
%       |......
%   Saddle(ind)
%       |.Up                : Upgoing saddle line coordinates [km]
%       |......
%   SeparationWidth         : Half-width of 'separation zone (no melt line zone)' [km]
%   Switch_LineCleanUp      : Switch to turn on/off LineCleanUp function
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Line                    : New melt line coordinates after trim [km]
%   LineDis                 : New distance along melt line after trim [km]
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   Separation              : Separation line coordinates [km]
%   iSaddle                 : Saddle index
%   SeparationSaddle        : Coordinates for saddle line portion of separation line [km]
%   Crest_z                 : Minimum depth of lid along each row [km]
%   iCrest                  : Index for lid crest along each row
%   Crest_xy                : x and y coordinate of lid crest [km]
%   iCrestSelected          : Index for crest not covered by upgoing saddle line
%   SeparationCrest         : Coordinates for crest portion of separation line [km]
%   SeparationLeft          : Coordinates for left boundary of 'separation zone' [km]
%   SeparationRight         : Coordinates for right boundary of 'separation zone' [km]
%   LineSegment             : Segment of melt line constructed for intersection determination [km]
%   Intersection
%       | Left, RIght       : Intersection information for melt line and separation line (left/right)
%   indIntersection
%       | Left, Right       : Intersection indicator
%   indOverlap
%       | Left, Right       : Overlap indicator
%   Row
%       |_Intersection
%           | Left Right    : Row index for intersection points
%       |_Overlap
%           | Left Right    : Row index for overlap points
%       |_Indicator
%           | Left Right    : Row index for intersection or overlap points
%   Column
%       |_Intersection      
%           | Left Right    : Column index for intersection points
%       |_Overlap
%           | Left Right    : Column index for overlap points
%       |_Indicator
%           | Left Right    : Column index for intersection or overlap points
%   iEnd
%       | Left, Right       : Index for the first intersection point
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   lineSegmentIntersect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if Switch_lineCleanUp;
    %% Separation Line Assemblage
    % % separation line consists of upgoing saddle lines and crest line of lid (when saddle lines are not present)
    Separation=[]; % initialize separation line
    for iSaddle=1:numel(Saddle); % assemble separation line using each saddle
        % % sort the upgoing saddle line descendingly based on y coordinate, (assemble separation line along y-direction)
        SeparationSaddle=sortrows(Saddle(iSaddle).Up,-2);
        if isempty(Separation); % assemble the first portion of separation line
            [Crest_z,iCrest]=min(Lid.Depth,[],2); % find minimum lid depth along each row
            Crest_xy=flipud([Grid.x(1,iCrest)',Grid.y(:,1)]); % find location of lid crest
            iCrestSelected=find(Crest_xy(:,2)>max(Saddle(iSaddle).Up(:,2))); % find the part of crest not covered by upgoing saddle line
            SeparationCrest=Crest_xy(iCrestSelected,:); % portion of crest selected to use as separation
            Separation=[SeparationCrest(1:end-1,:),SeparationCrest(2:end,:);SeparationSaddle(1:end-1,:),SeparationSaddle(2:end,:)]; % form separation line
        else
            if ~mod(iSaddle,2);
                % % between two saddle structures, crest portion is necessary
                if iSaddle~=numel(Saddle); % if not the last saddle, use crest to fill space between saddle lines
                    iCrestSelected=find((Crest_xy(:,2)<min(Saddle(iSaddle).Up(:,2)))&(Crest_xy(:,2)>max(Saddle(iSaddle+1).Up(:,2))));
                else % if at last saddle, use crest to fill space between saddle line and model boundary
                    iCrestSelected=find((Crest_xy(:,2)<min(Saddle(iSaddle).Up(:,2)))&(Crest_xy(:,2)>0));
                end
                SeparationCrest=Crest_xy(iCrestSelected,:);
                Separation=[Separation;Separation(end,end-1:end),SeparationSaddle(1,:);SeparationSaddle(1:end-1,:),SeparationSaddle(2:end,:);...
                        SeparationCrest(1:end-1,:),SeparationCrest(2:end,:)];
            else
                % % within the same saddle structure, saddle lines are close enough, no crest is needed
                Separation=[Separation;Separation(end,end-1:end),SeparationSaddle(1,:);SeparationSaddle(1:end-1,:),SeparationSaddle(2:end,:)];
            end
        end
    end

    SeparationLeft=[Separation(:,1)-SeparationWidth,Separation(:,2),Separation(:,3)-SeparationWidth,Separation(:,4)]; % left boundary of clean-up zone
    SeparationRight=[Separation(:,1)+SeparationWidth,Separation(:,2),Separation(:,3)+SeparationWidth,Separation(:,4)]; % right boundary of clean-up zone

    %% Line Clean-Up

    LineSegment=[Line(1:end-1,:),Line(2:end,:)]; % line segment for intersection determination

    % % determine if line segment intersects with separation lines
    IntersectionLeft=lineSegmentIntersect(LineSegment,SeparationLeft);
    indIntersectionLeft=max(max(IntersectionLeft.intAdjacencyMatrix));
    indOverlapLeft=max(max(IntersectionLeft.coincAdjacencyMatrix));
    IntersectionRight=lineSegmentIntersect(LineSegment,SeparationRight); 
    indIntersectionRight=max(max(IntersectionRight.intAdjacencyMatrix));
    indOverlapRight=max(max(IntersectionRight.coincAdjacencyMatrix));

    % % if intersect
    if indIntersectionLeft||indIntersectionRight||indOverlapLeft||indOverlapRight;
        [Row_IntersectionLeft,Column_IntersectionLeft]=find(IntersectionLeft.intAdjacencyMatrix==1);
        [Row_IntersectionRight,Column_IntersectionRight]=find(IntersectionRight.intAdjacencyMatrix==1);
        [Row_OverlapLeft,Column_OverlapLeft]=find(IntersectionLeft.coincAdjacencyMatrix==1);
        [Row_OverlapRight,Column_OverlapRight]=find(IntersectionRight.coincAdjacencyMatrix==1);
        % % find the first intersection point
        Row_IndicatorLeft=[Row_IntersectionLeft,Row_OverlapLeft];
        Row_IndicatorRight=[Row_IntersectionRight,Row_OverlapRight];
        Column_IndicatorLeft=[Column_IntersectionLeft,Column_OverlapLeft];
        Column_IndicatorRight=[Column_IntersectionRight,Column_OverlapRight];
        [iEndLeft,iiEndLeft]=min(Row_IndicatorLeft);
        [iEndRight,iiEndRight]=min(Row_IndicatorRight);
        if isempty(iEndLeft);
            iEndLeft=iEndRight+1;
        end
        if isempty(iEndRight);
            iEndRight=iEndLeft+1;
        end
        % % cut the line at intersection point
        % % it's unlikely the line intersects with both separation boundaries
        if iEndLeft<iEndRight;
            Line=[Line(1:iEndLeft,:);IntersectionLeft.intMatrixX(iEndLeft,Column_IndicatorLeft(iiEndLeft)),IntersectionLeft.intMatrixY(iEndLeft,Column_IndicatorLeft(iiEndLeft))];
            LineDis=[LineDis(1:iEndLeft);LineDis(iEndLeft)+IntersectionLeft.intNormalizedDistance1To2(iEndLeft,Column_IndicatorLeft(iiEndLeft))*(LineDis(iEndLeft)/abs(LineDis(iEndLeft)))];
        else
            Line=[Line(1:iEndRight,:);IntersectionRight.intMatrixX(iEndRight,Column_IndicatorRight(iiEndRight)),IntersectionRight.intMatrixY(iEndRight,Column_IndicatorRight(iiEndRight))];
            LineDis=[LineDis(1:iEndRight);LineDis(iEndRight)+IntersectionRight.intNormalizedDistance1To2(iEndRight,Column_IndicatorRight(iiEndRight))*(LineDis(iEndRight)/abs(LineDis(iEndRight)))];
        end
    end
end
return