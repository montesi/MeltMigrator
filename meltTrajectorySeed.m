function Seed=meltTrajectorySeed(Geometry,Grid,Res,Lid,T_MeltSeed,SpreadingDirection,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltTrajectorySeed.m
% Determine seed for melt lines by following contour of ExcessT
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |.PlateBoundaryType : Plate boundary type (1:ridge, 2:transform, 3:oblique segment
%       |......
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates [km]
%   Res
%       |.MeltSeed          : Distance between adjacent melt seeds [km]
%       |......
%   Lid
%       |.Depth             : Depth of lid [km]
%       |.ExcessT           : Excess temperature above solidus [degC]
%       |......
%   T_MeltSeed              : Desired countour of ExcessT for seed generation [degC]
%   SpreadingDirection      : Direction of seafloor spreading [deg]
%   indFigure               : Figure index
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Seed(ind)
%       |.x, .y             : Coordinates of melt line seeds
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   RidgeIndex              : Index of ridge segment along plate boundaries
%   TargetRidgeIndex        : Index of target (first and last) ridge segment
%   TargetPoint             : Coordinates for centers of target ridge segment [km]
%   Section
%       |.x, .y             : Coordinates for sections parallel to spreading direction and going through centers of ridge segment [km]
%       |.ExcessT           : Excess temperature above solidus along section [degC]
%   MaxT                    : Maximum ExcessT along section [degC]
%   iStart                  : Index for starting point for first melt seeds searching
%   SearchRange             : Search range for first melt seeds
%   InitialSeed
%       |_x, _y             : Coordinates for first melt seeds [km]
%   SeedDis                 : Distance along melt seed searching line [km]
%   SeedCoordinate1         : Coordinates for melt seeds from one searching direction [km]
%   SeedCoordinate2         : Coordinates for melt seeds from the other searching direction [km]
%   SeedCoordinate          : Coordinates for melt seeds [km]
%   SeedCoordinateDiff      : Distance between adjacent melt seeds [km]
%   iSeedDis                : Index for non-overlapping seeds
%   Seed_z                  : Depth of melt seeds [km]
%   InBoundaryIndex         : Index for point of interest within model boundary
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   crossToGrad
%   inBounds
%   lineSegmentIntersect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section Location

figure(indFigure); clf; hold on;
contourf(Grid.x,Grid.y,Lid.ExcessT,[-40:20:200]); colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');

% % find sections parallel to spreading direction and going through center of ridge axes 
RidgeIndex=find(Geometry.PlateBoundaryType==1); % locate ridge segments along plate boundaries
[SortRidgeLength,SortRidgeLengthIndex]=sort(Geometry.PlateBoundary.y(RidgeIndex)); % sort ridge segment based on y coordinate
% TargetRidgeIndex=RidgeIndex(SortRidgeLengthIndex(1,:)); % index of ridge segments
TargetRidgeIndex=RidgeIndex(SortRidgeLengthIndex(1,[1,end])); % index of ridge segments

for i=1:length(TargetRidgeIndex);
    TargetPoint(i,1)=mean([Geometry.PlateBoundary.x(TargetRidgeIndex(i)),Geometry.PlateBoundary.x(TargetRidgeIndex(i)+1)]); % center of ridge segment
    TargetPoint(i,2)=mean([Geometry.PlateBoundary.y(TargetRidgeIndex(i)),Geometry.PlateBoundary.y(TargetRidgeIndex(i)+1)]); % center of ridge segment
    
    % % locate a section across center of ridge segment, parallel to spreading direction
    if SpreadingDirection==0; % if spreading is along x-direction
        Section(i).x=[Geometry.ModelBoundary.x(1):Res.MeltSeed:Geometry.ModelBoundary.x(2)];
        Section(i).y=linspace(TargetPoint(i,2),TargetPoint(i,2),length(Section(i).x));
    elseif SpreadingDirection==90; % if spreading is along y-direction
        Section(i).y=[Geometry.ModelBoundary.y(1):Res.MeltSeed:Geometry.ModelBoundary.y(2)];
        Section(i).x=linspace(TargetPoint(i,1),TargetPoint(i,1),length(Section(i).y));
    else % if spreading is along other direction
        sampleLine=@(x) tand(SpreadingDirection).*x+TargetPoint(i,2)-tand(SpreadingDirection).*TargetPoint(i,1);
        Section(i).x=[Geometry.ModelBoundary.x(1):Res.MeltSeed:Geometry.ModelBoundary.x(2)];
        Section(i).y=sampleLine(Section(i).x);
        InBoundaryIndex=find((Section(i).y>=Geometry.ModelBoundary.y(1))&(Section(i).y<=Geometry.ModelBoundary.y(2)));
        Section(i).x=Section(i).x(InBoundaryIndex);
        Section(i).y=Section(i).y(InBoundaryIndex);
    end
    plot(Section(i).x,Section(i).y,'k','linewidth',2);
    Section(i).ExcessT=interp2(Grid.x,Grid.y,Lid.ExcessT,Section(i).x,Section(i).y); % get ExcessT data along section
end

%% Seeds Generation

options=odeset('reltol',1e-4); % odeset('reltol',1e-8), optimset('display','off')
% % along sections, find points with assigned ExcessT value and search for contour lines
for iSeed=1:length(Section); % two original seed-searching sections
    [MaxT,iStart]=max(Section(iSeed).ExcessT); % identify where melting is max along bottom boundary
    switch rem(iSeed,2); % each search starts from one side of ridge axis
        case 0
            SearchRange=[1:iStart];
        case 1
            SearchRange=[iStart:length(Section(iSeed).x)];
    end
    
    % % find first seeds along sections
    if mean(diff(Section(iSeed).y))==0; % if section is along x-direction
        InitialSeed_x=fzero(@(x) interp1(Section(iSeed).x(SearchRange),Section(iSeed).ExcessT(SearchRange),x)-T_MeltSeed,...
            [Section(iSeed).x(SearchRange(1)),Section(iSeed).x(SearchRange(end))],options);
        % % alternatively, use fsolve if optimization toolbox is available
%         InitialSeed_x=fsolve(@(x) interp1(Section(iSeed).x(SearchRange),Section(iSeed).ExcessT(SearchRange),x)-T_MeltSeed,...
%             [Section(iSeed).x(iStart)],options);
        InitialSeed_y=interp1(Section(iSeed).x,Section(iSeed).y,InitialSeed_x);
    else % if section is not along x-direction
        InitialSeed_y=fzero(@(y) interp1(Section(iSeed).y(SearchRange),Section(iSeed).ExcessT(SearchRange),y)-T_MeltSeed,...
            [Section(iSeed).y(SearchRange(1)),Section(iSeed).y(SearchRange(end))],options);
        % % alternatively, use fsolve if optimization toolbox is available
%         InitialSeed_y=fsolve(@(y) interp1(Section(iSeed).y(SearchRange),Section(iSeed).ExcessT(SearchRange),y)-T_MeltSeed,...
%             [Section(iSeed).y(iStart)],options);
        InitialSeed_x=interp1(Section(iSeed).y,Section(iSeed).x,InitialSeed_y);
    end
    scatter(InitialSeed_x,InitialSeed_y);
    
    % % starting from first seeds, find rest seeds
    [SeedDis,SeedCoordinate1]=ode113(@(l,P) crossToGrad(P,Grid,Geometry.ModelBoundary,Lid,1),[0:Res.MeltSeed:1000],...
        [InitialSeed_x,InitialSeed_y],options);
    [SeedDis,SeedCoordinate2]=ode113(@(l,P) crossToGrad(P,Grid,Geometry.ModelBoundary,Lid,2),[0:Res.MeltSeed:1000],...
        [InitialSeed_x,InitialSeed_y],options);
    clear SeedDis;
    
    % % find the seed contours located within the model domain, these seed contours may have looping problems, 
    % while seed contours cut by the model edges don't 
    EdgeFinder=[];
    for iEdgeFinder=1:size(SeedCoordinate1,1);
        Seed_z(iEdgeFinder)=interp2(Grid.x,Grid.y,Lid.Depth,SeedCoordinate1(iEdgeFinder,1),SeedCoordinate1(iEdgeFinder,2));
        EdgeFinder=[EdgeFinder,inBounds(SeedCoordinate1(iEdgeFinder,1),SeedCoordinate1(iEdgeFinder,2),Seed_z(iEdgeFinder),Geometry.ModelBoundary)];
    end
    for iEdgeFinder=1:size(SeedCoordinate2,1);
        Seed_z(iEdgeFinder)=interp2(Grid.x,Grid.y,Lid.Depth,SeedCoordinate2(iEdgeFinder,1),SeedCoordinate2(iEdgeFinder,2));
        EdgeFinder=[EdgeFinder,inBounds(SeedCoordinate2(iEdgeFinder,1),SeedCoordinate2(iEdgeFinder,2),Seed_z(iEdgeFinder),Geometry.ModelBoundary)];
    end
    % % clear the looping portion of seed contours
    if min(EdgeFinder)==1;
        SeedSegment1=[SeedCoordinate1(2:end-1,1),SeedCoordinate1(2:end-1,2),SeedCoordinate1(3:end,1),SeedCoordinate1(3:end,2)];
        SeedSegment2=[SeedCoordinate2(2:end-1,1),SeedCoordinate2(2:end-1,2),SeedCoordinate2(3:end,1),SeedCoordinate2(3:end,2)];
        SectionSegment=[Section(iSeed).x(1:end-1)',Section(iSeed).y(1:end-1)',Section(iSeed).x(2:end)',Section(iSeed).y(2:end)'];
        Intersection1=lineSegmentIntersect(SeedSegment1,SectionSegment);
        Intersection2=lineSegmentIntersect(SeedSegment2,SectionSegment);
        [IntInd1_Row,IntInd1_Col]=find(Intersection1.intAdjacencyMatrix==1);
        [IntInd2_Row,IntInd2_Col]=find(Intersection2.intAdjacencyMatrix==1);
        IntInd1=min(IntInd1_Row);
        [IntInd2,iIntInd2]=min(IntInd2_Row);
        SeedCoordinate1=SeedCoordinate1(1:IntInd1+1,:);
        SeedCoordinate2=[SeedCoordinate2(1:IntInd2+1,:);...
            Intersection2.intMatrixX(IntInd2_Row(iIntInd2),IntInd2_Col(iIntInd2)),Intersection2.intMatrixY(IntInd2_Row(iIntInd2),IntInd2_Col(iIntInd2))];
    end
    SeedCoordinate=[flipud(SeedCoordinate1);SeedCoordinate2];
    
    % % remove overlapping seeds
    SeedCoordinateDiff=diff(SeedCoordinate);
    SeedDis=sqrt(SeedCoordinateDiff(:,1).^2+SeedCoordinateDiff(:,2).^2);
    iSeedDis=find(SeedDis~=0); 
    SeedCoordinate=SeedCoordinate(iSeedDis,:);
    
    % % store seed information
    Seed_x=[];Seed_y=[];
    for iIn=1:size(SeedCoordinate,1);
        Seed_z(iIn)=interp2(Grid.x,Grid.y,Lid.Depth,SeedCoordinate(iIn,1),SeedCoordinate(iIn,2));
        InBoundaryIndex=inBounds(SeedCoordinate(iIn,1),SeedCoordinate(iIn,2),Seed_z(iIn),Geometry.ModelBoundary);
        if InBoundaryIndex;
            Seed_x=[Seed_x;SeedCoordinate(iIn,1)];
            Seed_y=[Seed_y;SeedCoordinate(iIn,2)];
        end
    end
    Seed(iSeed).x=Seed_x;
    Seed(iSeed).y=Seed_y;
    
    disp(sprintf('>>>    Melt seed contour from (%g,%g) to (%g,%g)',SeedCoordinate(1,:),SeedCoordinate(end,:)));
end

%% Seeds Plot

for iSeed=1:numel(Seed);
    plot(Seed(iSeed).x,Seed(iSeed).y,'w.');
end
axis equal; axis tight; box on;
title('ExcessT and Melt Line Seeds');

return