function [Shoulder,Crust,Cryptic] = meltExtraction(Geometry,Res,Shoulder,ExtractionWidth,ExtractionDepth,ExtractionSlope,SpreadingRate_Half)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltExtraction.m
% Determine extraction state of each patch of melt and calculate crustal thickness along plate boundaries
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.PlateBoundaryLength : Length of plate boundary segments [km]
%       |......
%   Res
%       |.Axis              : Plate boundary sampling resolution [km]
%   Shoulder(ind)
%       |.Swath(ind)
%           |.Polygon(ind)
%               |.x,.y      : Coordinates of polygon outline [km]
%           |.PolygonArea   : Area of each polygon [km^2]
%           |.PolygonCenter_z : Depth of polygon centers [km]
%           |.PolygonSlope  : Slope at polygon centers
%           |.MaxMeltFraction : Maximum melt fraction at sampling columns
%           |.AverageMeltFraction : Average melt fraction at sampling columns 
%           |.MeltFlux      : Melt flux at sampling columns [m/s]
%           |.MeltingDepth  : Melting depth at sampling columns [km]
%           |......
%       |......
%   ExtractionWidth         : Melt extraction width [km]
%   ExtractionDepth         : Melt extraction depth [km]
%   ExtractionSlope         : Melt transport critical slope
%   Spreading_Half          : Half-spreading rate [m/s]
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Shoulder
%       |.Distribution
%           |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km] 
%           |.Thickness     : Distribution of erupted crust along plate boundary from each shoulder [km]
%           |.AverageMeltFraction : Average melt fraction along plate boundary from each shoulder
%           |.MaxMeltFraction : Maximum melt fraction along plate boundary from each shoulder
%           |.MeltingDepth  : Melting depth along plate boundary from each shoulder [km]
%           |.CrypticThickness : Distribution of cryptic crust along plate boundary from each shoulder [km]
%           |.CrypticAverageMeltFraction : Average cryptic melt fraction along plate boundary from each shoulder
%           |.CrypticMaxMeltFraction : Maximum cryptic melt fraction along plate boundary from each shoulder
%           |.CrypticMeltingDepth : Cryptic melting depth along plate boundary from each shoulder [km]
%   Crust
%       |.Thickness         : Distribution of erupted crust along plate boundary [km] 
%       |.AverageMeltFraction : Average melt fraction along plate boundary
%       |.MaxMeltFraction   : Maximum melt fraction along plate boundary 
%       |.MeltingDepth      : Melting depth along plate boundary [km] 
%       |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%   Cryptic
%       |.Thickness         : Distribution of cryptic crust along plate boundary [km] 
%       |.AverageMeltFraction : Average cryptic melt fraction along plate boundary
%       |.MaxMeltFraction   : Maximum cryptic melt fraction along plate boundary
%       |.MeltingDepth      : Cryptic melting depth along plate boundary [km]
%       |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   nPlateBoundarySampling  : Sampling size along plate boundary
%   nPolygon                : Number of polygons within one swath
%   iShoulder               : Index of shoulder
%   iSwath                  : Index of swath
%   iPolygon                : Index of polygon
%   iExtraction             : Indicator of extraction
%   iStart                  : Index of starting point of extraction zone
%   iEnd                    : Index of end point of extraction zone
%   LocalSwath              : Current melt swath, short for Shoulder(iShoulder).Swath(iSwath) 
%   AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%   AverageMeltFraction     : Average melt fraction along plate boundary from each shoulder
%   MaxMeltFraction         : Maximum melt fraction along plate boundary from each shoulder
%   MeltingDepth            : Melting depth along plate boundary from each shoulder [km]
%   CrypticAverageMeltFraction : Average cryptic melt fraction along plate boundary from each shoulder
%   CrypticMaxMeltFraction  : Maximum cryptic melt fraction along plate boundary from each shoulder
%   CrypticMeltingDepth     : Cryptic melting depth along plate boundary from each shoulder [km]
%   LocalMeltFlux           : Accumulated melt flux during melt transport along swath [m/s] 
%   LocalAverageMeltFraction : Accumulated average melt fraction during melt transport along swath
%   LocalMaxMeltFraction    : Accumulated maximum melt fraction during melt transport along swath
%   LocalMeltingDepth       : Accumulated melting depth during melt transport along swath [km]
%   NewMeltFlux             : Melt flux from current polygon [m/s]
%   NewAverageMeltFraction  : Average melt fraction from current polygon
%   NewMaxMeltFraction      : Maximum melt fraction from current polygon
%   NewMeltingDepth         : Melting depth from current polygon [km]
%   LocalEruptedAverageMeltFraction : Average erupted melt fraction from each swath 
%   LocalEruptedMaxMeltFraction : Maximum erupted melt fraction from each swath
%   LocalEruptedMeltgingDepth   : Erupted melting depth from each swath [km]
%   LocalCrypticAverageMeltFraction : Average cryptic melt fraction from each swath 
%   LocalCrypticMaxMeltFraction : Maximum crypti melt fraction from each swath
%   LocalCrypticMeltingDepth    : Cryptic melting depth from each swath [km]
%   SegmentNumber           : Index of plate boundary segment polygon corner points associated with
%   AlongSegmentCoordinate  : Along plate boundary coordinates of polygon corner point projection [km]
%   DistanceToSegment       : Distance between polygon corner points and plate boundary [km]
%   ExtractionCoordinate    : Along plate boundary coordinates of extraction zone [km] 
%   LengthOnSegment         : Length of extraction zone [km]
%   Thickness               : Crustal thickness from current swath [km]
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   extractionDetermination
%   assignSegment
%   minDistBetweenTwoPolygons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% % calculate crustal thickness along 1D coordinate along plate boundary
% % initial 1D coordinate system
warning off all
nPlateBoundarySampling=ceil(sum(Geometry.PlateBoundaryLength)/Res.Axis);
AlongPlateBoundaryCoordinate=linspace(0,sum(Geometry.PlateBoundaryLength),nPlateBoundarySampling);
% % initial crust and cryptic crust structure
Crust.Thickness=AlongPlateBoundaryCoordinate*0;
Crust.AverageMeltFraction=AlongPlateBoundaryCoordinate*0;
Crust.MaxMeltFraction=AlongPlateBoundaryCoordinate*0;
Crust.MeltingDepth=AlongPlateBoundaryCoordinate*0;
Cryptic=Crust;

for iShoulder=1:numel(Shoulder)
    clear Distribution;
    % % 'Distribution' is crustal thickness profile structure for individual shoulder 
    Distribution.AlongPlateBoundaryCoordinate=AlongPlateBoundaryCoordinate;
    Distribution.Thickness=AlongPlateBoundaryCoordinate*0;
    Distribution.CrypticThickness=AlongPlateBoundaryCoordinate*0;
    
    % % initial local data structure
    AverageMeltFraction=AlongPlateBoundaryCoordinate*0;    
    MaxMeltFraction=AlongPlateBoundaryCoordinate*0;
    MeltingDepth=AlongPlateBoundaryCoordinate*0;
    CrypticAverageMeltFraction=AlongPlateBoundaryCoordinate*0;
    CrypticMaxMeltFraction=AlongPlateBoundaryCoordinate*0;
    CrypticMeltingDepth=AlongPlateBoundaryCoordinate*0;
    
    for iSwath=1:numel(Shoulder(iShoulder).Swath)
        LocalSwath=Shoulder(iShoulder).Swath(iSwath);
        nPolygons=numel(LocalSwath.PolygonCenter_x); % number of polygons
        % % initialize accumulated melts along swath
        LocalMeltFlux=0;
        LocalAverageMeltFraction=0; 
        LocalMaxMeltFraction=0; 
        LocalMeltingDepth=0;

        for iPolygon=1:nPolygons;
            % % calculate incoming melts from each polygon
            NewMeltFlux=LocalSwath.MeltFlux(iPolygon)*LocalSwath.PolygonArea(iPolygon); % melt flux
            NewAverageMeltFraction=LocalSwath.AverageMeltFraction(iPolygon)*NewMeltFlux; % average melt fraction
            NewMaxMeltFraction=LocalSwath.MaxMeltFraction(iPolygon)*NewMeltFlux; % maximum melt fraction
            NewMeltingDepth=LocalSwath.MeltingDepth(iPolygon)*NewMeltFlux; % depth of melting 

            % % accumulate melt along each swath, pass melt on to next polygon
            LocalMeltFlux=LocalMeltFlux+NewMeltFlux;
            LocalAverageMeltFraction=LocalAverageMeltFraction+NewAverageMeltFraction; 
            LocalMaxMeltFraction=LocalMaxMeltFraction+NewMaxMeltFraction;
            LocalMeltingDepth=LocalMeltingDepth+NewMeltingDepth;
            % % determine if current patch of melt is extracted
            iExtraction(iPolygon)=extractionDetermination(Geometry,LocalSwath,iPolygon,ExtractionWidth,ExtractionDepth,ExtractionSlope);
            
            % % when melt stops moving, find out where it ends
            if iExtraction(iPolygon)~=0;
                % % find out along segment coordinates of corner points of current polygon 
                [SegmentNumber1,AlongSegmentCoordinate1,DistanceToSegment1]=...
                    assignSegment(Geometry,LocalSwath.Polygon(iPolygon).x(1),LocalSwath.Polygon(iPolygon).y(1)); 
                [SegmentNumber2,AlongSegmentCoordinate2,DistanceToSegment2]=...
                    assignSegment(Geometry,LocalSwath.Polygon(iPolygon).x(2),LocalSwath.Polygon(iPolygon).y(2));
                [SegmentNumber3,AlongSegmentCoordinate3,DistanceToSegment3]=...
                    assignSegment(Geometry,LocalSwath.Polygon(iPolygon).x(3),LocalSwath.Polygon(iPolygon).y(3));
                ExtractionCoordinate=[AlongSegmentCoordinate1;AlongSegmentCoordinate2;AlongSegmentCoordinate3];
                % % find out plate boundary coordinates of melt patch (extraction zone coordinates)
                if min(ExtractionCoordinate)<=min(AlongPlateBoundaryCoordinate); % if assigned location is out of lower limit 
                    iStart=1;
                else
                    % % find out starting point of extraction zone 
                    iStart=max(find(AlongPlateBoundaryCoordinate<min(ExtractionCoordinate)));
                end
                if isnan(iStart);
                    iStart=1;
                end
                if max(ExtractionCoordinate)>=max(AlongPlateBoundaryCoordinate); % if assigned location is out of upper limit
                    iEnd=numel(AlongPlateBoundaryCoordinate);
                else
                    % % find out end point of extraction zone
                    iEnd=min(find(AlongPlateBoundaryCoordinate>max(ExtractionCoordinate)));
                end
                if isnan(iEnd)
                    iEnd=numel(AlongPlateBoundaryCoordinate);
                end
                LengthOnSegment=abs(AlongPlateBoundaryCoordinate(iEnd)-AlongPlateBoundaryCoordinate(iStart)); % length of melt-receiving segment
                Thickness=LocalMeltFlux/(2*SpreadingRate_Half*LengthOnSegment); % crustal thickness accumulated on segment
                
                % % assign melt to erupted/cryptic crust based on its extraction state 
                if iExtraction(iPolygon)==1; % extract as actual crust;
                    if LocalMeltFlux==0
                        LocalEruptedAverageMeltFraction=0;
                        LocalEruptedMaxMeltFraction=0;
                        LocalEruptedMeltingDepth=0;
                    else
                        LocalEruptedAverageMeltFraction=LocalAverageMeltFraction/LocalMeltFlux;
                        LocalEruptedMaxMeltFraction=LocalMaxMeltFraction/LocalMeltFlux;
                        LocalEruptedMeltingDepth=LocalMeltingDepth/LocalMeltFlux;
                    end
                    Distribution.Thickness(iStart:iEnd)=Distribution.Thickness(iStart:iEnd)+Thickness;
                    AverageMeltFraction(iStart:iEnd)=AverageMeltFraction(iStart:iEnd)+LocalEruptedAverageMeltFraction.*Thickness;
                    MaxMeltFraction(iStart:iEnd)=MaxMeltFraction(iStart:iEnd)+LocalEruptedMaxMeltFraction.*Thickness;
                    MeltingDepth(iStart:iEnd)=MeltingDepth(iStart:iEnd)+LocalEruptedMeltingDepth.*Thickness;
                elseif iExtraction(iPolygon)==2; % extract as cryptic crust
                    if LocalMeltFlux==0
                        LocalCrypticAverageMeltFraction=0;
                        LocalCrypticMaxMeltFraction=0;
                        LocalCrypticMeltingDepth=0;
                    else
                        LocalCrypticAverageMeltFraction=LocalAverageMeltFraction/LocalMeltFlux;
                        LocalCrypticMaxMeltFraction=LocalMaxMeltFraction/LocalMeltFlux;
                        LocalCrypticMeltingDepth=LocalMeltingDepth/LocalMeltFlux;
                    end
                    Distribution.CrypticThickness(iStart:iEnd)=Distribution.CrypticThickness(iStart:iEnd)+Thickness;
                    CrypticAverageMeltFraction(iStart:iEnd)=CrypticAverageMeltFraction(iStart:iEnd)+LocalCrypticAverageMeltFraction.*Thickness;
                    CrypticMaxMeltFraction(iStart:iEnd)=CrypticMaxMeltFraction(iStart:iEnd)+LocalCrypticMaxMeltFraction.*Thickness;
                    CrypticMeltingDepth(iStart:iEnd)=CrypticMeltingDepth(iStart:iEnd)+LocalCrypticMeltingDepth.*Thickness;
                end
                % % reset melt accumulation after extraction/refertilization
                LocalMeltFlux=0;
                LocalAverageMeltFraction=0;
                LocalMaxMeltFraction=0;
                LocalMeltingDepth=0;
            end
            % % store extraction information
            Shoulder(iShoulder).Swath(iSwath).Extraction=iExtraction;
        end
    end
    
    % % sum melt extraction/refertilization information from each swath within each shoulder 
    AverageMeltFraction=AverageMeltFraction./Distribution.Thickness;
    AverageMeltFraction(find(isnan(AverageMeltFraction)))=0;
    MaxMeltFraction=MaxMeltFraction./Distribution.Thickness;
    MaxMeltFraction(find(isnan(MaxMeltFraction)))=0;
    MeltingDepth=MeltingDepth./Distribution.Thickness;
    MeltingDepth(find(isnan(MeltingDepth)))=0;
    Distribution.AverageMeltFraction=AverageMeltFraction;
    Distribution.MaxMeltFraction=MaxMeltFraction;
    Distribution.MeltingDepth=MeltingDepth;
                    
    CrypticAverageMeltFraction=CrypticAverageMeltFraction./Distribution.CrypticThickness;
    CrypticAverageMeltFraction(find(isnan(CrypticAverageMeltFraction)))=0;
    CrypticMaxMeltFraction=CrypticMaxMeltFraction./Distribution.CrypticThickness;
    CrypticMaxMeltFraction(find(isnan(CrypticMaxMeltFraction)))=0;
    CrypticMeltingDepth=CrypticMeltingDepth./Distribution.CrypticThickness;
    CrypticMeltingDepth(find(isnan(CrypticMeltingDepth)))=0;
    Distribution.CrypticAverageMeltFraction=CrypticAverageMeltFraction;
    Distribution.CrypticMaxMeltFraction=CrypticMaxMeltFraction;
    Distribution.CrypticMeltingDepth=CrypticMeltingDepth;
    
    Shoulder(iShoulder).Distribution=Distribution;
    
    Crust.Thickness=Crust.Thickness+Distribution.Thickness;
    Crust.AverageMeltFraction=Crust.AverageMeltFraction+Distribution.AverageMeltFraction.*Distribution.Thickness;
    Crust.MaxMeltFraction=Crust.MaxMeltFraction+Distribution.MaxMeltFraction.*Distribution.Thickness;
    Crust.MeltingDepth=Crust.MeltingDepth+Distribution.MeltingDepth.*Distribution.Thickness;
    
    Cryptic.Thickness=Cryptic.Thickness+Distribution.CrypticThickness;
    Cryptic.AverageMeltFraction=Cryptic.AverageMeltFraction+Distribution.CrypticAverageMeltFraction.*Distribution.CrypticThickness;
    Cryptic.MaxMeltFraction=Cryptic.MaxMeltFraction+Distribution.CrypticMaxMeltFraction.*Distribution.CrypticThickness;
    Cryptic.MeltingDepth=Cryptic.MeltingDepth+Distribution.CrypticMeltingDepth.*Distribution.CrypticThickness;
end

% % sum melt extraction/refertilization information from each shoulder 
Crust.Thickness=abs(Crust.Thickness);
Cryptic.Thickness=abs(Cryptic.Thickness);
Crust.AverageMeltFraction=abs(Crust.AverageMeltFraction./Crust.Thickness);
Crust.AverageMeltFraction(find(isnan(Crust.AverageMeltFraction)))=0;
Crust.MaxMeltFraction=abs(Crust.MaxMeltFraction./Crust.Thickness);
Crust.MaxMeltFraction(find(isnan(Crust.MaxMeltFraction)))=0;
Crust.MeltingDepth=abs(Crust.MeltingDepth./Crust.Thickness);
Crust.MeltingDepth(find(isnan(Crust.MeltingDepth)))=0;
Crust.AlongPlateBoundaryCoordinate=AlongPlateBoundaryCoordinate;

Cryptic.AverageMeltFraction=abs(Cryptic.AverageMeltFraction./Cryptic.Thickness);
Cryptic.AverageMeltFraction(find(isnan(Cryptic.AverageMeltFraction)))=0;
Cryptic.MaxMeltFraction=abs(Cryptic.MaxMeltFraction./Cryptic.Thickness);
Cryptic.MaxMeltFraction(find(isnan(Cryptic.MaxMeltFraction)))=0;
Cryptic.MeltingDepth=abs(Cryptic.MeltingDepth./Cryptic.Thickness);
Cryptic.MeltingDepth(find(isnan(Cryptic.MeltingDepth)))=0;
Cryptic.AlongPlateBoundaryCoordinate=AlongPlateBoundaryCoordinate;

warning on all;

return