function [Parameter,Result,ResultSummary]=crustCalculation(NameTag,Geometry,Res,Shoulder,meltFunction,SpreadingDirection,SpreadingRate_Half,...
    TrimDistance,ExtractionWidth,ExtractionDepth,ExtractionSlope,MeltFractionCutoff,SmoothingWidth,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crustCalculation.m
% Estimate melt flux and crustal thickness at model
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   NameTag                 : Model name tag
%   Geometry
%       |......             : Geometry informaiton of model
%   Res
%       |......             : Resolution information of model
%   Shoulder(ind)
%       |......             : Melt shoulder information of model
%   MeltFunction            : Melt fraction function
%   SpreadingDirection      : Direction of spreading [deg] 0: x-direction, 90: y-direction 
%   Spreading_Half          : Half-spreading rate [m/s]
%   TrimDistance            : Distance to trim close to boundary [km]
%   ExtractionWidth         : Melt extraction width [km]
%   ExtractionDepth         : Melt extraction depth [km]
%   ExtractionSlope         : Melt transport critical slope
%   MeltFractionCutoff      : Retained melt fraction 
%   SmoothingWidth          : Crust-level melt redistribution length (horizontal dike propagation distance, or smoothing width) [km]
%   indFigure               : Figure index
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Parameter
%       |.ExtractionWidth   : Melt extraction width [km]
%       |.ExtractionDepth   : Melt extraction depth [km]
%       |.ExtractionSlope   : Melt transport critical slope
%       |.MeltFractionCutoff : Retained melt fraction 
%   Result(ind)
%       |.Shoulder
%           |.Swath(ind)
%               |.MeltFraction : Melt fraction at sampling columns
%               |.MaxMeltFraction : Maximum melt fraction at sampling columns
%               |.AverageMeltFraction : Average melt fraction at sampling columns 
%               |.MeltProduction : Melt production at sampling columns [1/s]
%               |.MeltFlux  : Melt flux at sampling columns [m/s]
%               |.MeltingDepth : Melting depth at sampling columns [km]
%               |......
%           |.Distribution
%               |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km] 
%               |.Thickness	: Distribution of erupted crust along plate boundary from each shoulder [km]
%               |.AverageMeltFraction : Average melt fraction along plate boundary from each shoulder
%               |.MaxMeltFraction : Maximum melt fraction along plate boundary from each shoulder
%               |.MeltingDepth : Melting depth along plate boundary from each shoulder [km]
%               |.CrypticThickness : Distribution of cryptic crust along plate boundary from each shoulder [km]
%               |.CrypticAverageMeltFraction : Average cryptic melt fraction along plate boundary from each shoulder
%               |.CrypticMaxMeltFraction : Maximum cryptic melt fraction along plate boundary from each shoulder
%               |.CrypticMeltingDepth : Cryptic melting depth along plate boundary from each shoulder [km]
%           |......
%       |.Crust
%           |.Thickness     : Distribution of erupted crust along plate boundary [km] 
%           |.AverageMeltFraction : Average melt fraction along plate boundary
%           |.MaxMeltFraction : Maximum melt fraction along plate boundary 
%           |.MeltingDepth	: Melting depth along plate boundary [km] 
%           |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%       |.Cryptic
%           |.Thickness     : Distribution of cryptic crust along plate boundary [km] 
%           |.AverageMeltFraction : Average cryptic melt fraction along plate boundary
%           |.MaxMeltFraction : Maximum cryptic melt fraction along plate boundary
%           |.MeltingDepth	: Cryptic melting depth along plate boundary [km]
%           |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%       |.History
%           |.Thickness     : Accummulated crustal thickness in the tile [km]
%           |.TileCenter_x, _y : Tile center coordinates [km]
%           |.AverageSegmentCrust : Averaged accummulated crustal thickness along plate boundaries [km]
%           |.AverageWholdDomain : Averaged crustal thickness in the whole domain [km]
%           |.AverageRidgeDomain : Averaged crustal thickness in the ridge domain [km]
%           |.AverageTransformDomain : Averaged crustal thickness in the transform domain [km]
%           |.AverageObliqueDomain : Averaged crustal thickness in the oblique domain [km]
%           |.AccretWidth_Ridge : Accretion zone width along ridge segment [km]
%           |.AccretWidth_Transform : Accretion zone width along transform faults [km]
%           |.SmoothingWidth : Smoothing width [km]
%           |.CellSize      : Size of tiles for crustal thickness calculation [km] 
%   ResultSummary           : A summary of key data from model
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   nCount                  : Total number of iteration
%   iCount                  : Iteration index
%   iFractionCut            : Index for melt fraction cutoff parameter iteration 
%   iExtractWidth           : Index for melt extraction width parameter iteration 
%   iExtractDepth           : Index for melt extraction depth parameter iteration 
%   iExtractSlope           : Index for critical melt transport slope parameter iteration 
%   iSmooth                 : Index for smoothing width parameter iteration 
%   Crust                   : Erupted crust informaiton from each iteration
%   Cryptic                 : Cryptic crust information from each iteration
%   SeafloorHistory         : Crustal thickness distribution in model domain from each iteration 
%   FiltedCrust             : Filted erupted crust information along plate boundary 
%   FiltedCryptic           : Filted cryptic crust information along plate boundary 
%   ReportData              : Model data output from each iteration
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   meltFluxCalculation
%   meltExtraction
%   plotCrustalAccretionProfile
%   crustalThicknessSmoothing
%   crustalHistory
%   dataOutput
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

ResultSummary=[];
iCount=0;
nCount=numel(MeltFractionCutoff)*numel(ExtractionWidth)*numel(ExtractionDepth)*numel(ExtractionSlope); % total number of iterations

% % initialize iteration
for iFractionCut=1:numel(MeltFractionCutoff);
    disp(sprintf('>>>    Working on MeltFractionCutoff=%g',MeltFractionCutoff(iFractionCut)));
    % % calculate melt flux 
    Shoulder=meltFluxCalculation(Shoulder,meltFunction,MeltFractionCutoff(iFractionCut));
    for iExtractWidth=1:numel(ExtractionWidth);
        disp(sprintf('>>>      Working on ExtractionWidth=%g',ExtractionWidth(iExtractWidth)));
        for iExtractDepth=1:numel(ExtractionDepth);
            disp(sprintf('>>>        Working on ExtractionDepth=%g',ExtractionDepth(iExtractDepth)));
            for iExtractSlope=1:numel(ExtractionSlope);
                iCount=iCount+1;
                disp(sprintf('>>>          Working on ExtractionSlope=%g, iteration %g/%g',ExtractionSlope(iExtractSlope),iCount,nCount));
                % % determine melt extraction and calculate crustal accretion profile 
                [Shoulder,Crust,Cryptic]=meltExtraction(Geometry,Res,Shoulder,...
                    ExtractionWidth(iExtractWidth),ExtractionDepth(iExtractDepth),ExtractionSlope(iExtractSlope),SpreadingRate_Half);
                % % store model result
                Parameter(iCount).ExtractionWidth=ExtractionWidth(iExtractWidth);
                Parameter(iCount).ExtractionDepth=ExtractionDepth(iExtractDepth);
                Parameter(iCount).ExtractionSlope=ExtractionSlope(iExtractSlope);
                Parameter(iCount).MeltFractionCutoff=MeltFractionCutoff(iFractionCut);
                Result(iCount).Shoulder=Shoulder;
                Result(iCount).Crust=Crust;
                Result(iCount).Cryptic=Cryptic;
                % % plot crustal accretion profile
                plotCrustalAccretionProfile(Geometry,Result(iCount).Shoulder,Result(iCount).Crust,indFigure);
                print(indFigure,'-dpdf','-r300',sprintf('%s_F%d_CrustProfile.pdf',NameTag,indFigure));

                for iSmooth=1:numel(SmoothingWidth);
                    % % smooth crust information
                    disp(sprintf('>>>            Working on SmoothingWidth=%g',SmoothingWidth(iSmooth)));
                    FiltedCrust=crustalThicknessSmoothing(Geometry,Crust,SmoothingWidth(iSmooth),indFigure+1);
                    print(indFigure+1,'-dpdf','-r300',sprintf('%s_F%d_CrustProfileFilted.pdf',NameTag,indFigure+1));
                    FiltedCryptic=crustalThicknessSmoothing(Geometry,Cryptic,SmoothingWidth(iSmooth),indFigure+2);
                    print(indFigure+2,'-dpdf','-r300',sprintf('%s_F%d_CrypticCrustProfileFilted.pdf',NameTag,indFigure+2));
                    
                    % % calculate crustal thickness map
                    SeafloorHistory=crustalHistory(Geometry,Crust,SpreadingDirection,TrimDistance,SmoothingWidth(iSmooth),5,10,1,indFigure+3);
                    print(indFigure+3,'-dpdf','-r300',sprintf('%s_F%d_CrustalThicknessMap.pdf',NameTag,indFigure+3));
                    print(indFigure+4,'-dpdf','-r300',sprintf('%s_F%d_AveragedCrustalThicknessProfile.pdf',NameTag,indFigure+4));
                    print(indFigure+5,'-dpdf','-r300',sprintf('%s_F%d_CrustalThicknessMap_3D.pdf',NameTag,indFigure+5));
                    Result(iCount).History(iSmooth)=SeafloorHistory;
                    
                    % % report model data
                    ReportData=dataOutput(SpreadingRate_Half,ExtractionWidth(iExtractWidth),ExtractionDepth(iExtractDepth),SmoothingWidth(iSmooth),SeafloorHistory);
                    ResultSummary=[ResultSummary;ReportData];
                end
            end
        end
    end
end
return