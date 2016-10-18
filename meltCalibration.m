function CalibrationFactor=meltCalibration(Model,Geometry,Grid,Vector,Res,Lid,meltFunction,lidTemperature,...
    SpreadingDirection,SpreadingRate_Half,MeltFractionCutoff,CrustalThickness_Reference,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltCalibration.m
% Calculate a factor to scale the crustal thickness to a reference value
% Hailong Bai
% June 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Model                   : 3D finite element model from COMSOL Multiphysics 4.x
%   Geometry
%       |.PlateBoundaryType : Indicator of plate boundary type. 1: ridge, 2: transform, 3: oblique segment
%       |.PlateBoundaryLength : Length of each plate boundary segment [km]
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |......
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates
%   Vector
%       |.z                 : 1D vectors of depth sampling coordinates
%       |......
%   Res
%       |.MeltCalibration   : Melt calibration resolution [km]
%       |......
%   Lid
%       |.T                 : Temperature of lid [degC]
%       |.Solidus           : Solidus temperature at lid [degC]
%       |......
%   MeltFunction            : Melt fraction function
%   LidTemperature          : Lid temperature function
%   SpreadingDirection      : Direction of seafloor spreading (relative to x axis) [deg]
%   SpreadingRate_Half      : Half-spreading rate [m/s]
%   MeltFractionCutoff      : Retained melt fraction
%   CrustalThickness_Reference : Reference crustal thickness [km]
%   indFigure               : Figure index
% -------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
% 	CalibrationFactor       : Factor used to scale crustal thickness
% -------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   RidgeIndex              : Indexes of the ridge segments
%   LongRidgeIndex          : Index of the longest ridge segment
%   TargetRidgeIndex        : Index of the target ridge segment
%   TargetPoint             : Center coordinates of the target ridge segment 
%   CalibrationSampling_
%       |x, y               : Sampling vectors of the target section
%   SampleLine              : Sampling line function in map view
%   InBoundaryIndex         : Index for the sampling points within the model domain
% -------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
% 	crustalThicknessCalculation2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% % find 2D section in the middle of the longest ridge segment
RidgeIndex=find(Geometry.PlateBoundaryType==1); % locate ridge segments along plate boundaries
LongRidgeIndex=find(Geometry.PlateBoundaryLength(RidgeIndex)==max(Geometry.PlateBoundaryLength(RidgeIndex))); % locate the longest ridge segment
TargetRidgeIndex=RidgeIndex(LongRidgeIndex(1)); % index of target ridge segment
TargetPoint=[mean([Geometry.PlateBoundary.x(TargetRidgeIndex),Geometry.PlateBoundary.x(TargetRidgeIndex+1)]),...
    mean([Geometry.PlateBoundary.y(TargetRidgeIndex),Geometry.PlateBoundary.y(TargetRidgeIndex+1)])]; % center of target ridge

if SpreadingDirection==0; % if spreading is along x-direction
    CalibrationSampling_x=[Geometry.ModelBoundary.x(1):Res.MeltCalibration:Geometry.ModelBoundary.x(2)];
    CalibrationSampling_y=linspace(TargetPoint(2),TargetPoint(2),length(CalibrationSampling_x));
elseif SpreadingDirection==90; % if spreading is along y-direction
    CalibrationSampling_y=[Geometry.ModelBoundary.y(1):Res.MeltCalibration:Geometry.ModelBoundary.y(2)];
    CalibrationSampling_x=linspace(TargetPoint(1),TargetPoint(1),length(CalibrationSampling_y));
else % if spreading is along other direction
    sampleLine=@(x) tand(SpreadingDirection).*x+TargetPoint(2)-tand(SpreadingDirection).*TargetPoint(1);
    CalibrationSampling_x=[Geometry.ModelBoundary.x(1):Res.MeltCalibration:Geometry.ModelBoundary.x(2)];
    CalibrationSampling_y=sampleLine(CalibrationSampling_x);
    InBoundaryIndex=find((CalibrationSampling_y>=Geometry.ModelBoundary.y(1))&(CalibrationSampling_y<=Geometry.ModelBoundary.y(2)));
    CalibrationSampling_x=CalibrationSampling_x(InBoundaryIndex);
    CalibrationSampling_y=CalibrationSampling_y(InBoundaryIndex);
end

figure(indFigure);
plot3(CalibrationSampling_x,CalibrationSampling_y,zeros(size(CalibrationSampling_x)),'b','linewidth',2);

CalibrationFactor=crustalThicknessCalculation2D(Model,Grid,Res,CalibrationSampling_x,CalibrationSampling_y,Vector.z,Lid,...
    meltFunction,lidTemperature,SpreadingRate_Half,MeltFractionCutoff)/CrustalThickness_Reference;

return
