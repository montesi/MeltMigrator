%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setParameters.m
% Set key parameters for MeltMigrator
% Hailong Bai
% October 2016
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   None
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   See variables defined below
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   None
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   meltFunctionRJ1981
%   meltFunctionMELTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Info

Switch_UseCOMSOL=0; % switch for COMSOL model usage, 1: using COMSOL, 2: not using COMSOL
% ModelName='MeltMigrationDemo_01.mph'; % name of the model being processed
ModelName='DemoModelResult.txt'; % name of the model being processed
NameTag='DEMO'; % name tag for the saved results
iFigure=1; % starting ID for figures

%% Spreading Info

SpreadingRate_Plate=2.0; % spreading rate of one plate [cm/yr]
SpreadingRate_etalP=2.0; % spreading rate of other plate [cm/yr]
SpreadingDirection=0; % spreading direction [deg], 0: x-direction, 90: y-direction

%% Domain Geometry Info

Geometry.ModelBoundary.x=[0,600]; % model domain limit in x direction [km]
Geometry.ModelBoundary.y=[0,550]; % model domain limit in y direction [km]
Geometry.ModelBoundary.z=[0,100]; % model domain limit in z direction [km]
TrimDistance=20; % trimming distance from model boundaries to avoid boundary effects [km]

%% Plate Boundary Info

Geometry.PlateBoundary.x=[200,200,300,300,400,400]; % plate boundary way points x coordinates [km]
Geometry.PlateBoundary.y=[550,400,400,250,150,0]; % plate boundary way points y coordinates [km]
Geometry.PlateBoundaryType=[1,2,1,3,1]; % plate boundary segment types, 1: ridge, 2: transform, 3: oblique segment

%% Melt Function Info

AdiabaticGradient=0.6475; % adiabatic gradient [degC/km];
solidus=@(z) 1100+3*z; % solidus function, Reid & Jackson, 1981
% solidus=@(z) (1120+4.3092*z-5.2224e-3*z.^2-AdiabaticGradient*z); % alternative solidus function, Hirschmann, 2000
% solidus=@(z) (1085.7+4.30*z-5.33e-3*z.^2); % alternative solidus function, Katz, 2003
meltFunction=@(z,T) meltFunctionRJ1981(z,T); % melt function, Reid & Jackson, 1981
% meltFunction=@(z,T) meltFunctionMELTS(z,T); % alternative melt function, using alphaMELTS
MeltFractionCutoff=0.01; % retained melt fraction

%% Melt Function Calibration Info

Switch_MeltCalibration=1; % switch for melt calibration
CrustalThickness_Reference=6; % reference crustal thickness [km]
% CalibrationFactor=1.7403; % calibration factor for melt function

%% Saddle Search Info

Switch_SaddleSelectionByMouse=0; % switch for saddle search starting points selection method, 0: to be specified by user, 1: to be selected by mousing using ginput
if ~Switch_SaddleSelectionByMouse % if not using mouse, specify the starting points
    SaddleInitialPoint=[250,400;350,200]; % starting points for saddle point search
end
SaddleSearchRadius=10; % radius of saddle point search
SaddleFlatCutoff=0.001; % criterion for saddle flatness
LidDepthLimit=60; % lower limit for permeability barrier [km]

%% Melt Migration & Extraction Info

T_MeltSeed=170; % temperature for melt seeds generation [degC]
ExtractionWidth=[4]; % extraction width [km]
ExtractionDepth=[20]; % extraction depth [km]
ExtractionSlope=[0.1]; % critical slope [km]
SmoothingWidth=[60]; % smoothing intensity [km]

%% Resolution Info

Res.x=5; % x resolution [km]
Res.y=5; % y resolution [km]
Res.nz=35; % sampling size, z direction
Res.zFactor=1.3; % z nonlinear factor
Res.MeltCalibration=10; % melt calibration resolution [km]
Res.Saddle=1; % saddle line resolution [km]
Res.MeltSeed=5; % melt line seed resolution [km]
Res.MeltTrajectory=5; % melt line resolution [km]
Res.MeltSwath=4; % melt swath resolution [km]
Res.Axis=1; % plate boundaries resolution for crustal thickness calculation [km]
