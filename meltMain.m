%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltMain.m
% Melt Migration Calculation Routine
% Postprocess mid-ocean ridge model to calculate melt trajectories, melt flux & crustal thickness
% Laurent Montesi, Hailong Bai after Jennifer Barry, Mark Behn, Laura Hebert
% Department of Geology, University of Maryland, College Park, MD, 20740
% Version 3.1, October 2016
%--------------------------------------------------------------------------
% PRE-EXISTING ------------------------------------------------------------
%   .mph                    : 3D finite element model from COMSOL Multiphysics 4.x;
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   Model                   : 3D finite element model from COMSOL Multiphysics 4.x
%   NameTag                 : Tag for current model
%   SpreadingRate_Plate     : Spreading rate of left-hand plate [cm/yr]
%   SpreadingRate_etalP     : Spreading rate of right-hand plate [cm/yr]
%   SpreadingRate_Full      : Full spreading rate [cm/yr]
%   SpreadingRate_Half      : Half-spreading rate [m/s]
%   SpreadingDirection      : Direction of seafloor spreading (relative to x axis) [deg]
%   AdiabaticGradient       : Adiabatic gradient [degC/km]
%   TrimDistance            : Distance to trim close to boundary [km]
%   Switch_MeltCalibration  : Switch for melt calibration. 1: perform, 0: do not perform
%   CrustalThickness_Reference : Reference crustal thickness [km]
% 	CalibrationFactor       : Factor used to scale crustal thickness
%   Switch_SaddleSelectionByMouse : Switch for saddle starting point selection method, 0: specified by user, 1: selected by mousing using ginput
%   SaddleInitialPoint      : Starting point for saddle point search [km]
%   SaddleSearchRadius      : Radius for saddle line search [km]
%   SaddleFlatCutoff        : Criterion for saddle flatness
%   LidDepthLimit           : Lower limit for permeability barrier [km]
%   T_MeltSeed              : Desired countour of ExcessT for seed generation [degC]
%   MeltFractionCutoff      : Retained melt fraction
%   ExtractionWidth         : Melt extraction width [km]
%   ExtractionDepth         : Melt extraction depth [km]
%   ExtractionSlope         : Melt transport critical slope
%   SmoothingWidth          : Crust-level melt redistribution length (horizontal dike propagation distance, or smoothing width) [km]
%   iFigure                 : Index of first figure
%   Solidus                 : Mantle solidus function
%   LidTemperature          : Lid temperature function
%   LidDepth                : Function to find depth for lid
%   MeltFunction            : Melt fraction function
%   Geometry
%       |.ModelBoundary     : Model geometry corner coordinates [km]
%       |.PlateBoundary     : Plate boundary segment coordinates [km]
%       |.PlateBoundaryType : Indicator of plate boundary type. 1: ridge, 2: transform, 3: oblique segment
%       |.PlateBoundaryLength_x, _y : Plate boundary segment length in x- and y-direction [km]
%       |.PlateBoundaryLength : Length of each plate boundary segment [km]
%       |.PlateBoundaryAngleRad : Angle of each plate boundary segment [rad]
%       |.PlateBoundaryAngleDeg : Angle of each plate boundary segment [deg]
%   Res
%       |.x, .y             : Sampling resolution in x- and y-direction [km]
%       |.nx, .ny, .nz      : Sampline size in x-, y- and z-direction
%       |.zFactor           : Nonlinear sampling factor in z-direction
%       |.MeltCalibration   : Melt calibration resolution [km]
%       |.Saddle            : Saddle line resolution [km]
%       |.MeltSeed          : Distance between adjacent melt seeds [km]
%       |.MeltTrajectory    : Distance between adjacent points on melt trajectory [km]
%       |.nMeltSwath        : Melt swath sampling size
%       |.MeltSwath         : Melt swath resolution [km]
%       |.Axis              : Plate boundaries resolution for crustal thickness calculation [km]
%   Vector
%       |.x, .y, .z         : 1D vectors of sampling coordinates
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates
%   Box
%       |.T                 : Temperature data of model [degC]
%       |.Velocity          : Velocity magnitude data of model [m/s]
%           |_x, _y, _z     : Velocity conponents in x-, y-, z-direciton [m/s]
%       |.SymmetricT        : Temperature data of reference symmetric model [degC]
%       |.SymmetricVelocity : Velocity magnitude data of reference symmetric model [m/s]
%           |_x, _y, _z     : Velocity conponents in x-, y-, z-direciton [m/s]
%       |.SymmetricMeltFraction : Melt fraction in symmetric model
%       |.SymmetricMeltProduction : Melt production in symmetric model [1/s]
%       |.SymmetricMeltFlux : Melt flux in symmetric model [m/s]
%       |.MeltFraction      : Melt fraction in current asymmetric model
%       |.MeltProduction    : Melt production in current asymmetric model [1/s] 
%       |.MeltFlux          : Melt flux in current asymmetric model [m/s]
%       |.DiffT             : Temperature difference between asymmetric and symmetric models [degC]
%       |.DiffVelocity      : Velocity difference between asymmetric and symmetric models [m/s]
%           |_x, _y, _z     : Velocity difference conponents in x-, y-, z-direciton [m/s]
%       |.DiffMeltFraction  : Melt fraction difference between asymmetric and symmetric models
%       |.DiffMeltProduction : Melt production difference between asymmetric and symmetric models [1/s]   
%       |.DiffMeltFlux      : Melt flux difference between asymmetric and symmetric models [m/s]
%   Lid
%       |.Depth             : Depth of lid [km]
%       |.DepthGrad         : Slope of lid in x or y direction
%       |.Slope             : Slope of lid
%       |.SlopeGrad         : Slope gradient in x or y direction
%       |.T                 : Temperature of lid [degC]
%       |.TGrad             : Horizontal thermal gradient [degC/km]
%       |.T_Solidus         : Solidus temperature at lid [degC]
%       |.ExcessT           : Excess temperature above solidus [degC]
%       |.ExcessTGrad       : Horizontal gradient of ExcessT [degC/km]
%   Saddle(ind)
%       |.Point             : Position of saddle point [km]
%       |.Down              : Downgoing saddle line coordinates [km]
%       |.DisDown           : Distance along downgoing saddle line [km]
%       |.Up                : Upgoing saddle line coordinates [km]
%       |.DisUp             : Distance along upgoing saddle line [km]
%   Seed(ind)
%       |.x, .y             : Coordinates of melt line seeds [km]
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
%       |.Swath(ind)
%           |.Polygon(ind)
%               |.x,.y      : Coordinates of polygon outline [km]
%           |.PolygonArea   : Area of each polygon [km^2]
%           |.SwathArea     : Area of swath [km^2]
%           |.PolygonCenter
%               |_x, _y, _z : Coordinates of polygon centers [km]
%               |_T         : Temperature at polygon centers [degC]
%           |.PolygonSlope  : Slope at polygon centers
%           |.Column
%               |_x, _y, _z : Coordinate at vertical sampling column below polygon [km]
%               |_T         : Temperature at sampling column [degC]
%               |_VerticalVelocity : Vertical velocity at sampling column [m/s]
%   Parameter
%       |.ExtractionWidth   : Melt extraction width [km]
%       |.ExtractionDepth   : Melt extraction depth [km]
%       |.ExtractionSlope   : Melt transport critical slope
%       |.MeltFractionCutoff : Retained melt fraction 
%   Result(ind)
%       |.Shoulder(ind)
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
%               |.Thickness : Distribution of erupted crust along plate boundary from each shoulder [km]
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
%           |.MeltingDepth  : Melting depth along plate boundary [km] 
%           |.AlongPlateBoundaryCoordinate : Coordinate along plate boundary [km]
%       |.Cryptic
%           |.Thickness     : Distribution of cryptic crust along plate boundary [km] 
%           |.AverageMeltFraction : Average cryptic melt fraction along plate boundary
%           |.MaxMeltFraction : Maximum cryptic melt fraction along plate boundary
%           |.MeltingDepth  : Cryptic melting depth along plate boundary [km]
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
% ATTENDING SCRIPTS -------------------------------------------------------
%   geometryTrim
%   meltFunctionRJ1981
%   meltFunctionMELTS
%   lidSample
%   plotLidInfo
%   meltCalibration
%       | crustalThicknessCalculation2D
%   asymmetryStudy
%   saddlePreparation
%       | alongSlopeSlopeGrad
%       | alongSlopeGrad
%   meltTrajectorySeed
%       | crossToGrad
%       | inBounds
%   meltTrajectory
%       | lineTruncate
%       | lineSegmentIntersect
%       | alongSlopeGrad
%       | lineCleanUp
%           | lineSegmentIntersect
%       | lineStore
%   meltSwath
%       | lineSampleDepth
%           | assignSegment
%       | clockwiseAreaCalc
%   crustCalculation
%       | meltFluxCalculation
%       | meltExtraction
%           | extractionDetermination
%               | minDistBetweenTwoPolygons
%           | assignSegment
%           | minDistBetweenTwoPolygons
%       | plotCrustalAccretionProfile
%       | crustalThicknessSmoothing
%           | fastSmooth
%       | crustalHistory
%           | fastSmooth
%           | minDistBetweenTwoPolygons
%           | plotCrustalThicknessMapIn3D
%       | dataOutput
%--------------------------------------------------------------------------
% SCRIPTS FROM MATLAB CENTRAL ---------------------------------------------
%   lineSegmentIntersect.m          Copyright 2010, U. Murat Erdem
%   minDistBetweenTwoPolygons.m     Copyright 2009, Guillaume Jacquenot
%   fastSmooth.m                    Copyright 2008, Tom O'Haver
%--------------------------------------------------------------------------
% NAMING RULES ------------------------------------------------------------
%   T                   : Temperature
%   P                   : Pressure
%   x                   : Length Coordinate
%   y                   : Width Coordinate
%   z                   : Depth Coordinate
%   Lid                 : Permeability Barrier
%   Res                 : Resolution
%   n                   : Number
%   i(or ind)           : Index
%   Max                 : Maximum
%   Min                 : Minimum
%   Dis                 : Distance
%   Grad                : Gradient
%   Trajectory, Line    : Melt trajectory, used in script interchangeably
%   Variables are named with the first letter in upper-case  (except for i)
%   Functions are named with the first letter in lower-case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Initialization

disp('>>> Melt migration calculation routine started');
setParameters;

% % load model
tic;
Model=mphload(ModelName);
fprintf('>>> Mid-ocean ridge model loaded in %g s\n',toc);

% % process spreading info
SpreadingRate_Full=SpreadingRate_Plate+SpreadingRate_etalP; % [cm/yr]
SpreadingRate_Half=0.5*SpreadingRate_Full*(3.17e-10); % [m/s]

%% Model & Segment Geometry [km]

Geometry=geometryTrim(Geometry,TrimDistance); % trim geometry to avoid boundary effects
% % plate boundary length
Geometry.PlateBoundaryLength_x=diff(Geometry.PlateBoundary.x); 
Geometry.PlateBoundaryLength_y=diff(Geometry.PlateBoundary.y); 
Geometry.PlateBoundaryLength=sqrt(Geometry.PlateBoundaryLength_x.^2+Geometry.PlateBoundaryLength_y.^2); 
% % plate boundary angle
Geometry.PlateBoundaryAngleRad=atan(Geometry.PlateBoundaryLength_y./Geometry.PlateBoundaryLength_x)+pi*(Geometry.PlateBoundaryLength_x<0); % [rad]
Geometry.PlateBoundaryAngleDeg=Geometry.PlateBoundaryAngleRad*180/pi; % [deg]
        
fprintf('>>> Model size: %d km x %d km x %d km\n',diff(Geometry.ModelBoundary.y),diff(Geometry.ModelBoundary.x),diff(Geometry.ModelBoundary.z));

% % plot geometry
figure(iFigure); clf; hold on;
plot3(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,zeros(size(Geometry.PlateBoundary.x)),'k','linewidth',2);
box on; axis equal;
set(gca,'boxstyle','full');
set(gca,'zdir','reverse');
xlim([Geometry.ModelBoundary.x]);
ylim([Geometry.ModelBoundary.y]);
zlim([Geometry.ModelBoundary.z]);
view(3); camproj('perspective');
title('Model Geometry and Mesh');

%% Resolution

Res.nx=ceil(diff(Geometry.ModelBoundary.x)/Res.x)+1; % sampling size, x direction
Res.ny=ceil(diff(Geometry.ModelBoundary.y)/Res.y)+1; % sampling size, y direction
Res.nMeltSwath=ceil(100/Res.MeltSwath); % melt swath sampling size

Vector.x=linspace(Geometry.ModelBoundary.x(1),Geometry.ModelBoundary.x(2),Res.nx);
Vector.y=linspace(Geometry.ModelBoundary.y(1),Geometry.ModelBoundary.y(2),Res.ny);
Vector.z=linspace(0,max(Geometry.ModelBoundary.z)^(1/Res.zFactor),Res.nz).^Res.zFactor; % nonlinear vertical sampling vector
[Grid.x,Grid.y]=meshgrid(Vector.x,Vector.y);
[Box.x,Box.y,Box.z]=meshgrid(Vector.x,Vector.y,Vector.z);
fprintf('>>> Sampling size: %d x %d x %d = %d points\n',size(Box.x),numel(Box.x));

% % plot mesh
figure(iFigure);
mesh(Grid.x,Grid.y,zeros(size(Grid.x)));
[GridTemp.x,GridTemp.z]=meshgrid(Vector.x,Vector.z);
mesh(GridTemp.x,zeros(size(GridTemp.x))+TrimDistance,GridTemp.z);
[GridTemp.y,GridTemp.z]=meshgrid(Vector.y,Vector.z);
mesh(zeros(size(GridTemp.y))+TrimDistance,GridTemp.y,GridTemp.z);
if ~Switch_SaddleSelectionByMouse
    for iSaddlePoint=1:size(SaddleInitialPoint,1)
        scatter(SaddleInitialPoint(iSaddlePoint,1),SaddleInitialPoint(iSaddlePoint,2),250,'fill','markeredgecolor','k','markerfacecolor','r','linewidth',1.5);
    end
end
clear GridTemp;

%% Lid Sampling

tic;
fprintf('>>> Sampling permeability barrier\n');
lidTemperature=@(z) (1240+(z)*1.9-AdiabaticGradient*z); % permeability barrier temperature [degC]
lidDepth=@(zs,Ts) fsolve(@(z)interp1(zs,Ts,z)-lidTemperature(z),0,optimset('display','off')); % permeability barrier depth [km]
[Lid,Box]=lidSample(Model,Box,Grid,Vector,Geometry,solidus,lidDepth); % sample lid
save(sprintf('%s_Lid.mat',NameTag),'Lid');
save(sprintf('%s_Box.mat',NameTag),'Box');
plotLidInfo(Geometry,Grid,Lid,iFigure+1); % plot lid information
print(iFigure+1,'-dpdf','-r300',sprintf('%s_F%d_Lid.pdf',NameTag,iFigure+1));
fprintf('>>> Permeability barrier sampled in %g s\n',toc);

%% Melting Calibration

if Switch_MeltCalibration
    tic;
    disp('>>> Calibrating melting function');
    CalibrationFactor=meltCalibration(Model,Geometry,Grid,Vector,Res,Lid,meltFunction,lidTemperature,...
        SpreadingDirection,SpreadingRate_Half,MeltFractionCutoff(1),CrustalThickness_Reference,iFigure);
    meltFunction=@(z,T) meltFunction(z,T)./CalibrationFactor; % calibrate melting function to generate reasonable crustal thickness 
    fprintf('>>> Obtained %g km of crust at reference section with CalibrationFactor=%g, in %g s\n',CrustalThickness_Reference,CalibrationFactor,toc);
end

%% Saddle Calculation

tic;
fprintf('>>> Performing saddle calculation\n');
if Switch_SaddleSelectionByMouse
    figure(iFigure+1);
    SaddleInitialPoint=ginput; % select starting points for saddle point search by mouse click on lid maps
end
Saddle=saddlePreparation(Geometry,Grid,Res,Lid,SaddleInitialPoint,SaddleSearchRadius,SaddleFlatCutoff,iFigure+4);
save(sprintf('%s_Saddle.mat',NameTag),'Saddle');
print(iFigure+4,'-dpdf','-r300',sprintf('%s_F%d_Saddle.pdf',NameTag,iFigure+4));
fprintf('>>> Saddle completed in %gs\n',toc);

%% Melt Line Seeds Generation

tic;
fprintf('>>> Generating melt seeds\n'); 
Seed=meltTrajectorySeed(Geometry,Grid,Res,Lid,T_MeltSeed,SpreadingDirection,iFigure+5);
save(sprintf('%s_Seed.mat',NameTag),'Seed');
print(iFigure+5,'-dpdf','-r300',sprintf('%s_F%d_Seed.pdf',NameTag,iFigure+5));
fprintf('>>> Melt seeds generated in %gs\n',toc);

%% Melt Line Calculation

tic;
fprintf('>>> Calculating melt trajectories\n'); 
Shoulder=meltTrajectory(Geometry,Grid,Res,Lid,Saddle,Seed,iFigure+6);
save(sprintf('%s_Shoulder.mat',NameTag),'Shoulder');
print(iFigure+6,'-dpdf','-r300',sprintf('%s_F%d_Shoulder.pdf',NameTag,iFigure+6));
fprintf('>>> Melt trajectories calculated in %gs\n',toc);

%% Melt Swath Generation

tic;
fprintf('>>> Creating melt swaths and polygons\n'); 
Shoulder=meltSwath(Model,Geometry,Grid,Res,Lid,Shoulder,LidDepthLimit,iFigure+7);
save(sprintf('%s_Shoulder.mat',NameTag),'Shoulder');
print(iFigure+7,'-dpdf','-r300',sprintf('%s_F%d_Swath.pdf',NameTag,iFigure+7));
fprintf('>>> Melt swaths and polygons created in %gs\n',toc);

%% Crustal Thickness Calculation

tic;
fprintf('>>> Calculating melt flux and crustal thickness\n'); 
[Parameter,Result,ResultSummary]=crustCalculation(NameTag,Geometry,Res,Shoulder,meltFunction,SpreadingDirection,SpreadingRate_Half,...
    TrimDistance,ExtractionWidth,ExtractionDepth,ExtractionSlope,MeltFractionCutoff,SmoothingWidth,iFigure+8);
save(sprintf('%s_Parameter.mat',NameTag),'Parameter');
save(sprintf('%s_Result.mat',NameTag),'Result');
save(sprintf('%s_ResultSummary.mat',NameTag),'ResultSummary');
fprintf('>>> Melt flux and crustal thickness calculated in %gs\n',toc);
disp('>>> Melt migration calculation routine completed');
