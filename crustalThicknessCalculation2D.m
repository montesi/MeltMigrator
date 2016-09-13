function CrustalThickness=crustalThicknessCalculation2D(Model,Grid,Res,Sampling_x,Sampling_y,Sampling_z,Lid,meltFunction,lidTemperature,...
    SpreadingRate_Half,MeltFractionCutoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crustalThicknessCalculation2D.m;
% Integrate melt produced in specified 2D section
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% June 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Model                   : 3D finite element structure from COMSOL Multiphysics 4.x; Spatial dimension in [km]
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates
%   Res
%       |.MeltCalibration   : Melt calibration resolution [km]
%       |......
%   Sampling_
%       |x, y, z            : Sampling vectors of the target section
%   Lid
%       |.T                 : Temperature of lid [degC]
%       |.Solidus           : Solidus temperature at lid [degC]
%       |......
%   MeltFunction            : Melt fraction function
%   LidTemperature          : Lid temperature function
%   SpreadingRate_Half      : Half-spreading rate [m/s]
%   MeltFractionCutoff      : Retained melt fraction
% -------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
% 	CrustalThickness        : Crustal thickness [km]
% -------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   Grid_             
%       |x, y, z            : 2D matrices of section sampling coordinates
%   Temperature             : Temperature of the section [degC]
%   VerticalVelocity        : Vertical velocity of the section [m/s]
%   IsLidTemperatureAboveSolidus
%   MeltFraction            : Melt fraction
%   MeltProduction          : Melt production [1/s]
%   MeltFlux                : Melt flux [m/s]
%   GradientMeltFraction_
%       |x, z               : Melt fraction gradient in x and z direction
%   SamplingSpacing         : Spacing along ridge axis [m]
% -------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
% 	None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

[Grid_x,Grid_z]=meshgrid(Sampling_x,Sampling_z);
[Grid_y,Grid_z]=meshgrid(Sampling_y,Sampling_z);

[Temperature,VerticalVelocity]=mphinterp(Model,{'T','w'},'coord',[Grid_x(:),Grid_y(:),-Grid_z(:)]',... 
    'unit',{'degC','m/s'},'dataset','dset1','solnum','end');
Temperature=reshape(Temperature,size(Grid_x));
VerticalVelocity=reshape(VerticalVelocity,size(Grid_x));

IsLidTemperatureAboveSolidus=interp2(Grid.x,Grid.y,Lid.T-Lid.T_Solidus,Sampling_x,Sampling_y)>0; % determine if lid temperature is above solidus

MeltFraction=max(meltFunction(Grid_z,Temperature),0); 
% [GradientMeltFraction_x,GradientMeltFraction_z]=gradient(MeltFraction,Sampling_x.*1000,Sampling_z.*1000);
% MeltProduction=max(-GradientMeltFraction_z.*VerticalVelocity,0);
% MeltFlux=sum(MeltProduction.*repmat(gradient(Grid_z(:,1)),1,size(MeltProduction,2)).*1000.*...
%     (MeltFraction>MeltFractionCutoff).*(Temperature>lidTemperature(Grid_z)),1).*IsLidTemperatureAboveSolidus; 
[GradientMeltFraction_x,GradientMeltFraction_z]=gradient(MeltFraction);
MeltProduction=max(-GradientMeltFraction_z.*VerticalVelocity,0).*(MeltFraction>MeltFractionCutoff); % melt production rate
MeltFlux=sum(MeltProduction).*IsLidTemperatureAboveSolidus; % Melt flux
SamplingSpacing=Res.MeltCalibration.*1000;
CrustalThickness=sum(MeltFlux.*SamplingSpacing/(2*SpreadingRate_Half))/1000;

return