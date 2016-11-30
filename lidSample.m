function [Lid,Box]=lidSample(Model,Box,Grid,Vector,Geometry,solidus,lidDepth,Switch_UseCOMSOL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidSample.m
% Sample COMSOL model and extract lid information
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% June 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Model                   : 3D model results file
%   Box
%       |.x, .y, .z         : 3D matrices of sampling coordinates
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates
%   Vector
%       |.x, .y, .z         : 1D vectors of sampling coordinates
%   Geometry
%       |.ModelBoundary_z   : Model depth limit
%       |......
%   Solidus                 : Mantle solidus function
%   LidDepth                : Function to find depth for lid
%   Switch_UseCOMSOL        : Switch for COMSOL model usage, 1: using COMSOL, 2: not using COMSOL
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
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
%   Box
%       |.T                 : Temperature data of model [degC]
%       |.Velocity          : Velocity magnitude data of model [m/s]
%           |_x, _y, _z     : Velocity conponents in x-, y-, z-direciton [m/s]
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Data Extraction

% % read 3D temperature and velocity data
if Switch_UseCOMSOL
    [Box.T,Box.Velocity,Box.Velocity_x,Box.Velocity_y,Box.Velocity_z]=mphinterp(Model,{'T','spf.U','u','v','w'},'coord',[Box.x(:)';Box.y(:)';-Box.z(:)'],...
        'unit',{'degC','m/s','m/s','m/s','m/s'},'dataset','dset1','solnum','end');
else
    interpolant_vx=scatteredInterpolant(Model{1},Model{2},Model{3},Model{4});
    interpolant_vy=scatteredInterpolant(Model{1},Model{2},Model{3},Model{5});
    interpolant_vz=scatteredInterpolant(Model{1},Model{2},Model{3},Model{6});
    interpolant_T=scatteredInterpolant(Model{1},Model{2},Model{3},Model{7});
    Box.T=interpolant_T(Box.x(:)',Box.y(:)',-Box.z(:)')-273.15;
    Box.Velocity_x=interpolant_vx(Box.x(:)',Box.y(:)',-Box.z(:)');
    Box.Velocity_y=interpolant_vy(Box.x(:)',Box.y(:)',-Box.z(:)');
    Box.Velocity_z=interpolant_vz(Box.x(:)',Box.y(:)',-Box.z(:)');
    Box.Velocity=sqrt(Box.Velocity_x.^2+Box.Velocity_y.^2+Box.Velocity_z.^2);
end
Box.T=reshape(Box.T,size(Box.x));
Box.Velocity=reshape(Box.Velocity,size(Box.x));
Box.Velocity_x=reshape(Box.Velocity_x,size(Box.x));
Box.Velocity_y=reshape(Box.Velocity_y,size(Box.x));
Box.Velocity_z=reshape(Box.Velocity_z,size(Box.x));

%% Lid Information Extraction

% % calculate the lid depth
for ix=1:length(Vector.x);
    for iy=1:length(Vector.y);
        Lid.Depth(iy,ix)=lidDepth(Vector.z,reshape(Box.T(iy,ix,:),size(Vector.z))); 
    end
end
Lid.Depth(find(Lid.Depth>max(Geometry.ModelBoundary.z)))=max(Geometry.ModelBoundary.z);
Lid.Depth=reshape(Lid.Depth,size(Grid.x));

% % read 2D information of lid
if Switch_UseCOMSOL
    [Lid.T,Lid.TGrad_x,Lid.TGrad_y]=mphinterp(Model,{'T','Tx','Ty'},'coord',[Grid.x(:)';Grid.y(:)';-Lid.Depth(:)'],...
        'unit','degC','dataset','dset1','solnum','end');
else
    Lid.T=interpolant_T(Grid.x(:)',Grid.y(:)',-Lid.Depth(:)')-273.15;
end
Lid.T=reshape(Lid.T,size(Grid.x));
% Lid.TGrad_x=reshape(Lid.TGrad_x,size(Grid.x)); % temperature gradient, x-direction
% Lid.TGrad_y=reshape(Lid.TGrad_y,size(Grid.x)); % temperature gradient, y-direction
Lid.T_Solidus=solidus(Lid.Depth); % solidus temperature at lid

[Lid.DepthGrad_x,Lid.DepthGrad_y]=gradient(Lid.Depth,Grid.x(1,2)-Grid.x(1,1),Grid.y(2,1)-Grid.y(1,1)); % lid depth gradient in x & y direction
Lid.Slope=(Lid.DepthGrad_x.^2+Lid.DepthGrad_y.^2).^(1/2); % lid slope
[Lid.SlopeGrad_x,Lid.SlopeGrad_y]=gradient(Lid.Slope,Grid.x(1,2)-Grid.x(1,1),Grid.y(2,1)-Grid.y(1,1)); % lid slope gradient in x & y direction

Lid.ExcessT=max(Box.T-solidus(Box.z),[],3); % excess temperature above solidus
[Lid.ExcessTGrad_x,Lid.ExcessTGrad_y]=gradient(Lid.ExcessT,Grid.x(1,2)-Grid.x(1,1),Grid.y(2,1)-Grid.y(1,1)); % excess temperature gradient in x & y direction

return