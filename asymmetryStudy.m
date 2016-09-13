function Box=asymmetryStudy(SymmetricModel,Geometry,Vector,Box,meltFunction,Depth,indFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% asymmetryStudy.m
% Compare current asymmetric model with a reference symmetric model
% Hailong Bai
% October 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   SymmetricModel          : Reference symmetric model from COMSOL Multiphysics 4.x
%   Geometry
%       |.PlateBoundary.x,y : Plate boundary coordinates [km]
%       |.ModelBoundary.x,y : Model boundary in x and y direction [km]
%       |......
%   Vector
%       |.x, .y, .z         : 1D vectors of sampling coordinates [km]
%   Box
%       |.T                 : Temperature data of current asymmetric model [degC]
%       |.Velocity          : Velocity magnitude data of current asymmetric model [m/s]
%           |_x, _y, _z     : Velocity conponents in x-, y-, z-direciton [m/s]
%   MeltFunction            : Melt fraction function
%   Depth                   : Depth of interest, where comparison results are plotted [km]
%   indFigure               : Figure index
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Box
%       |.SymmetricT        : Temperature data of reference symmetric model [degC]
%       |.SymmetricVelocity : Velocity magnitude data of reference symmetric model [m/s]
%           |_x, _y, _z     : Velocity conponents in x-, y-, z-direciton [m/s]
%       |.SymmetricMeltFraction
%                           : Melt fraction in symmetric model
%       |.SymmetricMeltProduction
%                           : Melt production in symmetric model [1/s]
%       |.SymmetricMeltFlux
%                           : Melt flux in symmetric model [m/s]
%       |.MeltFraction      : Melt fraction in current asymmetric model
%       |.MeltProduction    : Melt production in current asymmetric model [1/s] 
%       |.MeltFlux          : Melt flux in current asymmetric model [m/s]
%       |.DiffT             : Temperature difference between asymmetric and symmetric models [degC]
%       |.DiffVelocity      : Velocity difference between asymmetric and symmetric models [m/s]
%           |_x, _y, _z     : Velocity difference conponents in x-, y-, z-direciton [m/s]
%       |.DiffMeltFraction  : Melt fraction difference between asymmetric and symmetric models
%       |.DiffMeltProduction : Melt production difference between asymmetric and symmetric models [1/s]   
%       |.DiffMeltFlux      : Melt flux difference between asymmetric and symmetric models [m/s]
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%	DiffT                   : Temperature difference between asymmetric and symmetric models at specified depth [degC]
%	DiffVelocity            : Velocity difference between asymmetric and symmetric models at specified depth  [m/s]
%       |_x, _y, _z         : Velocity difference conponents in x-, y-, z-direciton at specified depth [m/s]
%	DiffMeltFraction        : Melt fraction difference between asymmetric and symmetric models at specified depth
%	DiffMeltProduction      : Melt production difference between asymmetric and symmetric models at specified depth [1/s]   
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Data Interpolation & Calculation

[Box.SymmetricT,Box.SymmetricVelocity,Box.SymmetricVelocity_x,Box.SymmetricVelocity_y,Box.SymmetricVelocity_z]=mphinterp(SymmetricModel,...
    {'T','spf.U','u','v','w'},'coord',[Box.x(:)';Box.y(:)';-Box.z(:)'],'unit',{'degC','m/s','m/s','m/s','m/s'},'dataset','dset1','solnum','end');
Box.SymmetricT=reshape(Box.SymmetricT,size(Box.x));
Box.SymmetricVelocity=reshape(Box.SymmetricVelocity,size(Box.x));
Box.SymmetricVelocity_x=reshape(Box.SymmetricVelocity_x,size(Box.x));
Box.SymmetricVelocity_y=reshape(Box.SymmetricVelocity_y,size(Box.x));
Box.SymmetricVelocity_z=reshape(Box.SymmetricVelocity_z,size(Box.x));
Box.SymmetricMeltFraction=(max(meltFunction(Box.z,Box.SymmetricT),0));
[GradientMeltFraction_x,GradientMeltFraction_y,GradientMeltFraction_z]=gradient(Box.SymmetricMeltFraction,Vector.x.*1000,Vector.y.*1000,Vector.z.*1000);
Box.SymmetricMeltProduction=max(-GradientMeltFraction_z.*Box.SymmetricVelocity_z,0);
for iRow=1:size(Box.x,1);
    for iCol=1:size(Box.x,2);
        Box.SymmetricMeltFlux(iRow,iCol)=sum(reshape(Box.SymmetricMeltProduction(iRow,iCol,:),size(Vector.z)).*gradient(Vector.z).*1000);
    end
end

Box.MeltFraction=(max(meltFunction(Box.z,Box.T),0));
[GradientMeltFraction_x,GradientMeltFraction_y,GradientMeltFraction_z]=gradient(Box.MeltFraction,Vector.x.*1000,Vector.y.*1000,Vector.z.*1000);
Box.MeltProduction=max(-GradientMeltFraction_z.*Box.Velocity_z,0);
for iRow=1:size(Box.x,1);
    for iCol=1:size(Box.x,2);
        Box.MeltFlux(iRow,iCol)=sum(reshape(Box.MeltProduction(iRow,iCol,:),size(Vector.z)).*gradient(Vector.z).*1000);
    end
end

Box.DiffT=Box.T-Box.SymmetricT;
Box.DiffVelocity=Box.Velocity-Box.SymmetricVelocity;
Box.DiffVelocity_x=Box.Velocity_x-Box.SymmetricVelocity_x;
Box.DiffVelocity_y=Box.Velocity_y-Box.SymmetricVelocity_y;
Box.DiffVelocity_z=Box.Velocity_z-Box.SymmetricVelocity_z;
Box.DiffMeltFraction=Box.MeltFraction-Box.SymmetricMeltFraction;
Box.DiffMeltProduction=Box.MeltProduction-Box.SymmetricMeltProduction;
Box.DiffMeltFlux=Box.MeltFlux-Box.SymmetricMeltFlux;

%% Comparison Data Interpolation

for iRow=1:size(Box.x,1);
    for iCol=1:size(Box.x,2);
        DiffT(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffT(iRow,iCol,:),size(Box.x,3),1),Depth);
        DiffVelocity(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffVelocity(iRow,iCol,:),size(Box.x,3),1),Depth);
        DiffVelocity_x(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffVelocity_x(iRow,iCol,:),size(Box.x,3),1),Depth);
        DiffVelocity_y(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffVelocity_y(iRow,iCol,:),size(Box.x,3),1),Depth);
        DiffVelocity_z(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffVelocity_z(iRow,iCol,:),size(Box.x,3),1),Depth);
        DiffMeltFraction(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffMeltFraction(iRow,iCol,:),size(Box.x,3),1),Depth);
        DiffMeltProduction(iRow,iCol)=interp1(reshape(Box.z(iRow,iCol,:),size(Box.x,3),1),reshape(Box.DiffMeltProduction(iRow,iCol,:),size(Box.x,3),1),Depth);
    end
end

%% Comparison Results Plot

% % melting results plot
figure(indFigure); clf;

subplot(2,2,1); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffT); shading interp; colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Temperature Difference [degC]');

subplot(2,2,2); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),Box.DiffMeltFlux); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Melt Flux Difference [m/s]');

subplot(2,2,3); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffMeltFraction*100); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
% set(gca,'Clim',[-10,10]);
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Melt Fraction Difference [%]');

subplot(2,2,4); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffMeltProduction); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
% set(gca,'Clim',[-2e-10,2e-10]);
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Melt Production Difference [1/s]');

% % velocity results plot
figure(indFigure+1); clf;

subplot(2,2,1); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffVelocity); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Velocity Difference [m/s]');

subplot(2,2,2); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffVelocity_x); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Velocity_x Difference [m/s]');

subplot(2,2,3); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffVelocity_y); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Velocity_y Difference [m/s]');

subplot(2,2,4); hold on;
pcolor(Box.x(:,:,1),Box.y(:,:,1),DiffVelocity_z); shading interp;  colorbar;
plot(Geometry.PlateBoundary.x,Geometry.PlateBoundary.y,'k.-');
axis equal; axis([Geometry.ModelBoundary.x,Geometry.ModelBoundary.y]);
title('Velocity_z Difference [m/s]');

return