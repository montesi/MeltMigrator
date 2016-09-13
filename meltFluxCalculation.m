function Shoulder=meltFluxCalculation(Shoulder,meltFunction,MeltFractionCutoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meltFluxCalculation.m
% Calculate melt flux in melt shoulders and swaths
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Shoulder(ind)
%       |.Swath(ind)
%           |.Column_z      : Depth of vertical sampling column below polygons [km]
%           |.Column_T      : Temperature at sampling columns [degC]
%           |.Column_VerticalVelocity : Vertical velocity at sampling columns [m/s]
%           |......
%       |......
%   MeltFunction            : Melt fraction function
%   MeltFractionCutoff      : Retained melt fraction 
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   Shoulder(ind)
%       |.Swath(ind)
%           |.MeltFraction  : Melt fraction at sampling columns
%           |.MaxMeltFraction : Maximum melt fraction at sampling columns
%           |.AverageMeltFraction : Average melt fraction at sampling columns 
%           |.MeltProduction : Melt production at sampling columns [1/s]
%           |.MeltFlux      : Melt flux at sampling columns [m/s]
%           |.MeltingDepth  : Melting depth at sampling columns [km]
%           |......
%       |......
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   iShoulder               : Shoulder index
%   iSwath                  : Swath index
%   iPolygon                : Polygon index
%   iDepth                  : Sampling column depth index
%   LocalSwath              : Local swath information, short for Shoulder(iShoulder).Swath(iSwath)
%   MeltFraction
%   MaxMeltFraction
%   AverageMeltFraction
%   VerticalMeltFractionGradient
%   HorizontalMeltFractionGradient
%   MeltProduction
%   MeltFlux
%   MeltingDepth
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

warning off all
for iShoulder=1:numel(Shoulder);
    for iSwath=1:numel(Shoulder(iShoulder).Swath);
        LocalSwath=Shoulder(iShoulder).Swath(iSwath);
        MeltFraction=[]; 
        for iPolygon=1:size(LocalSwath.Column_z,1);
            for iDepth=1:size(LocalSwath.Column_z,2)
                MeltFraction(iPolygon,iDepth)=max(meltFunction(LocalSwath.Column_z(iPolygon,iDepth),LocalSwath.Column_T(iPolygon,iDepth)),0);
            end
        end
        MaxMeltFraction=max(MeltFraction,[],2);  % maximum melt fraction at each column
        [VerticalMeltFractionGradient,HorizontalMeltFractionGradient]=gradient(MeltFraction);
        MeltProduction=max(-VerticalMeltFractionGradient.*LocalSwath.Column_VerticalVelocity,0).*(MeltFraction>MeltFractionCutoff); % melt production rate
        MeltFlux=sum(MeltProduction,2); % melt flux
        
        AverageMeltFraction=sum(MeltProduction.*MeltFraction,2)./MeltFlux; % average melt fraction
        AverageMeltFraction(find(isnan(AverageMeltFraction)))=0;
        MeltingDepth=sum(MeltProduction.*LocalSwath.Column_z,2)./MeltFlux; % depth of melting
        MeltingDepth(find(isnan(MeltingDepth)))=0;

        % % store melting info
        Shoulder(iShoulder).Swath(iSwath).MeltFraction=MeltFraction;
        Shoulder(iShoulder).Swath(iSwath).MaxMeltFraction=MaxMeltFraction;
        Shoulder(iShoulder).Swath(iSwath).AverageMeltFraction=AverageMeltFraction;
        Shoulder(iShoulder).Swath(iSwath).MeltProduction=MeltProduction;
        Shoulder(iShoulder).Swath(iSwath).MeltFlux=MeltFlux;
        Shoulder(iShoulder).Swath(iSwath).MeltingDepth=MeltingDepth;
    end
end
warning on all;

return