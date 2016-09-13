function Line=lineStore(LinePosition,LineDis,Grid,Lid,LineStart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lineStore.m
% Store current melt trajectory informaiton
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   LinePosition            : Melt trajectory coordinates [km]
%   LineDis                 : Distance along melt trajectory [km]
%   Lid
%       |.Depth             : Depth of lid [km]
%       |.DepthGrad         : Slope of lid in x or y direction
%       |.T                 : Temperature of lid [degC]
%       |.T_Solidus         : Solidus temperature at lid [degC]
%       |......
%   LineStart               : Coordinates of line starting point (seed) [km]
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
% Line
%       |.x, .y,            : Melt trajectory coordinates [km]
%       |.Distance          : Distance along trajectory [km]
%       |.Depth             : Melt trajectory depth [km]
%       |.DepthGrad         : Slope of trajectory in x or y direction
%       |.T                 : Temperature along trajectory [degC]
%       |.T_Solidus         : Solidus temperature along trajectory [degC]
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
% 	None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Line Information Storage

Line.x=LinePosition(:,1);
Line.y=LinePosition(:,2);
Line.Distance=LineDis;
Line.Depth=interp2(Grid.x,Grid.y,Lid.Depth,LinePosition(:,1),LinePosition(:,2));

Line.DepthGrad_x=interp2(Grid.x,Grid.y,Lid.DepthGrad_x,Line.x,Line.y);
Line.DepthGrad_y=interp2(Grid.x,Grid.y,Lid.DepthGrad_y,Line.x,Line.y);
Line.Slope=((Line.DepthGrad_x).^2+(Line.DepthGrad_y).^2).^(1/2);
Line.T=interp2(Grid.x,Grid.y,Lid.T,Line.x,Line.y);
Line.T_Solidus=interp2(Grid.x,Grid.y,Lid.T_Solidus,Line.x,Line.y);
Line.LineStart=LineStart;

return