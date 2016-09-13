function dP=crossToGrad(Position,Grid,Boundary,Lid,indDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crossToGrad.m;
% ODE to follow lid contour (perpendicular to slope)
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% June 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Position                : Current position coordinates
%   Grid             
%       |.x, .y             : 2D matrices of lid sampling coordinates
%   Boundary
%       |.x, .y, .z         : Model boundary in x, y, and z directions
%   Lid
%       |.Depth             : Depth of lid [km]
%       |.DepthGrad         : Slope of lid in x or y direction
%       |......
%   Cutoff                  : Cutoff criterion for flatness
%   indDir                  : Direction (1 for downslope, 2 for upslope)
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   dP                      : Unit vector; includes cutoff to define flat area
%--------------------------------------------------------------------------
% INTERNAL ---------------------------------------------------------------- 
%   C                       : Length of normalizing factor (0 if flatter than cutoff)
%   x, y, z                 : Current position
%   dzdx, dzdy              : Slope gradient of current position
%--------------------------------------------------------------------------
% ATTENDING SCRIPTS -------------------------------------------------------
%   inBounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

dP=zeros(2,1);
x=Position(1); y=Position(2); 
z=interp2(Grid.x,Grid.y,Lid.Depth,x,y);

InBoundaryIndex=inBounds(x,y,z,Boundary);
if InBoundaryIndex;
    dzdx=interp2(Grid.x,Grid.y,Lid.ExcessTGrad_x,x,y);
    dzdy=interp2(Grid.x,Grid.y,Lid.ExcessTGrad_y,x,y);
    
    C=(dzdx^2+dzdy^2).^(-1/2).*(-1)^(indDir);
    
    dP(1)=dzdy*C;
    dP(2)=-dzdx*C;
end
end