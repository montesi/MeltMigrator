function [TruncatedDis,TruncatedPosition]=lineTruncate(Distance,Position)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lineTruncate.m
% Truncate line to where P actually changes
% Laurent Montesi with Mark Behn, Laura Hebert
% Modified by Hailong Bai
% September 2015
%--------------------------------------------------------------------------
% INPUT -------------------------------------------------------------------
%   Distance                : Distance along line [km]
%   Position                : Line coordinates [km]
%--------------------------------------------------------------------------
% OUTPUT ------------------------------------------------------------------
%   TruncatedDis            : Truncated distance [km]
%   TruncatedPosition       : Truncated line coordinates [km]
%--------------------------------------------------------------------------
% INTERNAL ----------------------------------------------------------------
%   iTruncate               : Index for truncation
%--------------------------------------------------------------------------
% ATTENDING FUNCTIONS -----------------------------------------------------
%   None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

iTruncate=find((Position(:,1)~=Position(end,1))&(Position(:,2)~=Position(end,2)));
if isempty(iTruncate)
    iTruncate=1;
end
TruncatedPosition=Position(iTruncate,:);
TruncatedDis=Distance(iTruncate,:);

return