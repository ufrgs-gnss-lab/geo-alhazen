function [ev] = get_horizon_elevation_angle (Ha,R0)

% Get the elevation angle on the spherical horizon
%
% Input:
% Ha = receiver height (m)
% R0 = Earth radius (m)
%
% Output 
% ev = elevation angle on the horizon

if (nargin < 2) || isempty(R0),  R0 = get_earth_radius();  end

ev = asind(R0./(R0+Ha)) - 90;

% ev1 = acosd (R0./(R0+Ha));
