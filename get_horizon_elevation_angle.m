function [emin] = get_horizon_elevation_angle (Ha,Rs)

% Get the elevation angle on the spherical horizon
%
% Input:
% Ha = receiver height (m)
% Rs = Earth radius (m)
%
% Output 
% ev = elevation angle on the horizon

if (nargin < 2) || isempty(Rs),  Rs = get_earth_radius();  end

emin = asind(Rs./(Rs+Ha)) - 90;

