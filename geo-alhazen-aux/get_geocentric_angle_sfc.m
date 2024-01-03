function geo_ang_as = get_geocentric_angle_sfc (Ha,g,Rs)

% Returns geocentric angle between antenna and reflection point
%
% OUTPUT:
% - geo_ang_sfc: Geocentric angle antenna-reflection point
% 
% INPUT:
% - Ha: Antenna height (in meters)
% - g: Grazing angle (in degrees)
% - Rs: Earth surface radius (in meters)

if nargin==2 || isempty(Rs) || Rs == 0
    Rs = get_earth_radius(); %in meters
end

%% Antenna radius:
Ra = Rs+Ha;

%% Geocentric Angle Antenna-Reflection Point
geo_ang_as = acosd((Rs./Ra).*cosd(g))-g;