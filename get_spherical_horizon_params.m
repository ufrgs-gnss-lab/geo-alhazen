function [delay, graz_ang, arc_len, slant_dist, pos_spec, ehor] = get_spherical_horizon_params (Ha, Rs)

% Returns parameters of the specular reflection on the spherical horizon.
% All parameters are precisely derived from closed trigonometric 
% formulations or exact expected results. 
% 
% Input:
% Ha = antenna height above the surface (in meters)
% R0 = radius of the sphere (in meters)
% 
% Output:
% delay = interferometric delay (in meters)
% graz_ang = grazing angle (degrees)
% arc_len = arc length (in meters)
% slant_dist = slant distance (in meters)
% pos_spec = reflection point [Xhor, Yhor] (in meters)
% ehor = elevation angle (degrees)


if (nargin < 2) || isempty(Rs),  Rs = get_earth_radius();  end

%%
ehor = get_horizon_elevation_angle(Ha,Rs);
slant_dist = sqrt(2.*Rs.*Ha+Ha.^2);    
Yspec = Rs-Ha./(1+Ha./Rs);
Xspec = sqrt(Rs.^2-Yspec.^2);
pos_spec = [Xspec' (Yspec-Rs)'];
arc_len = Rs.*(asin(Xspec./Rs));
graz_ang = zeros (size(Ha));
delay = zeros (size(Ha));

end