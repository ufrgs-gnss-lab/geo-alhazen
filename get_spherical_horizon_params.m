function [ehor,Dihor,Rspechor,ghor,sldisthor,arclenhor] = get_spherical_horizon_params (Ha, Rs)

% Returns parameters of the specular reflection on the spherical horizon.
% All parameters are precisely derived from closed trigonometric 
% formulations or exact expected results. 
% 
% Input:
% Ha = antenna height above the surface (in meters)
% R0 = radius of the sphere (in meters)
% 
% Output:
% ehor = elevation angle (degrees)
% sldisthor = slant distance (in meters)
% Rspechor = reflection point [Xhor, Yhor] (in meters)
% arclenhor = arc length (in meters)
% ghor = grazing angle (degrees)
% Dihor = interferometric delay (in meters)

if (nargin < 2) || isempty(Rs),  Rs = get_earth_radius();  end

%%
ehor = get_horizon_elevation_angle(Ha,Rs);
sldisthor = sqrt(2.*Rs.*Ha+Ha.^2);    
Yhor = Rs-Ha./(1+Ha./Rs);
Xhor = sqrt(Rs.^2-Yhor.^2);
Rspechor = [Xhor' (Yhor-Rs)'];
arclenhor = Rs.*(asin(Xhor./Rs));
ghor = zeros (size(Ha));
Dihor = zeros (size(Ha));

end