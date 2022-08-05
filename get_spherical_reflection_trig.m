function [Di, sldist, g, gamma, der] = get_spherical_reflection_trig (e, Ha, Rs, varargin)

% GET_SPHERICAL_REFLECTION_TRIG_APPROX computes interferometric delay and
% slant distance on a spherical surface with trigonometric formulation
% based on an assumption of infinite distance of the transmitting satellite.
%
% INPUT:
% e: elevation angles (matrix; in degrees) 
% Ha: antenna height (matrix; in meters)
% Rs: Earth surface radius (scalar, in meters)
%
% OUTPUT:
% Di: interferometric delay or vaccum interferometric distance (matrix, in meters)
% Dl: slant distance between receiver and reflection point (matrix, in meters)

if isempty(Rs) || (nargin < 3);  Rs = get_earth_radius();  end

%%

[~, g, ~, ~, ~, ~, ~, ~, e_spec] = get_reflection_spherical (e, Ha, varargin);

Ra = Rs+Ha;
gamma = ((Ra.^2)/(Rs.^2))-cosd(g).^2;
sldist = Rs.*(sqrt(gamma)-sind(g));
der = e-e_spec; % or der = -(-e+e_spec)
Di = sldist.*(1-cosd(der));
% Di1 = sldist.*(cosd(2.*e+der)); % #WRONG!

end
