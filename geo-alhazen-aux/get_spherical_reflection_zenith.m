function [delay, graz_ang, arc_len, slant_dist, x_spec, y_spec] = get_spherical_reflection_zenith (Ha, Rs, frame)

% GET_SPHERICAL_REFLECTION_ZENITH: Returns parameters of the specular reflection on the spherical zenith.
% 
% All parameters are precisely derived from closed trigonometric formulations or exact expected results. 
% 
% INPUT:
% - Ha: antenna height above the surface (in meters)
% 
% OUTPUT:
% - delay: interferometric delay (in meters)
% - graz_ang: grazing angle (degrees)
% - arc_len: arc length (in meters)
% - slant_dist: slant distance (in meters)
% - x_spec, y_spec: specular reflection point coordinates (in meters)
% - ehor: elevation angle (degrees)
% 
% OPTIONAL INPUT:
% - R0: radius of the sphere (in meters)
% - frame: (char) coordinate reference frame ('local' - default - or 'quasigeo')

    if (nargin < 2) || isempty(Rs),  Rs = get_earth_radius();  end
    if (nargin < 3) || isempty(frame),  frame = 'local';  end

    delay = 2.*Ha;
    graz_ang = repmat(90,size(Ha));
    slant_dist = Ha;
    arc_len = zeros(size(Ha));
    y_spec = zeros(size(Ha));
    x_spec = zeros(size(Ha));
    
    if strcmpi(frame, 'local'),  return;  end
    % output will be in quasi-geocentric frame:
    y_spec = y_spec+Rs;
end