function geo_ang = get_geocentric_angle_vec (pos_ant, pos)

% Returns geocentric angle of a target (vectorial formulation)
%
% OUTPUT:
% - geo_ang: geocentric angle between antenna and target (vector, in degrees)
% 
% INPUT:
% - pos_ant: Antenna geocentric position (vector, in meters)
% - pos: Target geocentric position (vector, in meters)
% 
% Note: target may be surface or transmitter.

%% Geocentric Angle:
temp = dot_all(pos, pos_ant)./(norm_all(pos).*norm_all(pos_ant));
geo_ang = acosd(temp);

end