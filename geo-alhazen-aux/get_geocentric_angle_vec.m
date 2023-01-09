function geo_ang = get_geocentric_angle_vec (pos_ant, pos)

% Returns geocentric angle of a target (vectorial formulation)
%
% OUTPUT:
% - geo_ang: Target geocentric angle
% 
% INPUT:
% - pos_ant: Antenna position (vector, in meters)
% - pos: Target position (vector, in meters)
% 
% Note: target may be surface or transmitter.

%% Geocentric Angle:
temp = dot_all(pos, pos_ant)./(norm_all(pos).*norm_all(pos_ant));
geo_ang = 2.*(acosd(temp)-90);

end