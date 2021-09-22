function [g,gl,gr] = get_grazing_angle_vector (pos_ant,pos_sfc,pos_sat,pos_geocenter)

% Return grazing angle calculated by vector relation
% 
% INPUT:
% - pos_ant: position vector of antenna (x,y)
% - pos_sfc: position vector of reflection point (x,y)
% - pos_sat: position vector of antenna (x,y)
% - pos_geocenter: position vector of geocenter (x,y), optional
% 
% OUTPUT:
% - g: grazing angle (in degrees)
% - gl: left grazing angle (in degrees)
% - gr: right grazing angle (in degrees)

% relative position vectors, from surface to antenna or satellite:
pos_sat_sfc = pos_sat - pos_sfc;
pos_ant_sfc = pos_ant - pos_sfc;

% unit direction vectors, from surface to antenna or satellite:
dir_sat_sfc = pos_sat_sfc ./ norm_all(pos_sat_sfc);
dir_ant_sfc = pos_ant_sfc ./ norm_all(pos_ant_sfc);

% bistatic scattering angle at the surface, between antenna and satellite:
rst = acosd(dot_all(dir_ant_sfc, dir_sat_sfc));

% grazing angle:
g = 90 - 0.5*rst;

if (nargout < 2),  return;  end

% vertical direction at the surface reflection point:
pos_sfc_geocenter = pos_sfc - pos_geocenter;
dir_vert = pos_sfc_geocenter ./ norm_all(pos_sfc_geocenter);

% left and right grazing angles
gl = 90 - acosd(dot_all(dir_ant_sfc, dir_vert));
gr = 90 - acosd(dot_all(dir_vert, dir_sat_sfc));

end

