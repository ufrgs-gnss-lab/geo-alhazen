function [g,ga,gt] = get_grazing_angle_infinite (pos_ant, pos_spec, dir_sat, pos_geocenter)

% Return grazing angle calculated by vector relation
% 
% INPUT:
% - pos_ant: position vector of antenna (x,y)
% - pos_spec: position vector of reflection point (x,y)
% - pos_sat: direction vector of satellite (dx,dy)
% - pos_geocenter: position vector of geocenter (x,y), optional
% 
% OUTPUT:
% - g: grazing angle (in degrees)
% - ga: left grazing angle (in degrees)
% - gt: right grazing angle (in degrees)

% relative position vector, from antenna to surface:
pos_ant_sfc = pos_spec - pos_ant;

% distance from antenna to surface:
dist_ant_sfc = norm_all(pos_ant_sfc);

% unit direction vectors, from surface to antenna or satellite:
dir_ant_sfc = pos_ant_sfc ./ norm_all(pos_ant_sfc);

% bistatic scattering angle at the surface, between antenna and satellite:
rst = acosd(dot_all(dir_ant_sfc, dir_sat));

% grazing angle:
g = rst - 90;

if (nargout < 2),  return;  end

% vertical direction at the surface reflection point:
pos_sfc_geocenter = pos_spec - pos_geocenter;
dir_vert = pos_sfc_geocenter ./ norm_all(pos_sfc_geocenter);

% left and right grazing angles
ga = 90 - acosd(dot_all(dir_ant_sfc, dir_vert));
gt = 90 - acosd(dot_all(dir_vert, dir_sat_sfc));

end