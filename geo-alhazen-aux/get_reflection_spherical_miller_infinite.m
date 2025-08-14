function [graz_ang, geo_ang_as, x_spec, y_spec, dx_trans, dy_trans, di, slant_dist] = get_reflection_spherical_miller_infinite (e, Ha, Rs)

% This algorithm is based on Miller, W.J., Barnes, J.W., MacKenzie, S.M., 2021. 
% Solving the Alhazen–Ptolemy Problem: Determining Specular Points on Spherical 
% Surfaces for Radiative Transfer of Titan’s Seas. Planet. Sci. J. 2, 63. 
% https://doi.org/10.3847/PSJ/abe4dd

to = deg2rad(e);
Ra = Rs+Ha;
c = Rs/Ra;
s_obs = sin(2 * to);
c_obs = cos(2 * to);
c2 = c^2 + 0j;
c2_4 = c2 - 4;
c2_43 = c2_4^3;
c4 = c^4;

D0 = 1152 * c2_4 * (-1 + c2 - c_obs + c2 * c_obs) + 864 * c2 * s_obs^2;
D1 = sqrt(-4 * (160 - 128 * c2 + 4 * c4 + 96 * c_obs - 96 * c2 * c_obs)^3 + (-16 * c2_43 - D0)^2);
D2 = (16 * c2_43 + D0 - D1)^(1/3);
D3 = (40 - 32 * c2 + c4 + 24 * c_obs - 24 * c2 * c_obs) / (3 * 2^(2/3) * D2) + D2 / (24 * 2^(1/3));
D4 = sqrt(-0.25 * c2_4 + (1/12) * c2_4 + D3);

p1 = real(acos((0.5 * D4 + 0.5 * sqrt(-(1/3) * c2_4 - D3 + (c * s_obs) / (2 * D4)))));
p2 = real(acos((0.5 * D4 - 0.5 * sqrt(-(1/3) * c2_4 - D3 + (c * s_obs) / (2 * D4)))));

% Use the slope of the line between (-pi/2, pi/2) to ({l}, 0) 
% to deduce the correct branch
l = (pi() * (atan(1 / c) - to)) / (pi() + 2 * atan(1 / c));

if abs(p1 - l) >= abs(p2 - l)
    ts = p1;
else
    ts = p2;
end

ts = rad2deg(ts);

% Specular position
x_spec = Rs*cosd(ts);
y_spec = Rs*sind(ts)-Rs;
pos_spec = [x_spec y_spec];
pos_spec_geo = [x_spec y_spec+Rs];

% Antenna position
pos_ant = [0 Ha];
pos_ant_geo = [0 Ra];

% Transmitter direction
dx_trans = cosd(e);
dy_trans = sind(e);
dir_trans = [dx_trans dy_trans];

% Interferometric delay and slant distance
[di, slant_dist] = get_delay_infinite_trans (pos_spec, pos_ant, dir_trans); 

% Geocentric angle between antenna and specular point
geo_ang_as = 90-ts;

% Grazing angle
t2 = 2.*ts-e;
graz_ang = 90-(t2-ts);


