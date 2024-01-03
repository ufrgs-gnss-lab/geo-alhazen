function [graz_ang, geo_ang_as, x_spec, y_spec, x_trans, y_trans] = get_reflection_spherical_vuorinen (e, Ha, Ht, Rs)

% GET_REFLECTION_SPHERICAL_VUORINEN Calculates reflection on spherical 
% surface based on Vuorinen (2023) equations (internal document - filename cal20230511.tex)
%
% On the Ptolemy-Alhazen problem: source at infinite distance. 2023, May
% 11. Internal document. 
% 
% INPUT:
% - Ha: antenna/receiver height (in meters)
% - e: elevation angle (in radians)
% - Ht: Transmitter/satelitte height (in meters)
% - Rs: Earth radius (in meters)
% 
% OUTPUT:
% - x_spec, y_spec: reflection point in local frame (vectors, in meters)
% - x_trans, y_trans: transmitter point in local frame (vectors, in meters)
% - graz_ang: grazing angle of spherical reflection that satisfies Snell's Law (in degrees)
% - geo_ang_as: geocentric angle between antenna and reflection point (in degrees) 

% Antenna radius
Ra = Rs+Ha;

% Phase angle of the antenna
theta = deg2rad(90-e);

% Normalized antenna radius
ra = Ra./Rs;

% Polynomial coefficients
c4 = ra.*exp(-1i.*theta);
c3 = -1;
c2 = 0;
c1 = 1;
c0 = -ra.*exp(1i.*theta);

% Polynomial roots
ws = roots ([c4 c3 c2 c1 c0]);

% Phase angle of the reflection point
phis = angle(ws);

% Get correct root
ind = argmin(abs(phis-theta));
w = ws (ind);

% Reflection point on a quasigeocentric frame
pos_spec_complex = w./exp(-1i*(pi./2-theta)).*Rs-complex(0,Rs);
x_spec = real(pos_spec_complex);
y_spec = imag(pos_spec_complex);
pos_spec = [x_spec y_spec];

% Satellite position
[pos_trans, ~] = get_satellite_position (e, Ha, Ht, Rs, 2);
x_trans = pos_trans(1);
y_trans = pos_trans(2);

% Antenna position
pos_ant = [0 Ra];

% Grazing angle
graz_ang = get_grazing_angle_vector (pos_ant,pos_spec,pos_trans);

% Geocentric angle
geo_ang_as = deg2rad (theta-phis(ind));
