function [graz_ang, geo_ang_as, x_spec, y_spec, dx_trans, dy_trans] = get_reflection_spherical_vuorinen (e, Ha, Rs)
% GET_REFLECTION_SPHERICAL_VUORINEN Calculates reflection on spherical 
% surface based on Vuorinen (2023) equations (internal document - filename cal20230511.tex)
%
% On the Ptolemy-Alhazen problem: source at infinite distance. 2023, May
% 11. Internal document. 
% 
% INPUT:
% - Ha: antenna/receiver height (in meters)
% - e: elevation angle (in radians)
% - Rs: Earth radius (in meters)
% 
% NOTE: transmitter/satellite height is assumed infinite.
% 
% OUTPUT:
% - x_spec, y_spec: reflection point in local frame (vector, in meters)
% - dx_trans, dy_trans: transmitter direction in local frame (unit vector, in meters)
% - graz_ang: grazing angle of spherical reflection that satisfies Snell's Law (in degrees)
% - geo_ang_as: geocentric angle between antenna and reflection point (in degrees) 

% Antenna radius
Ra = Rs+Ha;

% Normalized antenna radius
ra = Ra./Rs;

% Phase angle
theta = deg2rad(90-e);

% Quartic polynomial coefficients
c4 = ra.*exp(-1i.*theta);
c3 = -1;
c2 = 0;
c1 = 1;
c0 = -ra.*exp(1i.*theta);

% Polynomial roots
ws = roots ([c4 c3 c2 c1 c0]);

% Candidate phase angles
phis = angle(ws);

% Phase angle at reflection point:
% (TODO: test interf. delay, not phase angle.)
ind = argmin(abs(phis-theta));
phi = phis(ind);
w = ws (ind);

% Geocentric angle at reflection point:
geo_ang_as = deg2rad (theta-phi);

% Reflection point in a quasigeocentric frame
pos_spec_complex = w./exp(-1i*(pi./2-theta)).*Rs-complex(0,Rs);
x_spec = real(pos_spec_complex);
y_spec = imag(pos_spec_complex);
pos_spec_geo = [x_spec y_spec];

% Satellite direction (unit vector)
% in either local or quasigeocentric frames:
dx_trans = cosd(e);
dy_trans = sind(e);
dir_trans = [dx_trans, dy_trans];

% Antenna geocentric position
pos_ant_geo = [0 Ra]; 

% Grazing angle
graz_ang = get_grazing_angle_infinite (pos_ant_geo, pos_spec_geo, dir_trans);