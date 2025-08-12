function [graz_ang, geo_ang_as, x_spec, y_spec, dx_trans, dy_trans, w0] = get_reflection_spherical_vuorinen (e, Ha, Rs)
% GET_REFLECTION_SPHERICAL_VUORINEN Calculates reflection on spherical 
% surface based on Vuorinen (2023) equations (internal document - filename cal20230511.tex)
%
% On the Ptolemy-Alhazen problem: source at infinite distance. 2023, May 11. Internal document. 
% 
% INPUT:
% - Ha: antenna/receiver height (in meters)
% - e: elevation angle (in radians)
% - Rs: Earth surface radius (in meters)
% 
% NOTE: transmitter/satellite height is assumed at infinite distance.
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

% Antenna geocentric rotated angle:
% (or satellite zenith angle):
phia = 90-e;

% Auxiliary variables:
fs = exp(1i.*deg2rad(phia));
fc = conj(fs);  % = exp(-1i.*deg2rad(thetapa));
fa = ra.*fs;

% Quartic polynomial coefficients
c4 = ra.*fc;
c3 = -1;
c2 = 0;
c1 = 1;
c0 = -ra.*fs;

% Polynomial roots
wk = roots ([c4 c3 c2 c1 c0]);

% Candidate geocentric rotated angles
phik = rad2deg(angle(wk));

% Reflection geocentric rotated angle:
phik (phik>90)=NaN;
ind = argmin(abs(phik-phia));
phi0 = phik(ind);
w0 = wk(ind);

% Geocentric angle at reflection point:
% (just remove the rotation)
geo_ang_as = phia-phi0;

% Reflection point in a quasigeocentric frame
phia_rad = deg2rad(phia);
pos_spec_complex = w0./exp(-1i*(pi./2-phia_rad)).*Rs-complex(0,Rs);
% x_spec = real(pos_spec_complex);
x_spec = phi0;
y_spec = imag(pos_spec_complex);
pos_spec_geo = [x_spec y_spec+Rs];

% Satellite direction (unit vector)
% in either local or quasigeocentric frames:
dx_trans = cosd(e);
dy_trans = sind(e);
dir_trans = [dx_trans, dy_trans];

% Grazing angle
graz_ang = 90-rad2deg(angle(fa./w0-1));

end 