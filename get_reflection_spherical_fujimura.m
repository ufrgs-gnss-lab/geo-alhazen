function [graz_ang, geo_ang, pos_spec, pos_trans] = get_reflection_spherical_fujimura (e, Ha, Ht, Rs)

% M. Fujimura et al. (2019)
% The Ptolemy–Alhazen Problem and Spherical Mirror Reflection
% Computational Methods and Function Theory,  v.19, p. 135–155

% Solution to spherical reflection with Fujimura et al. (2019) equations
%
% INPUT:
% e: satellite elevation angle (in degree)
% Ha: antenna/receiver height (in meters)
% Ht: Transmitter/satelitte height (in meters)
% Rs: Earth radius (in meters)
% debugit: debug flag - 0, 1, 2 (see code for details).
% 
% OUTPUT:
% pos_spec: reflection point location vector (x,y) (in meters), in local frame
% pos_trans: transmitter/satellite location vector (x,y) (in meters), in local frame
% graz_ang: grazing angle of spherical reflection that satisfies Snell's Law (in degrees)
% geo_ang: geocentric angle between receiver and reflection point (in degrees) 


% Geocentric radii to receiving antenna and transmitting satellite
Ra = Rs+Ha;
Rt = Rs+Ht;

% Functions to convert from position vector to complex number and vice versa:
posgeo2complex = @(pos) complex(pos(2), pos(1))./Rs;
complex2posgeo = @(cpx) [imag(cpx), real(cpx)].*Rs;  % in geo frame

% Geocenter position in local frame:
pos_geo = [0 -Rs];

% Receiving antenna position:
pos_ant = [0 Ha];  % in local frame
pos_ant_geo = pos_ant - pos_geo;  % in geocentric frame
z1 = posgeo2complex(pos_ant_geo);  % as complex number

% Transmitting satellite position:
[pos_trans, pos_trans_geo] = get_satellite_position (e, Ha, Ht, Rs, 2);
z2 = posgeo2complex(pos_trans_geo);

% Polynomial coefficients:
c4 = +conj(z1)*conj(z2);
c3 = -(conj(z1)+conj(z2));
c2 =  0;
c1 = (z1+z2);
c0 =-(z1*z2);
p = [c4 c3 c2 c1 c0];

% Polymomial roots solutions
us = roots(p);

% Fermat condition of minimum path lenght 
cond = abs(z1-us) + abs(z2-us);
[~,ind]= min(cond);

% Correct root
u = us(ind);

% Reflection point coordinates, in local frame:
pos_spec_geo = complex2posgeo(u);  % in geocentric frame
pos_spec = pos_spec_geo + pos_geo;  % in local frame

%% Geocentric angle between antenna and reflection point
geo_ang = rad2deg(angle(u));

%% Grazing Angle
graz_ang = atand(cotd(geo_ang)-(Rs./Ra)./sind(geo_ang));

end

