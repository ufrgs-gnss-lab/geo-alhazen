function [graz_ang, geo_ang, pos_spec, pos_trans] = get_reflection_spherical_fujimura (e, Ha, Ht, R0, debugit)

% M. Fujimura et al. (2019)
% The Ptolemy–Alhazen Problem and Spherical Mirror Reflection
% Computational Methods and Function Theory  v.19, p. 135–155

% Solution to spherical reflection with Fujimura et al. (2019) equations
%
% INPUT:
% e: satellite elevation angle (in degree)
% Ha: antenna/receiver height (in meters)
% Ht: Transmitter/satelitte height (in meters)
% R0: Earth radius (in meters)
% debugit: debug flag - 0, 1, 2 (see code for details).
% 
% OUTPUT:
% pos_spec: reflection point location vector (x,y) (in meters), in local frame
% pos_trans: transmitter/satellite location vector (x,y) (in meters), in local frame
% graz_ang: grazing angle of spherical reflection that satisfies Snell's Law (in degrees)
% geo_ang: geocentric angle between receiver and reflection point (in degrees) 

if (nargin < 5),  debugit = [];  end
if isempty(debugit),  debugit = 0;  end

% Geocentric radii to receiving antenna and transmitting satellite
Ra = R0+Ha;
Rt = R0+Ht;

% Functions to convert from position vector to complex number and vice versa:
posgeo2complex = @(pos) complex(pos(2), pos(1))./R0;
complex2posgeo = @(cpx) [imag(cpx), real(cpx)].*R0;  % in geo frame

% Geocenter position in local frame:
%pos_geo = [0 R0];  % WRONG!
pos_geo = [0 -R0];

% Receiving antenna position:
pos_ant = [0 Ha];  % in local frame
pos_ant_geo = pos_ant - pos_geo;  % in geocentric frame
z1 = posgeo2complex(pos_ant_geo);  % as complex number
if (debugit==1)
    z1b = (Ra./R0);
    disp(z1b-z1)
end

% Transmitting satellite position:
[pos_trans, pos_trans_geo] = get_satellite_position (e, Ha, Ht, R0, 2);
z2 = posgeo2complex(pos_trans_geo);
switch debugit
case 1
    pos_transb = R0*[imag(z2), real(z2)-1];
    disp(pos_transb-pos_trans)
    pos_trans = pos_transb;
case 2
    % Geocentric angle between antenna and satellite
    geo_ang_trans = get_geocentric_angle_trans (e, Ha, Ht, R0);
    z2c = (Rt./R0)*exp(1i*deg2rad(geo_ang_trans)); % as complex number
    pos_trans_geoc = complex2posgeo(z2c);  % in geocentric frame
    pos_transc = pos_trans_geoc + pos_geo;  % in local frame
     %disp(pos_transc-pos_trans)
     %disp(pos_trans_geoc-pos_trans_geo)
     %disp(z2c-z2)
    
    pos_trans = pos_transc;
    pos_trans_geo = pos_trans_geoc; %#ok<NASGU>
    z2 = z2c;
end

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
if (debugit==1)
    pos_specb = R0*[imag(u), real(u)-1];
    disp(pos_specb-pos_spec)
end

%% Geocentric angle between antenna and reflection point
geo_ang = rad2deg(angle(u));

%% Grazing Angle
graz_ang = atand(cotd(geo_ang)-(R0./Ra)./sind(geo_ang));

%%
if (debugit==1)
    graz_angb = get_grazing_angle_vector (pos_ant, pos_spec, pos_trans);
    geo_angb = get_geocentric_angle_spec (graz_angb, Ha, R0);
    disp(graz_angb-graz_ang)
    disp(geo_angb-geo_ang)
end

end

