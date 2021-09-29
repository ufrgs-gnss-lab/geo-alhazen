function [graz_ang, geo_ang, pos_spec, pos_trans] = get_reflection_spherical_martinneira (e, Ha, Ht, Rs)

% M. Martin-Neira (1993)
% A Passive Reflectometry and Interferometry System
% (PARIS): Application to Ocean Altimetry
% ESA Journal, v. 17

% Solution to spherical reflection with Martin-Neira (1993) equations
%
% INPUT:
% Ha: antenna/receiver height (in meters)
% e: elevation angle (in radians)
% Ht: Transmitter/satelitte height (in meters)
% Rs: Earth radius (in meters)
% 
% OUTPUT:
% pos_spec: reflection point location vector (x,y) (in meters)
% pos_trans: transmitter/satellite location vector (x,y) (in meters)
% graz_ang: grazing angle of spherical reflection that satisfies Snell's Law (in degree)
% geo_ang: Geocentric angle between receiver and reflection point (in radians) 

%% quasigeocentric origin
pos_geo = [0 -Rs];

%% transmitter/satellite position
[pos_trans, pos_trans_geo]  = get_satellite_position (e,Ha,Ht,Rs,2);

%% antenna's position
pos_ant = [0 Ha]; % antenna's position in local frame
pos_ant_geo = pos_ant - pos_geo;

%% Solution to root of quartic polynomial (gamma)
[c0,c1,c2,c3,c4,t0] = quartic_param (pos_ant_geo, pos_trans_geo, Rs); % Quartic coefficients and start value

cs = [c4 c3 c2 c1 c0];
ts = roots(cs); % Roots of quartic
[~,ind] = min(abs(ts-t0));
t = ts(ind);
gamma = atand(t).*2; % Rotation angle phi, p. 333

%% specular reflection point position
pos_spec_geo = [Rs*cosd(gamma) Rs*sind(gamma)]; % quasigeo. frame
pos_spec = pos_spec_geo + pos_geo;              % local frame

%% Grazing angle
graz_ang = 90-(gamma-e);

%% Geocentric angle between reflection point and subreceiver (phi1)
geo_ang = get_geocentric_angle (Ha,Ht,e,graz_ang,Rs);