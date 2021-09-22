function [vspec, vtrans, grazAng, phi1] = get_reflection_spherical_martinneira (H, e, R, Ht)

% M. Martin-Neira (1993)
% A Passive Reflectometry and Interferometry System
% (PARIS): Application to Ocean Altimetry
% ESA Journal, v. 17

% Solution to spherical reflection with Martin-Neira (1993) equations
%
% INPUT:
% H: antenna/receiver height (in meters)
% e: elevation angle (in degree)
% R: Earth radius (in meters)
% Ht: Transmitter/satelitte height (in meters)
% 
% OUTPUT:
% vspec: reflection point location vector (x,y) (in meters)
% vtrans: transmitter/satellite location vector (x,y) (in meters)
% grazAng: grazing angle of spherical reflection that satisfies Snell's Law (in degree)
% phi1: Geocentric angle between receiver and reflection point (in radians) 

%"v" before variables consists in a vector

if nargin==2 || isempty(R) || R == 0
    R = 6370e3; %in meters
end
 
% Average radius of GPS orbit is 26000 km (20000km + Earth radius)
if nargin==3 || isempty(Ht) || Ht == 0
    Ht = 20e6; %in meters - p. 335
end

%% Location vector of transmitter/satellite
vtrans = get_satellite_position (H,e,R,Ht,1);

%% Location vector of receiver
vrec = [0;R+H]; % Location vector of receiver

%% Solution to root of quartic polynomial (gamma)
[c0,c1,c2,c3,c4,t0] = quartic_param (vrec, vtrans, R); % Quartic coefficients and start value

cs = [c4 c3 c2 c1 c0];
ts = roots(cs); % Roots of quartic
[~,ind] = min(abs(ts-t0));
t = ts(ind);
gamma = atand(t).*2; % Rotation angle phi, p. 333

%% Position vector of specular reflection point 
vspec = [(R*cosd(gamma));(R*sind(gamma)-R)];

%% Change origin of system from subreceiver
vrec(2) = vrec(2)-R;
vtrans(2) = vtrans(2)-R;

%% Grazing angle
grazAng = 90-(gamma-e);

%% Geocentric angle between reflection point and subreceiver (phi1)
phi1 = get_geocentric_angle (H,Ht,e,grazAng,R);
