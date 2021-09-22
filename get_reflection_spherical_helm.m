function [vspec, vtrans, grazAng, phi1] = get_reflection_spherical_helm (H, e, R, Ht)

% A. Helm (2008)
% Ground-based GPS altimetry with the L1 OpenGPS receiver
% using carrier phase-delay observations of reflected GPS signals
% Scientific Technical Report 08/10
% doi: 10.2312/GFZ.b103-08104
% www.gfz-potsdam.de/bib/pub/str0810/0810.pdf

% Solution to spherical reflection with Helm (2008) equations
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
    Ht = 20e6; %in meters
end

%% Location vector of transmitter/satellite
vtrans = get_satellite_position (H,e,R,Ht,1);

%% Location vector of receiver
vrec = [0;R+H]; % Location vector of receiver

%% Solution to root of quartic polynomial (gamma)
% Solution to roots of quartic polynomial iterating gamma (without using all of roots)
% Iterative scheme based on modified Newthon-method - p.35-36
[c0,c1,c2,c3,c4,~,gamma0] = quartic_param (vrec, vtrans, R); % Quartic coefficients
    
i = 1;
i_max = 1000;
tns(1) = tand(gamma0/2);
gammas(1) = gamma0;
% cond = 1e-8*(180/pi);
cond = sqrt(eps());

while (i <= i_max)

    tn = tns(i);
    fn = c4.*tn.^4 + c3.*tn.^3 + c2.*tn.^2 + c1.*tn + c0;
    f1n = 4.*c4.*tn.^3 + 3.*c3.*tn.^2 + 2.*c2.*tn + c1;
    f2n = 12.*c4.*tn.^2 + 6.*c3.*tn + 2.*c2;
    kn = abs(f2n./f1n);
    tns(i+1) = tn-(fn./f1n);
    gammas(i+1) = 2.*atand((tns(i+1)));

    if (abs((gammas(i+1))-(gammas(i))) < cond)
        break
    end

    if (f1n==0 && kn>1) == true
        break
    end
    i = i+1;
end

if (i > i_max)
    disp('Did not converged!');
end

gamma = gammas(i);

%% Position vector of specular reflection point 
vspec = [(R*cosd(gamma));(R*sind(gamma)-R)];

%% Change origin of system from geocenter to subreceiver
vrec(2) = vrec(2)-R;
vtrans(2) = vtrans(2)-R;

%% Grazing angle
grazAng = 90-(gamma-e);

%% Geocentric angle between reflection point and subreceiver (phi1)
phi1 = get_geocentric_angle (H,Ht,e,grazAng,R);
