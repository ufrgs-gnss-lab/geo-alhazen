% function [delay, vspec, grazAng, slantDist, arcLen] = get_reflection_spherical_miller (H, e, R, Rt, de)
function [vspec, vtrans, grazAng, phi1] = get_reflection_spherical_miller (H, e, R, Ht)
% Main reference 
% A. R. Miller and E. Vegh (1993)
% "Exact Result for the Grazing Angle of Specular Reflection from a Sphere"
% SIAM Review, Vol. 35, No. 3 (Sep., 1993), pp. 472-480
% <https://doi.org/10.1137/1035091>
% <https://www.jstor.org/stable/2132427>

% Complementary reference
% A. R. Miller and E. Vegh (1990)
% "Computing the grazing angle of specular reflection"
% International Journal of Mathematical Education in Science and
% Technology, v. 21, n. 2
% <https://doi.org/10.1080/0020739900210213>

% Solution to spherical reflection with Miller and Vegh (1993) and 
% Miller and Vegh (1990) and equations
%
% INPUT:
% H: antenna/receiver height (in meters)
% e: elevation angle (in radians)
% R: Earth radius (in meters)
% Ht: Transmitter/satelitte height (in meters)
% 
% OUTPUT:
% vspec: reflection point location vector (x,y) (in meters)
% vtrans: transmitter/satellite location vector (x,y) (in meters)
% grazAng: grazing angle of spherical reflection that satisfies Snell's Law (in degree)
% phi1: Geocentric angle between receiver and reflection point (in radians) 

%"v" before variables consists in a vector
%%
if nargin==2 || isempty(R) || R == 0
    R = 6370e3; %in meters
end

% Average radius of GPS orbit is 20000 km
if nargin==3 || isempty(Ht) || Ht == 0
    Ht = 20e6; %in meters - p. 335
end

h1=H;  % antenna height
h2=Ht;  % satellite height
e = deg2rad(e); % elevation angle should be in radians
%% Quartic poynomial solution
% Slant distance between antenna and satellite (NOT between antenna and surface)
% (law of cosines at the antenna with nadir angle to satellite 90deg+e)
% (R+h2).^2 = d.^2 +(R+h1).^2 -2.*d.*(R+h1).*cos(pi/2+e)
% (quadratic polynomical, roots via Bhaskara's formula)
% d.^2 +(R+h1).^2 -2.*d.*(R+h1).*cos(pi/2+e) -(R+h2).^2=0
a=1;
b=-2*(R+h1).*cos(pi/2+e);
c=(R+h1).^2-(R+h2).^2;
D=b.^2-4.*a.*c;
d1=(-b+D.^0.5)./(2*a);
d2=(-b-D.^0.5)./(2*a);
d=max(d1,d2);
clear a b c D

% spheric-centric angle (phi) - Fig.1, p.473
% (law of cosines at the center of sphere)
phi=acos((d.^2 -(R+h1).^2 -(R+h2).^2)./(-2*(R+h1).*(R+h2)));

k1=R./(R+h1);  % mid p.473
k2=R./(R+h2);  % mid p.473

alpha=exp(1i.*phi).*(exp(1i.*phi) -k1.*k2);  % top p.473
beta=k1.^2 +k2.^2 -2.*k1.*k2.*exp(1i.*phi);  % top p.473
gamma=2*(k1.^2 +k2.^2 -k1.*k2.*cos(phi) -1);  % top p.473

n = numel(e);
zs = NaN(n,4);

for i=1:numel(e)
    % Roots of quartic polynomial
    zs(i,:) = roots([alpha(i), beta(i), gamma(i), conj(beta(i)), conj(alpha(i))]);  % top p.473
end

psi = z2psi(zs,phi, k1, k2);
grazAng = psi.*180./pi; % Grazing angle of spherical specular reflection

%% Geocentric angles
% Geocentric angle between receiver and reflection point
phi1 = get_geocentric_angle (H,Ht,rad2deg(e),grazAng,R);

%% Reflection point location vector
x = R.*sin(phi1);   % X coordinate of Reflection point
y = R.*cos(phi1)-R; % Y coordinate of Reflection point
vspec = [x y];

%% Location vector of transmitter/satellite
vtrans = get_satellite_position (H,rad2deg(e),R,Ht,0);

end
