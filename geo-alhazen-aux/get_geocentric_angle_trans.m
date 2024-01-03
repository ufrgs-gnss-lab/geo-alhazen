function [geo_ang_at, D] = get_geocentric_angle_trans (Ha,e,Ht,Rs)

% Returns geocentric angle between antenna and transmitter
%
% INPUT:
% - Ha: Antenna height (in meters)
% - Ht: Transmitter/satellite height (in meters)
% - e: Elevation angle of transmitter (in degrees)
% - Rs: Earth surface radius (in meters)

if nargin==2 || isempty(Ht) || Ht == 0
    Ht = get_satellite_height(); %in meters - p. 335
end

if nargin==3 || isempty(Rs) || Rs == 0
    Rs = get_earth_radius(); %in meters
end

%% Radii:
Ra = Rs+Ha;
Rt = Rs+Ht;

%% Normalize radii to values closer to 1 
% (to avoid overflow when squaring, R^2):
R0 = Rs;
Ra = Ra./R0;
Rt = Rt./R0;
Rs = Rs./R0;

%% Geocentric Angle Antenna-Transmitter
% Coefficients of quadratic polynomial
a=1;
b=-2*(Ra).*cosd(90+e);
c=(Ra).^2-(Rt).^2;

% Delta of quadratic
delta=b.^2-4.*a.*c;

% Roots
r1=(-b+delta.^0.5)./(2*a);
r2=(-b-delta.^0.5)./(2*a);

% Choose of correct root
D=max(r1,r2);

% Geocentric angle
geo_ang_at=acosd((D.^2 -(Ra).^2 -(Rt).^2)./(-2*(Ra).*(Rt)));
geo_ang_at=geo_ang_at';
%% Denormalization
D = D*R0;