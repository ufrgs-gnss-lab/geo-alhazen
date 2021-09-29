function [graz_ang, geo_ang, pos_spec, pos_trans] = get_reflection_spherical_miller (e, Ha, Ht, Rs)
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

%% convert degree to radians
e = deg2rad(e); % elevation angle should be in radians

%% Quartic poynomial solution
% spheric-centric angle (phi) - Fig.1, p.473
% (law of cosines at the center of sphere)
[~,geo_ang_at] = get_geocentric_angle (Ha,Ht,rad2deg(e),0,Rs);

k1=Rs./(Rs+Ha);  % mid p.473
k2=Rs./(Rs+Ht);  % mid p.473

alpha=exp(1i.*geo_ang_at).*(exp(1i.*geo_ang_at) -k1.*k2);  % top p.473
beta=k1.^2 +k2.^2 -2.*k1.*k2.*exp(1i.*geo_ang_at);  % top p.473
gamma=2*(k1.^2 +k2.^2 -k1.*k2.*cos(geo_ang_at) -1);  % top p.473

n = numel(e);
zs = NaN(n,4);

for i=1:numel(e)
    % Roots of quartic polynomial
    zs(i,:) = roots([alpha(i), beta(i), gamma(i), conj(beta(i)), conj(alpha(i))]);  % top p.473
end

psi = z2psi(zs,geo_ang_at, k1, k2);
graz_ang = psi.*180./pi; % Grazing angle of spherical specular reflection

%% Geocentric angles
% Geocentric angle between receiver and reflection point
geo_ang = get_geocentric_angle (Ha,Ht,rad2deg(e),graz_ang,Rs);

%% Reflection point location vector
pos_spec_geo = [(Rs.*sind(geo_ang)) (Rs.*cosd(geo_ang))];
pos_geo = [0 -Rs];
pos_spec = pos_spec_geo + pos_geo;   

%% Location vector of transmitter/satellite
pos_trans = get_satellite_position (rad2deg(e),Ha,Ht,Rs,0);

end
