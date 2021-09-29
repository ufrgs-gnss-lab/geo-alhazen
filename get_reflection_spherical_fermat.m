function [graz_ang, geo_ang, pos_spec, pos_trans] = get_reflection_spherical_fermat (e, Ha, Ht, Rs)

% Numerical method based on Fermat's principle of least time to solve the
% reflection on a spher
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


%% Interferometric delay 
% get distance between any two points (or pairs of points):
sumcol = @(arg) sum(arg, 2);  % sum columns
get_dist = @(pos1, pos2) sumcol((pos1-pos2).^2).^0.5;

% get reflection distance:
get_dist_reflect = @(pos_ant, pos_sat, pos_spec) ...
  get_dist(pos_ant, pos_spec) + get_dist(pos_spec, pos_sat);

% get direct distance:
get_dist_direct = @(pos_ant, pos_sat) ... 
  get_dist(pos_ant, pos_sat); %Unused

% get interferometric distance:
get_dist_interf = @(pos_spec, pos_ant, pos_sat) ...
  get_dist_reflect(pos_ant, pos_sat, pos_spec) - get_dist_direct(pos_ant, pos_sat); %Unused

%%
% get slant distance:
get_dist_slant_aux = @(e, R1, h21)...
  sqrt((R1 + h21)^2 - R1.^2.*cosd(e).^2) - R1.*sind(e); %Unused
get_dist_slant = @(e, R0, h1, h2) ...
  get_dist_slant_aux(e, R0+h1, h2-h1); % Unused 

% get antenna position:
get_pos_ant = @(h_ant) [0 h_ant];
 
%% Reflection point
get_y_sph = @(x, R0) sqrt(R0.^2 -x.^2) - R0;

get_pos_spec = @(x, R0) [x get_y_sph(x, R0)];

%%
%Antenna position
pos_ant = get_pos_ant(Ha); %in local frame

% Satellite position
pos_trans = get_satellite_position (e,Ha,Ht,Rs,0); %local frame

%% Iterations
opt = struct('MaxIter',1000, 'TolX',1e-6, 'TolFun',1e-6);
f = @(x) get_dist_reflect(pos_ant, pos_trans, get_pos_spec(x, Rs));
x_pla = Ha./tand(e); % Reflection point x axis on plane

if e~=0
    [x_sph] = fminsearch(f, x_pla,opt); % Reflection point x axis on sphere by Fermat's Principle
elseif e==0
    [x_sph] = fminsearch(f, 0, opt);
end

%%
% Reflection point position
pos_spec = get_pos_spec(x_sph, Rs); % Reflection point location

% Grazing Angle
graz_ang = real(get_grazing_angle_vector (pos_ant,pos_spec,pos_trans));

% Geocentric angle
geo_ang = get_geocentric_angle (Ha,Ht,e,graz_ang,Rs);