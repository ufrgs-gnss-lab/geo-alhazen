function [pos_sfc, pos_sat, grazAng, phi1] = get_reflection_spherical_fermat (H,e,R,Ht)

%%
if nargin==2 || isempty(R) || R == 0
    R = 6370e3; %in meters
end

% Average radius of GPS orbit is 20000 km
if nargin==3 || isempty(Ht) || Ht == 0
    Ht = 20e6;
end

%% Interferometric delay 
% get distance between any two points (or pairs of points):
sumcol = @(arg) sum(arg, 2);  % sum columns
get_dist = @(pos1, pos2) sumcol((pos1-pos2).^2).^0.5;

% get reflection distance:
get_dist_reflect = @(pos_ant, pos_sat, pos_sfc) ...
  get_dist(pos_ant, pos_sfc) + get_dist(pos_sfc, pos_sat);

% get direct distance:
get_dist_direct = @(pos_ant, pos_sat) ...
  get_dist(pos_ant, pos_sat);

% get interferometric distance:
get_dist_interf = @(pos_sfc, pos_ant, pos_sat) ...
  get_dist_reflect(pos_ant, pos_sat, pos_sfc) - get_dist_direct(pos_ant, pos_sat);


%%
% get slant distance:
get_dist_slant_aux = @(e, R1, h21)...
  sqrt((R1 + h21)^2 - R1.^2.*cosd(e).^2) - R1.*sind(e);
get_dist_slant = @(e, R0, h1, h2) ...
  get_dist_slant_aux(e, R0+h1, h2-h1);

% get antenna position:
get_pos_ant = @(h_ant) [0 h_ant];

% get satellite position:
get_pos_sat = @(elev_sat, R0, h_ant, h_sat) ...
    [h_sat.*(cosd((elev_sat)+asind(((sind(90+(elev_sat))).*(R0+h_ant)./h_sat))))
     h_sat.*(sind((elev_sat)+asind(((sind(90+(elev_sat))).*(R0+h_ant)./h_sat))))];
 
%%
% get sphere surface y coordinate 
% (origin is at antenna's foot -- NOT at the antenna!)
% (center of sphere is at x0=0, y0=-R0)
% formula of eccentric sphere is: (x-x0)^2+(y-y0)^2=R0^2
get_y_sph = @(x, R0) sqrt(R0.^2 -x.^2) - R0;

% get sphere surface position:
get_pos_sfc = @(x, R0) [x get_y_sph(x, R0)];

%%
R0 = R;
h_ant = H;
h_sat = Ht;
pos_ant = get_pos_ant(h_ant);

% Satellite position
pos_sat = get_satellite_position (h_ant,e,R0,h_sat,0);

opt = struct('MaxIter',1000, 'TolX',1e-6, 'TolFun',1e-6);
f = @(x) get_dist_reflect(pos_ant, pos_sat, get_pos_sfc(x, R0));
x_pla = h_ant./tand(e); % Reflection point x axis on plane

if e~=0
    [x_sph] = fminsearch(f, x_pla,opt); % Reflection point x axis on sphere by Fermat's Principle
elseif e==0
    [x_sph,~,exitflag] = fminsearch(f, 0, opt);
end

% if exitflag~=1
%     fprintf ('Did not converged')
% end

pos_sfc = get_pos_sfc(x_sph, R0); % Reflection point location

% Grazing Angle
grazAng = real(get_grazing_angle_vector (pos_ant,pos_sfc,pos_sat));

% Geocentric angle
phi1 = get_geocentric_angle (H,Ht,e,grazAng,R0);