function [delay, graz_ang, arc_len, slant_dist, pos_spec, pos_trans, elev_spec, delay_trig]= get_reflection_spherical(e, Ha, Ht, R0, algorithm, trajectory)

% Return results of reflection on spherical surface.
%
% INPUT:
% - es: elevation angles nnnr height (scalar, in meters)
% - Ht: Transmitter/satelitte height (scalar, in meters)
% - R0: Earth radius (scalar, in meters)
% - algorithm: (char), see source code for details
%   'fujimura'
%   'miller&vegh'
%   'martin-neira'
%   'helm'
%   'itu'
%   'fermat'
%   'line-sphere fermat'
%   'circle inversion'
% - trajectory: (char), see source code for details
%   'orbital'  % orbital trajectory (constant geocentric distance)
%   'horizontal'  % horizontal trajectory (constant y-axis coordinate)
%   'circular'  % circular trajectory around the receiving antenna (constant direct distance)
%
% OUTPUT:
% - delay: interferometric propagation delay on spherical surface (vector, in meters)
% - graz_ang: grazing angle of spherical reflection that satisfies Snell's law (vector, in degrees)
% - arc_len: arc lenght between subreceiver and reflection point (vector, in meters)
% - slant_dist: slant distance between receiver and reflection point (vector, in meters)
% - pos_spec: reflection point location vector (x,y) (matrix, in meters)
% - pos_trans: transmitter point location vector (x,y) (matrix, in meters)
% - elev_spec: elevation angle of spherical reflection (vector, in degrees)
% - delay_trig: trigonometric formulation of interferometric propagation delay on spherical surface (vector, in meters)

if (nargin < 1),  e = [];  end
if (nargin < 2),  Ha = [];  end
if (nargin < 3),  Ht = [];  end
if (nargin < 4),  R0 = [];  end
if (nargin < 5),  algorithm = [];  end
if (nargin < 6),  trajectory = [];  end

if isempty(e),   e  = 45;  end
if isempty(Ha),  Ha = 10;  end
if isempty(Ht),  Ht = get_satellite_height();  end
if isempty(R0),  R0 = get_earth_radius();  end
if isempty(algorithm),  algorithm = 'fujimura';  end
if isempty(trajectory),  trajectory = 'orbital';  end

siz = size(e);  e = e(:);  n = prod(siz);

%% transmitting satellite's trajectory:
pos_ant  = [0 Ha];
pos_foot = [0 0];
pos_geo  = [0 -R0];
pos_foot_geo = pos_foot - pos_geo;
switch lower(trajectory)
    case 'orbital'  % orbital trajectory (constant geocentric distance)
        Hts = repmat(Ht, [n 1]);
    case 'horizontal'  % horizontal trajectory (constant y-axis coordinate)
        vert_dist = (Ht-Ha)*ones(n,1);
        horiz_dist = vert_dist.*cotd(e(:));
        pos_trans_ant = [horiz_dist vert_dist];
        pos_trans = pos_ant + pos_trans_ant;
        pos_trans_geo = pos_trans + pos_foot_geo;
        Hts = norm_all(pos_trans_geo);
    case 'circular'  % circular trajectory around the receiving antenna (constant direct distance)
        %pos_trans = Ht.*[cosd(e), sind(e)];  % WRONG!
        direct_dist = Ht - Ha;
        pos_trans_ant = direct_dist.*[cosd(e(:)), sind(e(:))];
        pos_trans = pos_ant + pos_trans_ant;
        pos_trans_geo = pos_trans + pos_foot_geo;
        Hts = norm_all(pos_trans_geo);
    otherwise
        error('Unknown transmitting satellite trajectory "%s".', char(transmitter))
end

%% Choosing algorithm
switch lower(algorithm)
    case {'fujimura'}
        f = @get_reflection_spherical_fujimura;
    case {'miller','miller&vegh'}
        f = @get_reflection_spherical_miller;
    case {'martinneira','martin-neira'}
        f = @get_reflection_spherical_martinneira;
    case 'helm'
        f = @get_reflection_spherical_helm;
     case {'itu','itu1008'}
        f = @get_reflection_spherical_itu1008;
     case {'fermat','fermat''s principle'}
        f = @get_reflection_spherical_fermat;
    case {'line-sphere fermat','linespherefermat','fermatlsi'}
        f = @get_reflection_spherical_virtual_fermat;
    case {'circle inversion','circleinversion','circlelsi','circleinversionlsi'}
        f = @get_reflection_spherical_virtual_circinv;
    otherwise
        error('Unknown algorithm "%s"', char(algorithm));
end

%% pre-allocate varibles:
graz_ang  = NaN(n,1);
geo_ang   = NaN(n,1);
pos_spec  = NaN(n,2);
pos_trans = NaN(n,2);

assert(isscalar(Ha));
ev = get_horizon_elevation_angle (Ha,R0);

%%
for i=1:n
    if (e(i) < ev),  continue;  end
    [graz_ang(i), geo_ang(i), pos_spec(i,:), pos_trans(i,:)] = f(e(i), Ha, Hts(i), R0);        
end

%% Additional results

% Arc Length from subreceiver point to reflection point:
arc_len = deg2rad(geo_ang)*R0;

% Receiving antenna position vector in local frame:
pos_ant = [0 Ha]; 

% Relative position vectors:
pos_trans_spec = pos_trans - pos_spec;
pos_ant_spec   = pos_ant   - pos_spec;
pos_ant_trans  = pos_ant   - pos_trans;

% Slant distance from reflection point to receiver:
slant_dist = norm_all(pos_ant_spec);

% Interferometric propagation delay:
delay_direct = norm_all(pos_ant_trans);
delay_reflect_inc = norm_all(pos_trans_spec);
delay_reflect_out = norm_all(pos_ant_spec);
delay_reflect = delay_reflect_inc + delay_reflect_out;
delay = delay_reflect - delay_direct;

% Reflection elevation angle:
dir_rec_spec = pos_ant_spec./delay_reflect_out;
dir_vert = [0 1];
elev_spec = -90 + acosd(dot_all(-dir_vert, -dir_rec_spec));
debugit = 0;
if debugit
  elev_specb = e - 2*graz_ang;
  elev_specb(e==90) = -90;
  fpk e erb-elev_spec  % DEBUG
end

% trigonometric interferometric delay:
delay_trig = slant_dist.*(1-cosd(e-elev_spec));
if debugit
  delay_trigb = slant_dist.*(1-cosd(e-elev_specb));
  delay_trigc = slant_dist.*(1-cosd(2*graz_ang));
  fpk e [delay_trigc delay_trigb]-delay_trig
end

%% reshape output vectors as in input vector:
delay = reshape(delay, siz);
graz_ang = reshape(graz_ang, siz);
%geo_ang = reshape(geo_ang, siz);
arc_len = reshape(arc_len, siz);
slant_dist = reshape(slant_dist, siz);
elev_spec = reshape(elev_spec, siz);
delay_trig = reshape(delay_trig, siz);

%%

end

