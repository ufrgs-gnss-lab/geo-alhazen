function [graz_ang, geo_ang_as, x_spec, y_spec, dx_trans, dy_trans, di, slant_dist] = get_reflection_spherical_vuorinen (e, Ha, Rs, choice_method)
% GET_REFLECTION_SPHERICAL_VUORINEN Calculates reflection on spherical 
% surface based on Vuorinen (2023) equations (internal document - filename cal20230511.tex)
%
% On the Ptolemy-Alhazen problem: source at infinite distance. 2023, May
% 11. Internal document. 
% 
% INPUT:
% - Ha: antenna/receiver height (in meters)
% - e: elevation angle (in radians)
% - Rs: Earth radius (in meters)
% 
% NOTE: transmitter/satellite height is assumed infinite.
% 
% OUTPUT:
% - x_spec, y_spec: reflection point in local frame (vector, in meters)
% - dx_trans, dy_trans: transmitter direction in local frame (unit vector, in meters)
% - graz_ang: grazing angle of spherical reflection that satisfies Snell's Law (in degrees)
% - geo_ang_as: geocentric angle between antenna and reflection point (in degrees) 

if (nargin < 4), choice_method = 'heuristic';  end

% Antenna radius
Ra = Rs+Ha;

% Normalized antenna radius
ra = Ra./Rs;

% Phase angle
theta = deg2rad(90-e);

% Quartic polynomial coefficients
c4 = ra.*exp(-1i.*theta);
c3 = -1;
c2 = 0;
c1 = 1;
c0 = -ra.*exp(1i.*theta);

% Polynomial roots
ws = roots ([c4 c3 c2 c1 c0]);

% Candidate phase angles
phis = angle(ws);

% Satellite direction (unit vector)
% in either local or quasigeocentric frames:
dx_trans = cosd(e);
dy_trans = sind(e);
dir_trans = [dx_trans, dy_trans];

% Antenna's local position
pos_ant_local = [0 Ha]; 

%% Root choice process
if strcmp(choice_method,'heuristic')

    phis2 = phis;
    phis2(phis>deg2rad(90)) = NaN;
    ind = argmin(abs(phis2-theta));
    
elseif strcmp(choice_method,'rigorous')

    di0 = 3*(Ha);

    for i=1:numel (ws)

       w1 = ws(i);
       pos_spec1 = pos_spec_complex (w1, theta, Rs);
       [di1] = interf_delay_infinite_distance (pos_spec1, pos_ant_local, dir_trans);
       
       if di1<di0
           di0 = di1; 
           ind0 = i;
       end
    end
    ind = ind0;
end

w = ws(ind);
phi = phis(ind);
%%
% Geocentric angle at reflection point:
geo_ang_as = rad2deg(theta-phi);

% Reflection point in a quasigeocentric frame
[~, x_spec, y_spec] = pos_spec_complex (w, theta, Rs);
pos_spec = [x_spec y_spec];
pos_spec_geo = [x_spec y_spec+Rs];

% Interferometric delay
[di, slant_dist] = interf_delay_infinite_distance (pos_spec, pos_ant_local, dir_trans);

% Antenna geocentric position
pos_ant_geo = [0 Ra];

% Grazing angle
graz_ang = get_grazing_angle_infinite (pos_ant_geo, pos_spec_geo, dir_trans);

end 

function [pos_spec, x_spec, y_spec] = pos_spec_complex (w, theta, Rs)

    pos_spec_complex = w./exp(-1i*(pi./2-theta)).*Rs-complex(0,Rs);
    x_spec = real(pos_spec_complex);
    y_spec = imag(pos_spec_complex);
    pos_spec = [x_spec, y_spec];

end

function [di,slant_dist] = interf_delay_infinite_distance (pos_spec, pos_ant_local, dir_trans)

    slant_dist = norm_all(pos_ant_local - pos_spec);
    [~, delay_direct_aux] = proj_pt_line (pos_spec, pos_ant_local, dir_trans);
    di = slant_dist - delay_direct_aux;

end