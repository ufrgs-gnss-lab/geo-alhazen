%% demo1_rtr_sph: first testing of raytracing on a sphere
clear
% addpath(genpath('C:/MatLabFiles/atm-interf-sph-dev'));
addpath(genpath('C:/MatLabFiles/m'));
% setup_atm_interf ()
setup_spherical_reflection ()

%% Setup
e = [1:1:30]';
H = 50;
[~,n]=size(H);
[e_grid, H_grid]=ndgrid(e,H);
pos_ant_geod = zeros (n,3);
azim = 45;

%% Computation
[pos_spec, ~, ~, ~, ~, x_spec, y_spec] = get_osculating_spherical_reflection (pos_ant_geod, azim, e_grid, H_grid);

pos_spec_qgeoc = [pos_spec.qgeocart.x(:,1) pos_spec.qgeocart.y(:,1) pos_spec.qgeocart.z(:,1)];
pos_spec_geoc = [pos_spec.geocart.x(:,1) pos_spec.geocart.y(:,1) pos_spec.geocart.z(:,1)];
pos_spec_geod = [pos_spec.geod.lat(:,1) pos_spec.geod.long(:,1) pos_spec.geod.h(:,1)];


% %% Planar
% rtr.opt.surface_type = 'plane';
% rtr.opt.atm_planar = 1;
% [dt_p, da_p, dg_p, Ht_p, Ha_p, Hg_p, result_p, extra_p] = raytrace_interf_series2 (e, H, rtr);
% 
