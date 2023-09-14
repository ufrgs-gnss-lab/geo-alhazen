% DEMO_GEO_CARTESIAN: demonstration of coordinate conversion of reflection point from 2D local 
% cartesian to 3D in cartesian (local and geocentric) and geodetic frames
setup_spherical_reflection()

%% Setup
e = 10:2:30;
H = 10;
[~,n]=size(H);
[e_grid, H_grid]=ndgrid(e,H);
pos_ant_geod = zeros (n,3);
azim = 45;
ell = 'wgs84';

%% Computation
[pos_spec] = get_osculating_spherical_reflection (pos_ant_geod, azim, e_grid, H_grid, [], [], [], ell);

pos_spec_local = [pos_spec.localcart.x(:,1) pos_spec.localcart.y(:,1) pos_spec.localcart.z(:,1)];
pos_spec_geoc = [pos_spec.geocart.x(:,1) pos_spec.geocart.y(:,1) pos_spec.geocart.z(:,1)];
pos_spec_geod = [pos_spec.geod.lat(:,1) pos_spec.geod.long(:,1) pos_spec.geod.h(:,1)];
