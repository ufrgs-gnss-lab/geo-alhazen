function [pos_spec, delay, graz_ang, arc_len, slant_dist, x_spec, y_spec, x_trans, y_trans, elev_spec] = ...
get_osculating_spherical_reflection (base_geod, azim, e, Ha, Ht, algorithm, trajectory, ell_name, hs)
% GET_OSCULATING_SPHERICAL_REFLECTION Calculates specular reflection on
% an osculating sphere and retrieve the specular reflection point on
% local quasigeocentric cartesian, geocentric cartesian and geodetic
% coordinates

% INPUT:
% -base_geod: 3D geodetic coordinates of the base point or antenna's
% -foot [lat long height] (array, degrees and meter)
% -az: azimuth of the satellite (matrix, in degrees)
% -e, Ha, Ht, algorithm and trajectory, see get_spherical_reflection.m
%   documentation
% -ell_name: ellipsoid name (string). See get_ellipsoid.m documentation
% -hs: height between osculating sphere to surface

% OUTPUT
% -pos_spec: 3D position of the specular point on different frames (struct)
% -delay, graz_ang, arc_len, slant_dist, x_spec, y_spec, x_trans, y_trans,
%   elev_spec: see get_spherical_reflection.m documentation
    
   if (nargin<5) || isempty(Ht); Ht=[]; end
   if (nargin<6) || isempty(algorithm); algorithm=[];  end
   if (nargin<7) || isempty(trajectory);trajectory=[]; end
   if (nargin<8) || isempty(ell_name); ell_name = 'wgs84'; end
   if (nargin<9) || isempty(hs); hs=0; end
   
%%   
   ell = get_ellipsoid(ell_name); 
   Rs = get_radius_surface (base_geod(1), hs, ell);
   frame = 'quasigeo';
   
   %% Compute reflection on quasi-geocentric frame
   [delay, graz_ang, arc_len, slant_dist, x_spec, y_spec, x_trans, y_trans, elev_spec] = ...
            get_reflection_spherical(e, Ha, Ht, Rs, algorithm, trajectory, frame);
   pos_spec = struct();
   [~,n] = size(Ha);
   
   for i=1:n
       pt_local_cart2 = [x_spec(:,i) y_spec(:,i)];
       [pt_geoc_cart, pt_geod, pt_local_cart] = convert_from_local_cart3 (azim, pt_local_cart2, base_geod(i,:), ell);
       
       % Position on geocentric cartesian frame
       pos_spec.geocart.x(:,i) = pt_geoc_cart(:,1);
       pos_spec.geocart.y(:,i) = pt_geoc_cart(:,2);
       pos_spec.geocart.z(:,i) = pt_geoc_cart(:,3);
       
       % Position on local quasigeocentric cartesian frame
       pos_spec.qgeocart.x(:,i) = pt_local_cart(:,1);
       pos_spec.qgeocart.y(:,i) = pt_local_cart(:,2);
       pos_spec.qgeocart.z(:,i) = pt_local_cart(:,3);
       
       % Position on geodectic frame
       pos_spec.geod.lat(:,i) = pt_geod(:,1);
       pos_spec.geod.long(:,i) = pt_geod(:,2);
       pos_spec.geod.h(:,i) = pt_geod(:,3);
   end
   
   %% Transform the origin  #test    
%    [pos_cart_base0, Rg] = get_sphere_osculating (base_geod, ell); % cartesian base point
%    lat0 = deg2rad(pos_cart_base0(:,1));
%    long0 = deg2rad(pos_cart_base0(:,2));
%    h0 = pos_cart_base0(:,3);
%    
%    n0 = compute_prime_vertical_radius (lat0, ell);
%    pos_cart_center0 = pos_cart_base0 - Rs.*n0; % cartesian center (through eq. 36 of dissertation)
%    assignin('base','pos_center',pos_cart_center0);
%    assignin('base','pos_ant',pos_cart_base0);
   
   
end