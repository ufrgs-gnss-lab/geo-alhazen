function [pt_geoc_cart, pt_geod, pt_local_cart] = convert_from_local_cart3 (az, pt_local_cart1, pos_ant_local, base_geod, ell)
% CONVERT_FROM_LOCAL_CART convert the local quasigeocentric cartesian,
% originally on 2D frame alligned to the satellite azimuth
%
% INPUT:
% - az: satellite azimuth
% - pt_local_cart1: specular point position vector in local frame
% - pos_ant_local_cart: antenna position vector in local frame
% - base_geod: 3D geodetic coordinates of the base point or antenna's
% - ell: ellipsoid name
%
% OUTPUT
% - pt_geoc_cart: 3D geocentric cartesian coordinates of the specular point
% - pt_geod: 3D geodetic coordinates of the specular point
% - pt_local_cart: 3D local cartesian coordinates of the specular point


    n = size(pt_local_cart1,1);
    pt_local_cart2 = [pt_local_cart1(:,1) zeros(n,1) pt_local_cart1(:,2)];
    
    % Rotate Z-axis with the azimuth
    pt_local_cart = (myrotate([0,0,-az], pt_local_cart2));
    
    % Convert from local coordinates to geodetic and geocentric cartesian
    [pt_geoc_cart, pt_geod] = convert_from_local_cart2 (pt_local_cart-pos_ant_local, base_geod, ell);
    
end