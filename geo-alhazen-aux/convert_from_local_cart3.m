function [pt_geoc_cart, pt_geod, pt_local_cart] = convert_from_local_cart3 (azim, pt_local_cart2, base_geod, ell)
% CONVERT_FROM_LOCAL_CART convert the local quasigeocentric cartesian,
% originally on 2D frame alligned to the satellite azimuth
 
    n = size(pt_local_cart2,1);
    
    % Define 2D position [X Y] to 3D [X 0 Y]
    pt_local_cart3 = [pt_local_cart2(:,1) zeros(n,1) pt_local_cart2(:,2)];
    
    % Rotate the local position over Z axis with azimuth
    pt_local_cart = myrotate([0,0,90-azim], pt_local_cart3);
    
    % Convert to geocentric cartesian and geodetic coordinates
    [pt_geoc_cart, pt_geod] = convert_from_local_cart2 (pt_local_cart, base_geod, ell);
end