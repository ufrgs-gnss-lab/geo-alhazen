setup_spherical_reflection()

Ha_dom = 100;
e_dom = 0:1:89.9;
Rs = get_earth_radius();
[e_grid, Ha_grid] = ndgrid (e_dom, Ha_dom);

[~, g, ~, ~, X_spec, Y_spec, X_trans, Y_trans] = get_reflection_spherical(e_grid, Ha_grid, [], [], [], [], 'quasigeo');
        
pos_spec(1,:) = X_spec;
pos_spec(2,:) = Y_spec; 

pos_trans(1,:) = X_trans;
pos_trans(2,:) = Y_trans; 

pos_ant (1,:) = zeros(size(e_dom));
pos_ant (2,:) = Rs+Ha_dom;

geo_ang_as_vec = get_geocentric_angle_vec (pos_ant', pos_spec');
geo_ang_at_vec = get_geocentric_angle_vec (pos_ant', pos_trans');

[geo_ang_as, geo_ang_at] = get_geocentric_angle (Ha_dom,[],e_dom,g,[]);
dif_geo_ang_as = geo_ang_as-geo_ang_as_vec;
dif_geo_ang_at = geo_ang_at'-geo_ang_at_vec;

%%
figure
plot (e_dom, dif_geo_ang_as)
xlim ([0 90])
xlabel ('Elevation angle (degrees)')
ylabel ('Differences (degrees)')
title ('Geocentric angle antenna-spec. point')

figure
plot (e_dom, dif_geo_ang_at)
xlim ([0 90])
xlabel ('Elevation angle (degrees)')
ylabel ('Differences (degrees)')
title ('Geocentric angle antenna-transmitter')