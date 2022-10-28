%% simpler demo to retrieve interf. delay for elevation variation

Ha_dom = 10;
ehor = get_spherical_reflection_horizon_elev (Ha_dom);
e_dom = ehor:1:90;
[e_grid, Ha_grid] = ndgrid (e_dom, Ha_dom);
frame = 'quasigeo';

[Di] = get_reflection_spherical (e_grid, Ha_grid, [], [], [], [], frame);

figure
plot (e_dom, Di)
xlabel ('Elevation angle (degrees)')
ylabel ('Interf. delay (m)')
grid on
legend ('10-m antenna height', 'Location', 'best')
xlim ([ehor 90])