%% simpler demo to retrieve interf. delay for height variation

e_dom = 10;
Ha_dom = 0:5:500;
[e_grid, Ha_grid] = ndgrid (e_dom, Ha_dom);
frame = 'quasigeo';

[Di] = get_reflection_spherical (e_grid, Ha_grid, [], [], [], [], frame);

figure
plot (Ha_dom, Di)
xlabel ('Antenna height (m)')
ylabel ('Interf. delay (m)')
grid on
legend ('10° elevation angle', 'Location', 'best')
