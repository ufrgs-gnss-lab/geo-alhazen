% Demonstrate retunrning spherical reflection variables with a height grid
% and a elevation grid as input

setup_spherical_reflection () % Initial setup

e_dom = 1:5:90;
Ha_dom = [10 50 100 200 300 500 1000];
[e_grid, Ha_grid] = meshgrid(e_dom, Ha_dom);

[Di, g, arclen, sldist, xspec, yspec] = get_reflection_spherical (e_grid, Ha_grid);


%%
figure, imagesc(e_dom, Ha_dom, Di)
xlabel('Sat. Elev. (deg.)')
ylabel('Ant. Height (h)')
title('Interf. delay (m)')
grid on
set(gca(), 'YDir','normal')

%%
figure, surf(e_dom, Ha_dom, Di)
xlabel('Sat. Elev. (deg.)')
ylabel('Ant. Height (h)')
title('Interf. delay (m)')
grid on
axis vis3d

