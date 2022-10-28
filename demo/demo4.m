%% Demonstrate the trigonometric approximation formulas

setup_spherical_reflection () % Initial setup

e_dom = 1:1:90;
Ha_dom = [10 50 100 200 300 500 1000];
[e_grid, Ha_grid] = meshgrid(e_dom, Ha_dom);
Rs = get_earth_radius(); 

[Di1, g1, ~, sldist1] = get_reflection_spherical (e_grid, Ha_grid);
[Di2, sldist2, g2] = get_spherical_reflection_trig (e_grid, Ha_grid, Rs);

%%
figure
plot (e_dom, Di1-Di2)
xlabel ('elevation angle (degree)')
ylabel ('\deltaD_i')
grid on

figure
plot (e_dom, sldist1-sldist2)
xlabel ('elevation angle (degree)')
ylabel ('\deltaD_l')
grid on