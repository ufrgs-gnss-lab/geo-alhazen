% Demonstrate the relationship between height and spherical horizon 
% elevation angle

setup_spherical_reflection () % Initial setup

Ha_max = 1000;
Ha_step = 1;
Has = 0:Ha_step:Ha_max;
ehors = get_spherical_reflection_horizon_elev (Has);

figure
plot(Has, ehors, '--r','LineWidth',3)
ylim([min(ehors) 0])
ylabel ('Elevation angle to spherical horizon (degrees)')
xlabel ('Antenna height (m)')
grid on
set(gca,'FontSize',18)
