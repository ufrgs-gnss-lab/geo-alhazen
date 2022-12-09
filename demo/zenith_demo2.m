% Graphical results for the paper
setup_spherical_reflection() % Initial setup
clearvars

%% Input values
Has = 10:10:1000;
ezen = repmat(90,size(Has)); % Minimum elevation angle
algs = {'fujimura','martinneira','helm','millerandvegh','fermat'};
frame = 'quasigeo';

%% Pre-allocate data
n = numel(Has);
m = numel(algs);
tmp = NaN(n,m);
Di = tmp;
g = tmp;
arclen = tmp;
sldist = tmp;
X_spec = tmp;
Y_spec = tmp;

%% Computation of parameters for each algorithm
for i=1:m
    algorithm = algs{i};
    [Di(:,i), g(:,i), arclen(:,i), sldist(:,i), X_spec(:,i), Y_spec(:,i)]...
            = get_reflection_spherical (ezen(:), Has(:), [], [], algorithm, [], frame);
end

%% Expected values on spherical horizon
[Diref, gref, arclenref, sldistref, X_specref, Y_specref] ... 
                   = get_spherical_reflection_zenith (Has(:), [], frame);

%% Differences from expectation
dif_Di = Di - Diref;
dif_g = g - gref;
dif_X = X_spec - X_specref;
dif_Y = Y_spec - Y_specref;
dif_sd = sldist - sldistref;
dif_al = arclen - arclenref;

%% Plot setup
set(groot,'defaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)
leg = {'Fujimura','Martín-Neira','Helm','Miller','Fermat'};

%% Subplot
figure
for i=1:6  
    switch i
        case 1
            subplot(3,2,1)
            plot(Has', (dif_Di))
            ylabel('Interf. delay (m)')
        case 2
            subplot(3,2,2)
            plot(Has', (dif_g))
            ylabel(sprintf('Grazing angle\n(degrees)'))
        case 3
            subplot(3,2,3)
            plot(Has', (dif_al))
            ylabel('Arc length (m)')
        case 4
            subplot(3,2,4)
            plot(Has', (dif_sd))
            ylabel('Slant distance (m)')
        case 5
            subplot(3,2,5)
            plot(Has', (dif_X))
            ylabel('X coord. (m)')
            xlabel('Antenna Heigth (m)')
        case 6
            subplot(3,2,6)
            plot(Has', (dif_Y))
            ylabel('Y coord. (m)')   
            xlabel('Antenna Heigth (m)')
    end    
    grid on
end
legend (leg, 'Location', 'none', 'NumColumns',5)