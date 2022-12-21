%% Comparison to Fujimura: Results at 45° for main manuscript result
setup_spherical_reflection()

%% Input values
Ha_dom = 10:10:1000;
e_dom = [60 45 30 10 5 0];
[e_grid, Ha_grid] = ndgrid (e_dom, Ha_dom);
algs = {'fujimura','martinneira','helm','millerandvegh','fermat'};
frame = 'quasigeo';

%% Computation of parameters for each algorithm
for i=1:numel(algs)
    algorithm = algs{i};
    [Di, g, arclen, sldist, X_spec, Y_spec]...
            = get_reflection_spherical (e_grid, Ha_grid, [], [], algorithm, [], frame);
        
        switch i
            case 1
                fuji = get_structure (Di, g, arclen, sldist, X_spec, Y_spec);   % Fujimura
            case 2
                mn = get_structure (Di, g, arclen, sldist, X_spec, Y_spec);     % Martin-Neira
            case 3
                helm = get_structure (Di, g, arclen, sldist, X_spec, Y_spec);   % Helm
            case 4
                mil = get_structure (Di, g, arclen, sldist, X_spec, Y_spec);    % Miller & Vegh
            case 5
                fer = get_structure (Di, g, arclen, sldist, X_spec, Y_spec);    % Fermat
        end
end

%% Comparison
[Di_c60, Di_c45, Di_c30, Di_c10, Di_c5, Di_c0, comp_algs] = get_comp_elev (fuji.Di, mn.Di, helm.Di, mil.Di, fer.Di, e_dom);
[g_c60, g_c45, g_c30, g_c10, g_c5, g_c0] = get_comp_elev (fuji.g, mn.g, helm.g, mil.g, fer.g, e_dom);
[arclen_c60, arclen_c45, arclen_c30, arclen_c10, arclen_c5, arclen_c0] = get_comp_elev (fuji.arclen, mn.arclen, helm.arclen, mil.arclen, fer.arclen, e_dom);
[sldist_c60, sldist_c45, sldist_c30, sldist_c10, sldist_c5, sldist_c0] = get_comp_elev (fuji.sldist, mn.sldist, helm.sldist, mil.sldist, fer.sldist, e_dom);
[X_spec_c60, X_spec_c45, X_spec_c30, X_spec_c10, X_spec_c5, X_spec_c0] = get_comp_elev (fuji.X_spec, mn.X_spec, helm.X_spec, mil.X_spec, fer.X_spec, e_dom);
[Y_spec_c60, Y_spec_c45, Y_spec_c30, Y_spec_c10, Y_spec_c5, Y_spec_c0] = get_comp_elev (fuji.Y_spec, mn.Y_spec, helm.Y_spec, mil.Y_spec, fer.Y_spec, e_dom);

%% RMSE table
rmse_45 (1,:) = sqrt (sum(g_c45',1).^2/numel(Ha_dom));
rmse_45 (2,:) = sqrt (sum(Di_c45',1).^2/numel(Ha_dom));
rmse_45 (3,:) = sqrt (sum(X_spec_c45',1).^2./numel(Ha_dom));
rmse_45 (4,:) = sqrt (sum(Y_spec_c45',1).^2 /numel(Ha_dom));
rmse_45 (5,:) = sqrt (sum(sldist_c45',1).^2/numel(Ha_dom));
rmse_45 (6,:) = sqrt (sum(arclen_c45',1).^2/numel(Ha_dom));

%% Subpplot
set(groot,'defaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)
leg = {'Martín-Neira','Helm','Miller','Fermat'};

figure
for i=1:6  
    switch i
        case 1
            subplot(3,2,1)
            plot(Ha_dom', (Di_c45))
            ylabel('Interf. delay (m)')
        case 2
            subplot(3,2,2)
            plot(Ha_dom', (g_c45))
            ylabel(sprintf('Grazing angle\n(degrees)'))
        case 3
            subplot(3,2,3)
            plot(Ha_dom', (arclen_c45))
            ylabel('Arc length (m)')
        case 4
            subplot(3,2,4)
            plot(Ha_dom', (sldist_c45))
            ylabel('Slant distance (m)')
        case 5
            subplot(3,2,5)
            plot(Ha_dom', (X_spec_c45))
            ylabel('X coord. (m)')
            xlabel('Antenna Heigth (m)')
        case 6
            subplot(3,2,6)
            plot(Ha_dom', (Y_spec_c45))
            ylabel('Y coord. (m)')   
            xlabel('Antenna Heigth (m)')
    end    
    grid on
end
legend (leg, 'Location', 'none', 'NumColumns',5)

%% Get structure
function [stc_alg] = get_structure (Di, g, arclen, sldist, X_spec, Y_spec)
stc_alg = struct();
stc_alg.Di = Di;
stc_alg.g = g;
stc_alg.arclen = arclen;
stc_alg.sldist = sldist;
stc_alg.X_spec = X_spec;
stc_alg.Y_spec = Y_spec;
end

%% Function to compare for a given elevation angle
function [comp60, comp45, comp30, comp10, comp5, comp0, comp_algs] = get_comp_elev (fuji, mn, helm, mil, fer, e_dom)
       
    c1 = fuji-mn;
    c2 = fuji-helm;
    c3 = fuji-mil;
    c4 = fuji-fer;
    
    for i=1:numel(e_dom)
        switch i
            case 1
                comp60 = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
            case 2
                comp45 = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
            case 3
                comp30 = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
            case 4
                comp10 = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
            case 5
                comp5  = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
            case 6
                comp0  = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
        end
    end
    
    comp_algs = {'Fujimura - Martin-Neira', 'Fujimura - Helm','Fujimura - Miller','Fujimura - Fermat'};

end