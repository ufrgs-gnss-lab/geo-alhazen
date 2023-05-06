% Comparison to Fujimura (3D Plots)
setup_spherical_reflection()
%% Input values

% Ha_dom = [10 50 100 200 300 500 1000];
Ha_dom = 10:10:1000;
e_dom = [90 45 30 10 5 0];
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
[Di_c90, Di_c45, Di_c30, Di_c10, Di_c5, Di_c0, comp_algs] = get_comp_elev (fuji.Di, mn.Di, helm.Di, mil.Di, fer.Di, e_dom);
[g_c90, g_c45, g_c30, g_c10, g_c5, g_c0] = get_comp_elev (fuji.g, mn.g, helm.g, mil.g, fer.g, e_dom);
[arclen_c90, arclen_c45, arclen_c30, arclen_c10, arclen_c5, arclen_c0] = get_comp_elev (fuji.arclen, mn.arclen, helm.arclen, mil.arclen, fer.arclen, e_dom);
[sldist_c90, sldist_c45, sldist_c30, sldist_c10, sldist_c5, sldist_c0] = get_comp_elev (fuji.sldist, mn.sldist, helm.sldist, mil.sldist, fer.sldist, e_dom);
[X_spec_c90, X_spec_c45, X_spec_c30, X_spec_c10, X_spec_c5, X_spec_c0] = get_comp_elev (fuji.X_spec, mn.X_spec, helm.X_spec, mil.X_spec, fer.X_spec, e_dom);
[Y_spec_c90, Y_spec_c45, Y_spec_c30, Y_spec_c10, Y_spec_c5, Y_spec_c0] = get_comp_elev (fuji.Y_spec, mn.Y_spec, helm.Y_spec, mil.Y_spec, fer.Y_spec, e_dom);

%% Plot comparison Di
t = [1:1:4];

figure
subplot (2,2,1)
    imagesc(t, Ha_dom, Di_c90')
    title ('e=90°')
    colorbar
subplot (2,2,2)
    imagesc(t, Ha_dom, Di_c45')
    title ('e=45°')
    colorbar
subplot (2,2,3)
    imagesc(t, Ha_dom, Di_c30')
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    imagesc(t, Ha_dom, Di_c10')
    title ('e=10°')      
    colorbar
sgtitle ('Interf. Delay')

figure
subplot (2,2,1)
    surf(Ha_dom, t', Di_c90)
    title ('e=90°')
    colorbar
subplot (2,2,2)
    surf(Ha_dom, t', Di_c45)
    title ('e=45°')
    colorbar
subplot (2,2,3)
    surf(Ha_dom, t', Di_c30)
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    surf(Ha_dom, t', Di_c10)
    title ('e=10°')      
    colorbar
sgtitle ('Interf. Delay')

%% Plot comparison g
t = [1:1:4];

figure
subplot (2,2,1)
    imagesc(t, Ha_dom, g_c90')
    title ('e=90°')
    colorbar
subplot (2,2,2)
    imagesc(t, Ha_dom, g_c45')
    title ('e=45°')
    colorbar
subplot (2,2,3)
    imagesc(t, Ha_dom, g_c30')
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    imagesc(t, Ha_dom, g_c10')
    title ('e=10°')      
    colorbar
sgtitle ('Grazing angle')

figure
subplot (2,2,1)
    surf(Ha_dom, t', g_c90)
    title ('e=90°')
    colorbar
subplot (2,2,2)
    surf(Ha_dom, t', g_c45)
    title ('e=45°')
    colorbar
subplot (2,2,3)
    surf(Ha_dom, t', g_c30)
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    surf(Ha_dom, t', g_c10)
    title ('e=10°')      
    colorbar
sgtitle ('Grazing angle')

%% Plot comparison arc length
t = [1:1:4];

figure
subplot (2,2,1)
    imagesc(t, Ha_dom, arclen_c90')
    title ('e=90°')
    colorbar
subplot (2,2,2)
    imagesc(t,  Ha_dom, arclen_c45')
    title ('e=45°')
    colorbar
subplot (2,2,3)
    imagesc(t,  Ha_dom, arclen_c30')
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    imagesc(t, Ha_dom, arclen_c10')
    title ('e=10°')      
    colorbar
sgtitle ('Arc length')

figure
subplot (2,2,1)
    surf(Ha_dom, t', arclen_c90)
    title ('e=90°')
    colorbar
subplot (2,2,2)
    surf(Ha_dom, t', arclen_c45)
    title ('e=45°')
    colorbar
subplot (2,2,3)
    surf(Ha_dom, t', arclen_c30)
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    surf(Ha_dom, t', arclen_c10)
    title ('e=10°')      
    colorbar
sgtitle ('Arc length')

%% Plot comparison slant distance
t = [1:1:4];

figure
subplot (2,2,1)
    imagesc(t, Ha_dom, sldist_c90')
    title ('e=90°')
    colorbar
subplot (2,2,2)
    imagesc(t, Ha_dom, sldist_c45')
    title ('e=45°')
    colorbar
subplot (2,2,3)
    imagesc(t, Ha_dom, sldist_c30')
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    imagesc(t, Ha_dom, sldist_c10')
    title ('e=10°')      
    colorbar
sgtitle ('Slant distance')

figure
subplot (2,2,1)
    surf(Ha_dom, t', sldist_c90)
    title ('e=90°')
    colorbar
subplot (2,2,2)
    surf(Ha_dom, t', sldist_c45)
    title ('e=45°')
    colorbar
subplot (2,2,3)
    surf(Ha_dom, t', sldist_c30)
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    surf(Ha_dom, t', sldist_c10)
    title ('e=10°')      
    colorbar
sgtitle ('Slant distance')

%% Plot comparison X
t = [1:1:4];

figure
subplot (2,2,1)
    imagesc(t', Ha_dom, X_spec_c90')
    title ('e=90°')
    colorbar
subplot (2,2,2)
    imagesc(t', Ha_dom, X_spec_c45')
    title ('e=45°')
    colorbar
subplot (2,2,3)
    imagesc(t', Ha_dom, X_spec_c30')
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    imagesc(t', Ha_dom, X_spec_c10')
    title ('e=10°')      
    colorbar
sgtitle ('X Coordinate')

figure
subplot (2,2,1)
    surf(Ha_dom, t', X_spec_c90)
    title ('e=90°')
    colorbar
subplot (2,2,2)
    surf(Ha_dom, t', X_spec_c45)
    title ('e=45°')
    colorbar
subplot (2,2,3)
    surf(Ha_dom, t', X_spec_c30)
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    surf(Ha_dom, t', X_spec_c10)
    title ('e=10°')      
    colorbar
sgtitle ('X Coordinate')

%% Plot comparison Y
t = [1:1:4];

figure
subplot (2,2,1)
    imagesc(t, Ha_dom, Y_spec_c90')
    title ('e=90°')
    colorbar
subplot (2,2,2)
    imagesc(t, Ha_dom, Y_spec_c45')
    title ('e=45°')
    colorbar
subplot (2,2,3)
    imagesc(t, Ha_dom, Y_spec_c30')
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    imagesc(t, Ha_dom, Y_spec_c10')
    title ('e=10°')      
    colorbar
sgtitle ('Y Coordinate')

figure
subplot (2,2,1)
    surf(Ha_dom, t', Y_spec_c90)
    title ('e=90°')
    colorbar
subplot (2,2,2)
    surf(Ha_dom, t', Y_spec_c45)
    title ('e=45°')
    colorbar
subplot (2,2,3)
    surf(Ha_dom, t', Y_spec_c30)
    title ('e=30°')     
    colorbar
subplot (2,2,4)
    surf(Ha_dom, t', Y_spec_c10)
    title ('e=10°')      
    colorbar
sgtitle ('Y Coordinate')

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
function [comp90, comp45, comp30, comp10, comp5, comp0, comp_algs] = get_comp_elev (fuji, mn, helm, mil, fer, e_dom)
       
    c1 = fuji-mn;
    c2 = fuji-helm;
    c3 = fuji-mil;
    c4 = fuji-fer;
    
    for i=1:numel(e_dom)
        switch i
            case 1
                comp90 = [c1(i,:); c2(i,:); c3(i,:); c4(i,:)];
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
    
    comp_algs = {'1=Fuj-MN', '2=Fuj-H','3=Fuj-Mil','4=Fuj-Fer'};

end