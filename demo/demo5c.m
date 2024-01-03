%% Response to 2 last comments of 2nd reviewer

setup_spherical_reflection() % Initial setup

%% Input values
Has = [10 50 100 200 300 500 1000];
ehors = get_spherical_reflection_horizon_elev (Has); % Minimum elevation angle
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
xspec = tmp;
yspec = tmp;

%% Computation of parameters for each algorithm
for i=1:m
    algorithm = algs{i};
    [Di(:,i), g(:,i), arclen(:,i), sldist(:,i), X_spec(:,i), Y_spec(:,i), X_trans(:,i), Y_trans(:,i)]...
        = get_reflection_spherical (ehors(:), Has(:), [], [], algorithm, [], frame);
    
    pos_ant = [zeros(numel(Has),1) 6370e3+Has(:)]; 
    pos_spec = [X_spec(:,i) Y_spec(:,i)];
    pos_trans = [X_trans(:,i) Y_trans(:,i)];
    
    pos_trans_spec = pos_trans - pos_spec;
    pos_ant_spec   = pos_ant   - pos_spec;
    pos_ant_trans  = pos_ant   - pos_trans;
    delay_direct(:,i) = norm_all(pos_ant_trans);
    delay_reflect_inc(:,i) = norm_all(pos_trans_spec);
    delay_reflect_out(:,i) = norm_all(pos_ant_spec);
    delay_reflect(:,i) = delay_reflect_inc(:,i) + delay_reflect_out(:,i);
    delay(:,i) = delay_reflect(:,i) - delay_direct(:,i);
    
end

%% Reference
[Diref, gref, arclenref, sldistref, X_specref, Y_specref] ... 
                   = get_spherical_horizon_params (Has(:), [], frame);

pos_spec_ref = [X_specref Y_specref];

pos_trans_spec_ref = pos_trans - pos_spec_ref;
pos_ant_spec_ref   = pos_ant   - pos_spec_ref;
pos_ant_trans_ref  = pos_ant   - pos_trans;
delay_direct_ref = norm_all(pos_ant_trans_ref);
delay_reflect_inc_ref = norm_all(pos_trans_spec_ref);
delay_reflect_out_ref = norm_all(pos_ant_spec_ref);
delay_reflect_ref = delay_reflect_inc_ref + delay_reflect_out_ref;
delay = delay_reflect_ref - delay_direct_ref;
    
%% Compute differences
dif_Di = Di - Diref;
dif_g = g - gref;
dif_X = X_spec - X_specref;
dif_Y = Y_spec - Y_specref;
dif_sd = sldist - sldistref;
dif_al = arclen - arclenref;

dif_direct = delay_direct - delay_direct_ref;
dif_reflect = delay_reflect - delay_reflect_ref;

per_dif_Di = dif_Di./Di.*100;
perc_dif_g = dif_g./g.*100;
perc_dif_X = dif_X./X_spec.*100;
perc_dif_Y = dif_Y./Y_spec.*100;
perc_dif_sd = dif_sd./sldist.*100;
perc_dif_al = dif_al./arclen.*100;

dif_direct = delay_direct - delay_direct_ref;
dif_reflect = delay_reflect - delay_reflect_ref;

