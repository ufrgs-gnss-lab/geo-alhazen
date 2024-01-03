e = 30;
Ha = 10;

R0 = get_earth_radius ();
Ra = R0 + Ha;
theta = deg2rad(90-e);
ra = Ra./R0;

[~, ~, ~, ~, x_spec, y_spec] = get_spherical_reflection (e, Ha, [], [], [], [], 'local');
pos_spec_ref = [x_spec y_spec];

%%
% Coefficients
c4 = ra.*exp(-1i.*theta);
c3 = -1;
c2 = 0;
c1 = 1;
c0 = -ra.*exp(1i.*theta);

ws = roots ([c4 c3 c2 c1 c0]);

test = ws.*conj(ws);
cond = abs(test)-1;
[~,ind] = min(abs(cond));

w1 = ws(ind);
pos_spec_qgeor = [real(w1) imag(w1)].*R0;

phi = atand (pos_spec_qgeor(2)./pos_spec_qgeor(1));

x = R0.*cosd(e).*cosd(phi)-R0.*sind(e).*sind(phi);
y = R0.*sind(e).*cosd(phi)+R0.*cosd(e).*sind(phi);


%%

% pos_ant_qgeor = Ra.*[tand(theta) cosd(theta)]; 
% 
% pos_ant_spec = pos_ant_qgeor-pos_ref_qgeor;
% norm_pos_ant_spec = norm_all(pos_ant_spec);
% dir_ant_spec = pos_ant_spec./norm_pos_ant_spec;
% 
% pos_ref_qgeo = norm_pos_ant_spec.*dir_ant_spec;


