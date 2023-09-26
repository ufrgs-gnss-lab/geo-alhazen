setup_spherical_reflection()

e = 5:1:30;
Ha = [10 50];
Ht = get_satellite_height();
Rs = get_earth_radius();

for i=1:numel(e)
    for j=1:numel(Ha)
        e_aux = e(i);
        Ha_aux = Ha(j);
        [graz_ang_v(i,j), geo_ang_v(i,j), x_spec_v(i,j), y_spec_v(i,j), x_trans_v(i,j), y_trans_v(i,j)] = get_reflection_spherical_vuorinen (e_aux, Ha_aux, Ht, Rs);
        [graz_ang_f(i,j), geo_ang_f(i,j), x_spec_f(i,j), y_spec_f(i,j), x_trans_f(i,j), y_trans_f(i,j)] = get_reflection_spherical_fujimura (e_aux, Ha_aux, Ht, Rs);
    end
end

dg = graz_ang_v-graz_ang_f;
dgeo_ang = geo_ang_v-geo_ang_f;
dx = x_spec_v-x_spec_f;
dy = y_spec_v-y_spec_f;
dxt = x_trans_v-x_trans_f;
dyt = y_trans_v-y_trans_f;