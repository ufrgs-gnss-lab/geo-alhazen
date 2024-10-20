function [di,slant_dist] = get_delay_infinite_trans (pos_spec, pos_ant, dir_trans)

% specular and antenna positions must be on same frame

    slant_dist = norm_all(pos_ant - pos_spec);
    [~, delay_direct_aux] = proj_pt_line (pos_spec, pos_ant, dir_trans);
    di = slant_dist - delay_direct_aux;

end