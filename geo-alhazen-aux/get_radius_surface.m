function Rs = get_radius_surface (lat, hs, ell)

    if (nargin < 1) || isempty(lat),  lat = 0;  end
    if (nargin < 2) || isempty(hs),  hs = 0;  end
    if (nargin < 3) || isempty(ell),  ell = get_ellipsoid();  end
    Rc = get_radius_gaussian(lat, ell);
    Rs = Rc + hs;
   
end