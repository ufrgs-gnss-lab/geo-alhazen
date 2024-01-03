function Re = get_radius_global (ell, type)
% definitions: https://en.wikipedia.org/wiki/Earth_radius#Global_radii  
   
   a = ell.a; % semi-major axis
   b = ell.b; % semi-minor aixs
   e = ell.e; % eccentricity
   switch lower(type)
   case 'mean'
      Re = (2*a+b)./3;
   case 'authalic'
      Re = sqrt((a^2/2)+(b^2/2).*(atanh(e)/e));
   case {'meridional mean', 'meridional'}
      Re = ((a^(3/2)+b^(3/2))/2)^(2/3);
   end
end
