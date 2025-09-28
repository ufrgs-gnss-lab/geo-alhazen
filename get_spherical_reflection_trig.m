function [Di, sldist, g, gamma, der, geo_ang_as] = get_spherical_reflection_trig (e, Ha, Rs, optnum)

% GET_SPHERICAL_REFLECTION_TRIG computes interferometric delay and
% slant distance on a spherical surface with a hybrid trigonometric
% formulation that is essentially based on assumptions of infinite distance 
% of the transmitting satellite and uses variables from a finite satellite
% distance method.
%
% INPUT:
% e: elevation angles (matrix; in degrees) 
% Ha: antenna height (matrix; in meters)
% Rs: Earth surface radius (scalar, in meters)
% optnum: optional inputs for the numerical spherical reflection (see
%         get_spherical_reflection for details)(struct)
%       - Ht: transmitter height;
%       - algorithm: algorithm to compute spherical reflection
%       - trajectory: transmitter trajectory
%       - frame: coordinate reference frame;

% OUTPUT:
% Di: interferometric delay or vaccum interferometric distance (matrix, in meters)
% sldist: slant distance between receiver and reflection point (matrix, in
% meters)
% g: grazing angle of the reflection (matrix, in degrees);
% gamma: ratio of antenna radius w.r.t the surface radius (matrix, in degrees);
% der: difference direct elevation angle to the elevation angle of the
% reflection (matrix, in degrees);
% geo_ang_as: geocentric angle between antenna and reflection point (matrix, in degrees);

if (nargin < 3) || isempty(Rs);  Rs = get_earth_radius();  end
if (nargin < 4) || isempty(optnum);  optnum = struct();  end
%%
[g,e_spec, geo_ang_as] = get_spherical_finite (e, Ha, Rs, optnum);
Ra = Rs+Ha;
gamma = Ra.^2./Rs.^2-cosd(g).^2;
sldist = Rs.*(sqrt(gamma)-sind(g));
der = -e-e_spec;
Di = sldist.*(1-cosd(e-e_spec));

end

function [g,e_spec, geo_ang_as] = get_spherical_finite (e, Ha, Rs, optnum)
if isfieldempty (optnum, 'Ht'),  optnum.Ht = [];  end
if isfieldempty (optnum, 'algorithm'),  optnum.algorithm = 'millerinf';  end
if isfieldempty (optnum, 'trajectory'),  optnum.trajectory = [];  end
if isfieldempty (optnum, 'frame'),  optnum.frame = [];  end
Ht = optnum.Ht;
algorithm = optnum.algorithm;
trajectory = optnum.trajectory;
frame = optnum.frame;

[~, g, ~, ~, ~, ~, ~, ~, e_spec, ~, geo_ang_as] = ...
    get_reflection_spherical (e, Ha, Ht, Rs, algorithm, trajectory, frame); %Based on finite satellite distance

end