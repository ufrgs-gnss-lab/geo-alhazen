function [Di, sldist, g, gamma, der] = get_spherical_reflection_trig (e, Ha, Rs, optnum)

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
% Dl: slant distance between receiver and reflection point (matrix, in meters)

if (nargin < 3) || isempty(Rs);  Rs = get_earth_radius();  end
if (nargin < 4) || isempty(optnum);  optnum = struct();  end
%%
[g,e_spec] = get_spherical_finite (e, Ha, Rs, optnum);
Ra = Rs+Ha;
gamma = (Ra./Rs).^2-cosd(g).^2;
sldist = Rs.*(sqrt(gamma)-sind(g));
der = e-e_spec; % or der = -(-e+e_spec) (master thesis incorrect)
Di = sldist.*(1-cosd(der));
% Di1 = sldist.*(1-cosd(2.*e+der)); % #WRONG! (master thesis incorrect)

end

function [g,e_spec] = get_spherical_finite (e, Ha, Rs, optnum)
if isfieldempty (optnum, 'Ht'),  optnum.Ht = [];  end
if isfieldempty (optnum, 'algorithm'),  optnum.algorithm = [];  end
if isfieldempty (optnum, 'trajectory'),  optnum.trajectory = [];  end
if isfieldempty (optnum, 'frame'),  optnum.frame = [];  end
Ht = optnum.Ht;
algorithm = optnum.algorithm;
trajectory = optnum.trajectory;
frame = optnum.frame;

[~, g, ~, ~, ~, ~, ~, ~, e_spec] = ...
    get_reflection_spherical (e, Ha, Ht, Rs, algorithm, trajectory, frame); %Based on finite satellite distance

end