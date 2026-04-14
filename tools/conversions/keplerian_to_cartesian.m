function [r_eci, v_eci] = keplerian_to_cartesian(sma, ecc, incl, raan, argp, ta, MU)
% KEPLERIAN_TO_CARTESIAN converts keplerian orbital elements to carteisan
% ECI state vector
%
% INPUTS
%   sma  - semi-major axis         [km]
%   ecc  - eccentricity            [ndim]
%   incl - inclination             [rad]
%   raan - right asc. of asc. node [rad]
%   argp - argument of peri        [rad]
%   ta   - true anomaly            [rad]
%   MU   - GM of central body      [km^3/s^2]
%
% OUTPUTS
%   r - 3x1 ECI position [km]
%   v - 3x1 ECI velocity [km/s]
rotation_matrices;

% cartesian state in orbit (lvlh) frame
p = sma*(1-ecc^2);
hmag = sqrt(p*MU);
rmag = p/(1 + ecc*cos(ta));

r_lvlh = rmag * [cos(ta); sin(ta); 0];
v_lvlh = (MU / hmag) * [-sin(ta); ecc + cos(ta); 0];

% rotate to eci
rot = R3(argp) * R1(incl) * R3(raan); % not lvlh_to_eci - take transpose
r_eci = rot' * r_lvlh;
v_eci = rot' * v_lvlh;