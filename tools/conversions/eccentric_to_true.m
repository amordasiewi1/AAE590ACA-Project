function ta = eccentric_to_true(ea, ecc)
% ECCENTRIC_TO_TRUE converts eccentric anomaly to true anomaly (only for
% closed orbit)
%
% INPUTS
%   ea  - eccentric anomaly [rad]
%   ecc - eccentricity      [ndim]
%
% OUTPUTS
%   ta - true anomaly [rad]

if ecc >= 1
    warning('eccentric_to_true() only valid for closed orbits (ecc < 1)')
end

ta = 2*atan(sqrt((1+ecc)./(1-ecc)).*tan(ea/2));