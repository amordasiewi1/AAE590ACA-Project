function ma = eccentric_to_mean(ea)
% ECCENTRIC_TO_MEAN converts given eccentric anomaly to mean anomaly via
% kepler's equation
%
% INPUTS
%   ea  - eccentric anomaly [rad]
%   ecc - eccentricity      [ndim]
%
% OUTPUTS
%   ma - mean anomaly [rad]

ma = ea - e.*sin(ea);