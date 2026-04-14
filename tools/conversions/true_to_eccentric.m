function ea = true_to_eccentric(ta, ecc)
% TRUE_TO_ECCENTRIC converts true anomaly to eccentric anomaly
%
% INPUTS
%   ta  - true anomaly [rad]
%   ecc - eccentricity [ndim]
%
% OUTPUTS
%   ea - eccentric anomaly

ea = 2*atan(sqrt((1-ecc)./(1+ecc)).*tan(ta/2));