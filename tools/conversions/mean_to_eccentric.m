function [ea, niter_out] = mean_to_eccentric(ma,ecc,tol)
% MEAN_TO_ECCENTRIC solves Kepler's equation to compute eccentric enomaly from mean
% anomaly
%
% INPUTS
%   ma - mean anomaly [rad]
%   e  - eccentricity [ndim]
%
% OUTPUTS
%   ea - eccentric anomaly [rad]
arguments
    ma
    ecc
    tol = 1e-12;
end

ea_curr = ma;
niter = 0;
niter_out = 0;

% trivial case
if ecc == 0
    ea = ma;
    return
end

% numerically sovlve, newton-raphson
while true
    niter = niter + 1;
    num = ea_curr - ecc .* sin(ea_curr) - ma;
    den = 1 - ecc .* cos(ea_curr);
    ea_next = ea_curr - num./den;
    
    if isnan(abs(ea_curr - ea_next))
        break
    end
    if all(abs(ea_curr - ea_next) < tol)
        break
    else
        ea_curr = ea_next;
    end
end

ea = ea_next;

if nargout > 1
    niter_out = niter;
end