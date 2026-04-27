function psi = asteroid_min_energy_shooting(lambda0, x0, target_xf, MU, B, ...
    t0, tf, ode_opts)
% initial state augmented with initial costate guess
y0 = [x0; lambda0];

% integrate
[~,Y] = ode45(@(t,y) asteroid_min_energy_ode(t,y,MU,B),...
    [t0 tf], y0, ode_opts);
xf = Y(end,1:6)';

% boundary condition constraint
psi = xf - target_xf;