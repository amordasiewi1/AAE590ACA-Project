%% CONSTANTS
MU_dim = 132712440017.99; % [km^3/s^2] grav param of the sun
AU_dim = 149597870.7;     % [km]

%% NON-DIMENSIONALIZATION
LU = AU_dim;
TU = sqrt(LU^3 / MU_dim);
MU = (TU^2/LU^3) * MU_dim;

%% TRANSFER PARAMETERS - FIXED TIME
t0 = 0;
tf_dim = 2 * (365*24*60*60); % [s]
tf = (1/TU) * tf_dim;

%% DEPARTURE PARAMETERS - EARTH
% keplerian elements
ear.sma  = 1.0 * AU_dim;
ear.ecc  = 0;
ear.incl = 0;
ear.raan = 0;
ear.argp = 0;
ear.ta   = 0;

% conversion to cartesian elements
[ear.r0_dim, ear.v0_dim] = keplerian_to_cartesian(...
    ear.sma, ear.ecc, ear.incl, ear.raan, ear.argp, ear.ta, MU_dim);

% non-dimensionalize
ear.r0 = (1 / LU) * ear.r0_dim;
ear.v0 = (TU / LU) * ear.v0_dim;

%% TARGET (ARRIVAL) PARAMETERS - ASTEROID
% keplerian elements
ast.sma  = 2.54015762 * AU_dim;       % [km]
ast.ecc  = 0.14886466;            %
ast.incl = deg2rad(15.26079519);  % [rad]
ast.raan = deg2rad(16.45332772);  % [rad]
ast.argp = deg2rad(26.37519590);  % [rad]

ast.mm = sqrt(MU_dim/ast.sma^3);
ast.ma = deg2rad(27.41742926) + ast.mm*(tf_dim); % [rad]
ast.ea = mean_to_eccentric(ast.ma,ast.ecc);
ast.ta = eccentric_to_true(ast.ea,ast.ecc);

% conversion to cartesian elements
[ast.rf_dim, ast.vf_dim] = keplerian_to_cartesian(...
    ast.sma,ast.ecc,ast.incl,ast.raan,ast.argp,ast.ta,MU_dim);

% non-dimensionalize
ast.rf = (1 / LU) * ast.rf_dim;
ast.vf = (TU / LU) * ast.vf_dim;

%% ODE PARAMS
ode_tol = 1e-12;
ode_opts = odeset(reltol=ode_tol, abstol=ode_tol);

%% INITIAL COSTATE SHOOTING
% fsolve options
fsolve_tol = 1e-8;
fsolve_opts = optimoptions('fsolve',...
    Display='iter',...
    FunctionTolerance=fsolve_tol,...
    OptimalityTolerance=fsolve_tol,...
    StepTolerance=fsolve_tol);

% spacecraft system control
B = [zeros(3,3); eye(3)];
% u_max = 0.1; % same as ps4
thrust = 0.8 / 1000; % [kN] tpyical eprop thurst: 10-300 mN
m_sc = 500; % [kg]
accel_dim = thrust / m_sc;
u_max = (TU^2 / LU) * accel_dim;

% initial state
x0 = [ear.r0; ear.v0];

% final target state
target_xf = [ast.rf; ast.vf];

% initial costate guess from minimum energy case
rho_vec= [1, 0.1, 1e-2, 1e-3];
lambda0_guess = load('min_energy_tf2_lambda0.mat').lambda0;
lambda0_all = zeros(length(lambda0_guess),length(rho_vec));

% solve
for i = 1:length(rho_vec)
    rho = rho_vec(i);
    lambda0_guess = fsolve(@(lambda0) asteroid_min_fuel_shooting(...
        lambda0, x0, target_xf, MU, B, u_max, rho, t0, tf, ode_opts), ...
        lambda0_guess, fsolve_opts);
    lambda0_all(:,i) = lambda0_guess;
end

%% INTEGRATE SOLUTION
% integration params
N = 10000;
t = linspace(t0,tf,N);
y0 = [x0; lambda0_all(:,end)];
rho = rho_vec(end);

% integrate and extract soln
traj = ode45(@(t,y) asteroid_min_fuel_ode(t,y,MU,B,u_max,rho),...
    [t0 tf], y0, ode_opts);
Y = deval(traj,t);
r = Y(1:3,:);
v = Y(4:6,:);
lambda = Y(7:12,:);

% recompute control history
u = zeros(3,N);
for j = 1:N
    p = -B'*lambda(:,j);
    uhat = p/norm(p);
    S = norm(p) - 1;
    gamma = (u_max/2)*(1 + tanh(S/rho));
    u(:,j) = gamma * uhat;
end

% compute fuel cost
fuel_cost = trapz(t, vecnorm(u,2,1));

%% PLOT
% 'propagate' earth/asteroid orbits
M = 10000;
theta = linspace(0,2*pi,M);
ear.r = (1/LU) * ear.sma * [sin(theta); cos(theta); zeros(1,M)];
ast.r = zeros(3,M);
for k = 1:M
    [ast.r(:,k), ~] = keplerian_to_cartesian(...
        ast.sma,ast.ecc,ast.incl,ast.raan,ast.argp,theta(k),MU_dim);
    ast.r(:,k) = (1/LU) * ast.r(:,k);
end

% trajectory plot
f1 = figure(1); hold on; grid on; axis equal;
plot3(ear.r(1,:),ear.r(2,:),ear.r(3,:),':',...
    color=[0.5 0.5 0.5],...
    linewidth=1,...
    HandleVisibility='off')
plot3(ast.r(1,:),ast.r(2,:),ast.r(3,:),':',...
    color=[0.5 0.5 0.5],...
    linewidth=1,...
    HandleVisibility='off')
plot3(r(1,:),r(2,:),r(3,:),'k',...
    linewidth=2,...
    displayname='Trajectory');
xlabel('$$r_1$$ [AU]')
ylabel('$$r_2$$ [AU]')
zlabel('$$r_3$$ [AU]')

% control history plot
f2 = figure(2); hold on; grid on;
plot(t,u(1,:),...
    displayname='$$u_1$$')
plot(t,u(2,:),...
    displayname='$$u_2$$')
plot(t,u(3,:),...
    displayname='$$u_3$$')
plot(t,vecnorm(u,2,1),...
    displayname='$$||u||_2$$')
xlabel('Time [ndim]')
ylabel('$$u_i$$ [ndim]')
legend

%% trajectory plot with thrust / coast coloring
f4 = figure(4); hold on; grid on; axis equal;

% reference orbits
plot3(ear.r(1,:),ear.r(2,:),ear.r(3,:),':', ...
    color=[0.5 0.5 0.5], ...
    linewidth=1, ...
    HandleVisibility='off')
plot3(ast.r(1,:),ast.r(2,:),ast.r(3,:),':', ...
    color=[0.5 0.5 0.5], ...
    linewidth=1, ...
    HandleVisibility='off')

% thrust magnitude
umag = vecnorm(u,2,1);

% thrusting threshold
tol = 1e-10;
isThrust = umag > tol;

% plot trajectory by thrust/coast segments
N = length(t);
idx = 1;

thrust_plotted = false;
coast_plotted = false;

while idx < N
    mode = isThrust(idx);
    j = idx;
    while j < N && isThrust(j) == mode
        j = j + 1;
    end
    if mode
        if ~thrust_plotted
            plot3(r(1,idx:j), r(2,idx:j), r(3,idx:j), ...
                'r', linewidth=2, HandleVisibility='off')
            thrust_plotted = true;
        else
            plot3(r(1,idx:j), r(2,idx:j), r(3,idx:j), ...
                'r', linewidth=2, HandleVisibility='off')
        end
    else
        if ~coast_plotted
            plot3(r(1,idx:j), r(2,idx:j), r(3,idx:j), ...
                'k', linewidth=2, HandleVisibility='off')
            coast_plotted = true;
        else
            plot3(r(1,idx:j), r(2,idx:j), r(3,idx:j), ...
                'k', linewidth=2, HandleVisibility='off')
        end
    end

    idx = j;
end

% thrust arrows
skip = 100;                 % increase to reduce arrows
arrow_scale = 0.5;        % visual scaling in AU-ish units
idxq = 1:skip:N;
idxq = idxq(isThrust(idxq));   % only plot arrows where thrusting

% normalize thrust direction
uhat = u(:,idxq) ./ vecnorm(u(:,idxq),2,1);

quiver3(r(1,idxq), r(2,idxq), r(3,idxq), ...
        arrow_scale*uhat(1,:), ...
        arrow_scale*uhat(2,:), ...
        arrow_scale*uhat(3,:), ...
        0, ...
        'b', ...
        linewidth=1.2, ...
        MaxHeadSize=1.5, ...
        displayname='Thrust direction')

xlabel('$$r_1$$ [AU]')
ylabel('$$r_2$$ [AU]')
zlabel('$$r_3$$ [AU]')