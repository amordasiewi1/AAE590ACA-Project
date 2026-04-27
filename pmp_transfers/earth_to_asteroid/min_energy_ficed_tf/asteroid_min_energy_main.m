
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

% initial state
x0 = [ear.r0; ear.v0];

% final target state
target_xf = [ast.rf; ast.vf];

% guess
lambda0_guess = zeros(6,1);

% solve
lambda0 = fsolve(@(lambda0) asteroid_min_energy_shooting(...
    lambda0, x0, target_xf, MU, B, t0, tf, ode_opts), ...
    lambda0_guess, fsolve_opts);

%% INTEGRATE SOLUTION
% integration params
N = 10000;
t = linspace(t0,tf,N);
y0 = [x0; lambda0];

% integrate and extract soln
traj = ode45(@(t,y) asteroid_min_energy_ode(t,y,MU,B),...
    [t0 tf], y0, ode_opts);
Y = deval(traj,t);
r = Y(1:3,:);
v = Y(4:6,:);
lambda = Y(7:12,:);

% recompute control history
u = zeros(3,N);
for j = 1:N
    u(:,j) = -0.5*B'*lambda(:,j);
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

