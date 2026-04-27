function ydot = asteroid_min_fuel_ode(t,y,MU,B,u_max,rho)
% extract states
x = y(1:6);
r = x(1:3);
v = x(4:6);
lambda = y(7:12);

rmag = norm(r);

% optimal control
p = -B'*lambda;
uhat = p/norm(p);
S = norm(p) - 1;
gamma = (u_max/2)*(1 + tanh(S/rho));
u = gamma * uhat;

% spacecraft state derivative
f0 = [v;
      -(MU/rmag^3) * r];
xdot = f0 + B*u;

% 2bp state jacobian
dadr_2bp = (MU / rmag^3) * (3*(r*r')/rmag^2 - eye(3));
A = [zeros(3,3) eye(3);
     dadr_2bp   zeros(3,3)];

% costate derivative
lambdadot = -A'*lambda;

% return
ydot = [xdot; lambdadot];