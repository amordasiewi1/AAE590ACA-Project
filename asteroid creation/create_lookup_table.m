%% Initialize
clc; clear; close all;

%% Declare Variables
asteroid_laydown_file = load("asteroid_laydown.mat");
all_targets = asteroid_laydown_file.asteroid_params; 
num_targets = size(all_targets, 1);
mining_values = all_targets(:, 7); 

% Physics Constants
mu_sun = 132712440017.99;       
a_earth = 149597898;            
DU = a_earth; 
TU = sqrt(DU^3 / mu_sun);
max_thrust = 0.05; 
umax_canonical = max_thrust / 1000 / (DU / TU^2); 

% Mission Timing
mission_duration_days = 3650; 
time_step_days = 30;
departure_epochs = 0:time_step_days:mission_duration_days;
num_epochs = length(departure_epochs);
max_transfer_days = 1200; % Maximum time allowed for a single leg
t_cap_canonical = (max_transfer_days * 24 * 3600) / TU;

% Pre-allocate Tensors
dv_tensor_kms = zeros(num_targets, num_targets, num_epochs);
tof_tensor_days = zeros(num_targets, num_targets, num_epochs); % New: Track actual timing

%% 1. Parallel Integration Phase
% We no longer pre-compute departure_states/arrival_states globally.
% We compute the STARTING state at the moment of departure and integrate the TARGET dynamically.
disp('Building 3D Time-Expanded Cost Tensor with Dynamic Target Integration...');
tic;

all_targets(:, 3:6) = deg2rad(all_targets(:, 3:6)); % CRITICAL: Match ACO units

parfor t_idx = 1:num_epochs
    local_dv = zeros(num_targets, num_targets);
    local_tof = zeros(num_targets, num_targets);
    t_dep_days = departure_epochs(t_idx);
    t_dep_TU = (t_dep_days * 24 * 3600) / TU
    
    for start_idx = 1:num_targets
        % Calculate starting asteroid position at THIS specific departure epoch
        x_start_base = all_targets(start_idx, 1:6)';
        n_start = sqrt(1 / x_start_base(1)^3);
        x0_curr = x_start_base;
        x0_curr(6) = x_start_base(6) + n_start * t_dep_TU;
        
        for end_idx = 1:num_targets
            if start_idx == end_idx 
                continue; 
            end
            
            % Target asteroid base parameters
            x_target_base = all_targets(end_idx, 1:6)';
            
            % DYNAMIC INTEGRATION: Target position is updated INSIDE this function
            [dv_val, ~, actual_tof_days] = rk4_dynamic_rendezvous(...
                x0_curr, x_target_base, t_dep_days, ...
                umax_canonical, t_cap_canonical, DU, TU, mu_sun);
            
            % If it converged before the cap, store values. Otherwise, NaN.
            if actual_tof_days >= (max_transfer_days - 1)
                local_dv(start_idx, end_idx) = NaN;
                local_tof(start_idx, end_idx) = NaN;
            else
                local_dv(start_idx, end_idx) = dv_val;
                local_tof(start_idx, end_idx) = actual_tof_days;
            end
        end
    end
    dv_tensor_kms(:, :, t_idx) = local_dv;
    tof_tensor_days(:, :, t_idx) = local_tof;
end
toc;

% Save both tensors. The ACO script MUST use tof_tensor_days to know when it arrived.
save('asteroid_3D_tensor.mat', 'dv_tensor_kms', 'tof_tensor_days', 'departure_epochs', 'mining_values');
disp('3D Tensor Complete with TOF tracking!');

%% Supporting Functions
function [dv_total_kms, x_final, tof_days] = rk4_dynamic_rendezvous(x0, x_target_base, t_dep_days, u_max, t_cap, DU, TU, mu_sun)
    dt = (0.5 * 24 * 3600) / TU; % 12-hour steps
    x = x0;
    dv_DU = 0;
    t_elapsed_TU = 0;
    
    n_target = sqrt(1 / x_target_base(1)^3); % Target mean motion (canonical)
    
    while t_elapsed_TU < t_cap
        % 1. Dynamically update the Target asteroid position to the CURRENT integration time
        % We add the elapsed integration time to the departure epoch
        t_current_TU = (t_dep_days * 24 * 3600 / TU) + t_elapsed_TU;
        x_target_now = x_target_base;
        x_target_now(6) = x_target_base(6) + n_target * t_current_TU;
        
        % 2. Check for Convergence (Cartesian Position & Velocity)
        % CLAMP spacecraft state to physical elliptical bounds to prevent imaginary numbers
        a_sc = max(1e-3, x(1));
        e_sc = max(1e-4, min(x(2), 0.99));
        
        % Convert Mean Anomaly (x(6)) to True Anomaly for Spacecraft
        E_sc = Keplers(mod(x(6), 2*pi), e_sc);
        ta_sc = 2 * atan(sqrt((1 + e_sc) / (1 - e_sc)) * tan(E_sc / 2));
        
        % Convert Mean Anomaly to True Anomaly for Target (Targets are natively safe)
        E_tgt = Keplers(mod(x_target_now(6), 2*pi), x_target_now(2));
        ta_tgt = 2 * atan(sqrt((1 + x_target_now(2)) / (1 - x_target_now(2))) * tan(E_tgt / 2));
        
        % Calculate Cartesian states using bounded elements and True Anomaly
        [r_sc, v_sc] = KeplerianToCartesian(a_sc, e_sc, x(3), x(4), x(5), ta_sc, 1);
        [r_tgt, v_tgt] = KeplerianToCartesian(x_target_now(1), x_target_now(2), x_target_now(3), x_target_now(4), x_target_now(5), ta_tgt, 1);
        
        dr = norm(r_sc - r_tgt);
        dv = norm(v_sc - v_tgt);
        
        % Tolerance: ~0.001 AU (150,000 km) and ~1e-3 Canonical Velocity (~30 m/s)
        if dr < 0.08 && dv < 0.05
            break;
        end
        
        % 3. RK4 Step
        [dx1, u1] = get_derivatives(x, x_target_now, u_max);
        [dx2, u2] = get_derivatives(x + 0.5*dt*dx1, x_target_now, u_max);
        [dx3, u3] = get_derivatives(x + 0.5*dt*dx2, x_target_now, u_max);
        [dx4, u4] = get_derivatives(x + dt*dx3, x_target_now, u_max);
        
        x = x + (dt/6)*(dx1 + 2*dx2 + 2*dx3 + dx4);
        dv_DU = dv_DU + (dt/6)*(u1 + 2*u2 + 2*u3 + u4);
        t_elapsed_TU = t_elapsed_TU + dt;
    end
    
    dv_total_kms = dv_DU * (DU / TU);
    x_final = x;
    tof_days = (t_elapsed_TU * TU) / (24 * 3600);
end

function [x_dot, u_mag] = get_derivatives(x, x_target, u_max)
     a = max(1e-3, x(1));                       
     e = max(1e-4, min(x(2), 0.90));           
     i = max(1e-4, min(x(3), pi - 1e-4));       
     omg = x(4); w = x(5); M = x(6);
     
     dx = zeros(6,1);
     dx(1:3) = x(1:3) - x_target(1:3);
     dx(4:6) = atan2(sin(x(4:6)-x_target(4:6)), cos(x(4:6)-x_target(4:6)));
     
     K = diag([0.5, 1.0, 1.0, 1.0, 1.0, 0.1]);
     
     E = Keplers(mod(M, 2*pi), e);
     ta = 2*atan(sqrt((1+e)/(1-e)) * tan(E / 2)); 
     
     p = a*(1-e^2); h = sqrt(p); r = p / (1 + e*cos(ta));
     b = a * sqrt(1-e^2);
     
     Bs = 1/h * [2*a^2*e*sin(ta), 2*a^2*p/r, 0;
                 p*sin(ta), (p+r)*cos(ta) + r*e, 0;
                 0, 0, r*cos(ta+w);
                 0, 0, r * sin(ta+w) / sin(i);
                 -p * cos(ta) / e, (p+r)*sin(ta) / e, -r*sin(ta+w) / tan(i);
                 (b*p*cos(ta)/(a*e)) - 2*b*r/a, -b*(p+r)*sin(ta)/(a*e), 0];
                
     u = -pinv(Bs) * (K * dx);
     if norm(u) > u_max, u = u / norm(u) * u_max; end
    
     x_dot = Bs * u; 
     x_dot(6) = x_dot(6) + sqrt(1 / (a^3)); 
     u_mag = norm(u);
end

function [r_vec, v_vec] = KeplerianToCartesian(a, e, i, omg, w, ta, mu)
    % Pure radian formulation using the Perifocal (PQW) Frame
    p = a * (1 - e^2);
    r_mag = p / (1 + e * cos(ta));
    
    % Position and Velocity in Perifocal Frame
    r_pqw = [r_mag * cos(ta); r_mag * sin(ta); 0];
    v_pqw = sqrt(mu / p) * [-sin(ta); e + cos(ta); 0];
    
    % Rotation matrix from Perifocal to Inertial (3-1-3)
    cw = cos(w); sw = sin(w);
    co = cos(omg); so = sin(omg);
    ci = cos(i); si = sin(i);
    
    Q = [ co*cw - so*sw*ci, -co*sw - so*cw*ci,  so*si;
          so*cw + co*sw*ci, -so*sw + co*cw*ci, -co*si;
          sw*si,             cw*si,             ci ];
          
    r_vec = Q * r_pqw;
    v_vec = Q * v_pqw;
end

function E_current = Keplers(M,e)
    tol = 10^-8;
    E_0 = M;
    for i=1:10000
       % Set initial loop condition
       if (i==1)
           E_current = E_0;
       end

       % Calculate next step using NR
       E_next = E_current - (E_current - e*sin(E_current) - M) / (1 - e * cos(E_current));
       
       % Check if tol is met
       if norm(E_next - E_current) < tol
           %fprintf("%f", E_current * 180/pi)
           E_current = E_next;
           return
       end

       % Set next loop iteration
       E_current = E_next;
    end
    Fail

end