%% Propagate Optimal Path & Visualize Control History (RK4 Version)
clc; clear; close all;

%% 1. Load Data
try
    load('aco_results.mat'); 
catch
    error('Could not find aco_results.mat. Please save the optimal path variables from your ACO script first.');
end

load('asteroid_3D_tensor.mat'); % dv_tensor_kms, departure_epochs, mining_values, tof_tensor_days
asteroid_laydown_file = load("asteroid_laydown.mat");
targets = asteroid_laydown_file.asteroid_params;

%% 2. Setup Physics Constants
mu_sun = 132712440017.99;       
a_earth = 149597898;            
DU = a_earth; 
TU = sqrt(DU^3 / mu_sun);
max_thrust = 0.05; 
umax_canonical = max_thrust / 1000 / (DU / TU^2); 

% Normalize angles
targets(:, 3:6) = deg2rad(targets(:, 3:6)); 

%% 3. Flight Data Extraction & Propagation
fprintf('Propagating Flight Data using original RK4 logic...\n');
t_all = [];
x_all = [];
u_all = [];
target_cart = []; 

% INITIALIZE STARTING STATE ONLY ONCE
s_start = global_best_path(1);
e_idx_start = global_best_epochs(1); 
t_start_mission_TU = (departure_epochs(e_idx_start) * 24 * 3600) / TU;

n_start = sqrt(1 / targets(s_start,1)^3);
x_curr = targets(s_start, 1:6)';
x_curr(6) = x_curr(6) + n_start * t_start_mission_TU;

intercept_points = struct();

for leg = 1:actual_hops
    s_ast = global_best_path(leg);
    e_ast = global_best_path(leg+1);
    e_idx = global_best_epochs(leg);
    
    % Timing
    leg_tof_days = tof_tensor_days(s_ast, e_ast, e_idx);
    t_dep_days = departure_epochs(e_idx);
    t_arr_days = t_dep_days + leg_tof_days;
    
    target_elem_base = targets(e_ast, 1:6)';
    
    % --- DYNAMIC RK4 INTEGRATION (Exact logic from create_lookup_table.m) ---
    [T_out_TU, X_out, leg_u, leg_target_pos] = simulate_rk4_leg(...
        x_curr, target_elem_base, t_dep_days, leg_tof_days, umax_canonical, DU, TU);

    % Convert Time to Days for tracking
    leg_t = (T_out_TU * TU) / (24 * 3600);

    % Terminal Sanity Check
    [sc_x_f, sc_y_f, sc_z_f] = elementsToCartesian(X_out(end,:));
    sc_final_pos = [sc_x_f, sc_y_f, sc_z_f]; 
    target_final_pos = leg_target_pos(end, :);
    
    miss_distance_AU = norm(sc_final_pos - target_final_pos);
    miss_distance_km = miss_distance_AU * DU;
    
    fprintf('Leg %d to Ast %d: Miss Distance = %.6f AU (%.2f km)\n', ...
        leg, e_ast, miss_distance_AU, miss_distance_km);
    
    intercept_points(leg).sc = sc_final_pos;
    intercept_points(leg).target = target_final_pos;
    
    t_all = [t_all; leg_t];
    x_all = [x_all; X_out];
    u_all = [u_all; leg_u];
    target_cart = [target_cart; leg_target_pos];
    
    x_curr = X_out(end, :)';
    
    % Mining / Coasting Phase
    if leg < actual_hops
        next_dep_days = departure_epochs(global_best_epochs(leg+1));
        wait_days = next_dep_days - t_arr_days;
        
        if wait_days > 0
            wait_t = linspace(t_arr_days, next_dep_days, round(wait_days)+1)';
            wait_t = wait_t(2:end); 
            
            wait_x = zeros(length(wait_t), 6);
            wait_target_pos = zeros(length(wait_t), 3);
            n_curr_can = sqrt(1 / x_curr(1)^3); 
            n_target = sqrt(1 / target_elem_base(1)^3);
            
            for ws = 1:length(wait_t)
                dt_coast_TU = (wait_t(ws) - t_arr_days) * 24 * 3600 / TU;
                wait_x(ws, :) = x_curr';
                wait_x(ws, 6) = x_curr(6) + n_curr_can * dt_coast_TU;
                
                xf_target = target_elem_base;
                xf_target(6) = target_elem_base(6) + n_target * ( (wait_t(ws)*24*3600)/TU );
                [wait_target_pos(ws,1), wait_target_pos(ws,2), wait_target_pos(ws,3)] = elementsToCartesian(xf_target);
            end
            
            t_all = [t_all; wait_t];
            x_all = [x_all; wait_x];
            u_all = [u_all; zeros(length(wait_t), 1)];
            target_cart = [target_cart; wait_target_pos];
            
            x_curr = wait_x(end, :)';
        end
    end
end

%% 4. Visualization
% Figure 1: Control Input History
figure('Name', 'Control Input History', 'Color', 'w', 'Position', [100, 100, 900, 350]);
hold on; grid on;
y_max = max_thrust * 1.2;
for leg = 1:actual_hops-1
    s_ast = global_best_path(leg);
    e_ast = global_best_path(leg+1);
    e_idx = global_best_epochs(leg);
    leg_tof = tof_tensor_days(s_ast, e_ast, e_idx);

    t_dep_days = departure_epochs(global_best_epochs(leg));
    t_arr_days = t_dep_days + leg_tof;
    t_next = departure_epochs(global_best_epochs(leg+1));
    fill([t_arr_days t_next t_next t_arr_days], [0 0 y_max y_max], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
plot(t_all, u_all, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Applied Thrust');
plot([t_all(1), t_all(end)], [max_thrust, max_thrust], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Max Thrust Limit');
xlabel('Mission Elapsed Time (Days)');
ylabel('Thrust Magnitude (m/s^2)');
title('Spacecraft Control Effort Over Time (Green = Mining/Coasting)');
legend('Location', 'northeast');
ylim([0, y_max]);
xlim([t_all(1), t_all(end)]);

% Figure 2: Intercept Geometry Validation
figure('Name', 'Intercept Geometry', 'Color', 'w'); hold on; axis equal; view(3);
grid on;
for leg = 1:actual_hops
    plot3(intercept_points(leg).sc(1), intercept_points(leg).sc(2), intercept_points(leg).sc(3), ...
          'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot3(intercept_points(leg).target(1), intercept_points(leg).target(2), intercept_points(leg).target(3), ...
          'kx', 'MarkerSize', 12, 'LineWidth', 1.5);
    
    line([intercept_points(leg).sc(1) intercept_points(leg).target(1)], ...
         [intercept_points(leg).sc(2) intercept_points(leg).target(2)], ...
         [intercept_points(leg).sc(3) intercept_points(leg).target(3)], 'Color', 'r');
end
title('Spatial Overlap at Intercept Moments');
xlabel('X (AU)'); ylabel('Y (AU)'); zlabel('Z (AU)');
legend('Spacecraft Arrival', 'Target Asteroid', 'Miss Distance Vector');

%% Helper Functions

function [T_history_TU, X_history, U_history, Target_history] = simulate_rk4_leg(x0, target_elem_base, t_dep_days, leg_tof_days, u_max, DU, TU)
    dt_TU = (0.5 * 24 * 3600) / TU; % 12-hour steps exactly like lookup table
    x = x0;
    t_elapsed_TU = 0;
    n_target = sqrt(1 / target_elem_base(1)^3);
    
    % Preallocate
    num_steps_est = ceil(leg_tof_days / 0.5) + 10;
    T_history_TU = zeros(num_steps_est, 1);
    X_history = zeros(num_steps_est, 6);
    U_history = zeros(num_steps_est, 1);
    Target_history = zeros(num_steps_est, 3);
    
    idx = 1;
    while true
        t_current_TU = (t_dep_days * 24 * 3600 / TU) + t_elapsed_TU;
        
        % Target position
        x_target_now = target_elem_base;
        x_target_now(6) = target_elem_base(6) + n_target * t_current_TU;
        
        % Check Convergence (Original Logic)
        err = x - x_target_now;
        err(4:6) = atan2(sin(err(4:6)), cos(err(4:6)));
        
        if abs(err(1)) < 1e-3 && abs(err(6)) < 0.02 
            break;
        end
        
        % Failsafe slightly above known TOF in case of rounding diffs
        if t_elapsed_TU > ((leg_tof_days + 5) * 24 * 3600 / TU) 
            break;
        end
        
        % Log State Before Step
        T_history_TU(idx) = t_current_TU;
        X_history(idx, :) = x';
        [tX, tY, tZ] = elementsToCartesian(x_target_now);
        Target_history(idx, :) = [tX, tY, tZ];
        
        % RK4 Step
        [dx1, u1] = get_derivatives(x, x_target_now, u_max);
        [dx2, u2] = get_derivatives(x + 0.5*dt_TU*dx1, x_target_now, u_max);
        [dx3, u3] = get_derivatives(x + 0.5*dt_TU*dx2, x_target_now, u_max);
        [dx4, u4] = get_derivatives(x + dt_TU*dx3, x_target_now, u_max);
        
        x = x + (dt_TU/6)*(dx1 + 2*dx2 + 2*dx3 + dx4);
        
        % Store thrust magnitude applied during this step
        U_history(idx) = norm(u1) * 1000 * (DU / TU^2); 
        
        t_elapsed_TU = t_elapsed_TU + dt_TU;
        idx = idx + 1;
    end
    
    % Log Final State
    T_history_TU(idx) = (t_dep_days * 24 * 3600 / TU) + t_elapsed_TU;
    X_history(idx, :) = x';
    
    x_target_now = target_elem_base;
    x_target_now(6) = target_elem_base(6) + n_target * T_history_TU(idx);
    [tX, tY, tZ] = elementsToCartesian(x_target_now);
    Target_history(idx, :) = [tX, tY, tZ];
    
    [~, u_final] = get_derivatives(x, x_target_now, u_max);
    U_history(idx) = norm(u_final) * 1000 * (DU / TU^2);
    
    % Trim arrays to actual size
    T_history_TU = T_history_TU(1:idx);
    X_history = X_history(1:idx, :);
    U_history = U_history(1:idx);
    Target_history = Target_history(1:idx, :);
end

function [X, Y, Z] = elementsToCartesian(elem)
    a = elem(1); e = elem(2); i = elem(3); omg = elem(4); w = elem(5); M = elem(6);
    E = M; for k=1:5, E = E - (E - e*sin(E) - M)/(1 - e*cos(E)); end
    r = a * (1 - e*cos(E));
    ta = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2));
    X = r*(cos(omg)*cos(w+ta) - sin(omg)*sin(w+ta)*cos(i));
    Y = r*(sin(omg)*cos(w+ta) + cos(omg)*sin(w+ta)*cos(i));
    Z = r*(sin(w+ta)*sin(i));
end

function [x_dot, u_mag] = get_derivatives(x, x_target, u_max)
     a = max(1e-3, x(1));                       
     e = max(1e-4, min(x(2), 0.90));           
     i = max(1e-4, min(x(3), pi - 1e-4));       
     omg = x(4); w = x(5); M = x(6);
     
     dx = zeros(6,1);
     dx(1:3) = x(1:3) - x_target(1:3);
     dx(4:6) = atan2(sin(x(4:6)-x_target(4:6)), cos(x(4:6)-x_target(4:6)));
     
     K = diag([10000, 5000, 2000, 1000, 1000, 5000]);
     
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

function E = Keplers(M, e)
    E = M; 
    for k = 1:10
        E = E - (E - e*sin(E) - M)/(1 - e*cos(E));
    end
end