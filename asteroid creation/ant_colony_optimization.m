%% Ant Colony Optimization for Asteroid Tour (Knapsack Approach - Variable Hops)
clc; 
clear; 
close all;

%% 1. Load Data & Set Constraints
load('asteroid_3D_tensor.mat'); % dv_tensor_kms, departure_epochs, mining_values
asteroid_laydown_file = load("asteroid_laydown.mat");
targets = asteroid_laydown_file.asteroid_params;

% PHYSICS CONSTANTS (Required for Animation)
mu_sun = 132712440017.99;       
a_earth = 149597898;            
DU = a_earth; 
TU = sqrt(DU^3 / mu_sun);
max_thrust = 0.05; 
umax_canonical = max_thrust / 1000 / (DU / TU^2); 
num_targets = size(dv_tensor_kms, 1);
num_epochs = length(departure_epochs);

% --- CRITICAL UNIT NORMALIZATION ---
% The targets dataset stores 'a' in AU, and angles in DEGREES. 
% We MUST convert angles to radians for orbital mechanics math.
% Because 'a' is already in AU, it is ALREADY natively in Canonical Distance Units (DU).
targets(:, 3:6) = deg2rad(targets(:, 3:6));

% --- MISSION BUDGET CONSTRAINTS ---
mission_dv_budget = 60;   % km/s total budget for the ENTIRE mission
max_leg_dv = 25;          % km/s limit per single leg
%transfer_days = 200; 
min_stay_days = 20; 

% Clean the tensor
clean_tensor = dv_tensor_kms;
clean_tensor(clean_tensor == 0) = NaN;
clean_tensor(clean_tensor > max_leg_dv) = NaN;

% Clean the TOF tensor using the same NaN mask as the DV tensor
clean_tof_tensor = tof_tensor_days;
clean_tof_tensor(isnan(clean_tensor)) = NaN;

%% 2. ACO Hyperparameters
num_ants = 100;            
num_iterations = 100;     
alpha_aco = 1.0; % Renamed to avoid shadowing MATLAB's alpha() function             
beta = 1.25;               
rho = 0.25;                
Q = 100;                  

% Initialize Pheromone Matrix
tau = ones(num_targets, num_targets) * 100; % Start high to encourage early exploration

% Global tracking
global_best_fitness = 0;
global_best_path = [];
global_best_epochs = [];
global_best_cost = 0;
global_best_profit = 0;
fitness_history = zeros(num_iterations, 1);

disp(['Launching ', num2str(num_ants), ' budget-constrained ants...']);

%% 3. Main ACO Loop (Variable Hops)
h = waitbar(0, 'Colony exploring asteroid belt...');
for iter = 1:num_iterations
    ant_paths = cell(num_ants, 1);
    ant_epochs = cell(num_ants, 1);
    ant_fitness = zeros(num_ants, 1);
    
    for ant = 1:num_ants
        % 1. Initialization: Start at a random asteroid and epoch
        current_ast = randi(num_targets);
        current_epoch_idx = randi(floor(num_epochs/6)); 
        
        path = current_ast;
        epochs = []; % Stores the departure epoch index for each transit leg
        total_dv = 0;
        total_profit = mining_values(current_ast);
        dv_remaining = mission_dv_budget;
        
        % The clock starts at the arrival/start time of the first asteroid
        current_arrival_day = departure_epochs(current_epoch_idx);
        
        while true
            % 2. Find the earliest valid departure time (Arrival + Stay)
            earliest_next_dep = current_arrival_day + min_stay_days;
            valid_start_idx = find(departure_epochs >= earliest_next_dep, 1, 'first');
            
            if isempty(valid_start_idx), break; end 
            
            probabilities = zeros(num_targets, 1);
            best_epochs_for_targets = zeros(num_targets, 1);
            dvs_for_targets = zeros(num_targets, 1);
            
            for next_ast = 1:num_targets
                if ismember(next_ast, path), continue; end
                
                valid_f_epochs = valid_start_idx:num_epochs;
                dv_opts = squeeze(clean_tensor(current_ast, next_ast, valid_f_epochs));
                [min_dv, loc_idx] = min(dv_opts);
                
                if ~isnan(min_dv) && min_dv <= dv_remaining
                    % Get the absolute epoch index from the lookup table
                    abs_epoch_idx = valid_f_epochs(loc_idx);
                    best_epochs_for_targets(next_ast) = abs_epoch_idx;
                    dvs_for_targets(next_ast) = min_dv;
                    
                    % Heuristic: Reward high profit, low DV, and shorter TOF
                    leg_tof = clean_tof_tensor(current_ast, next_ast, abs_epoch_idx);
                    eta = mining_values(next_ast) / (min_dv + 0.1);
                    probabilities(next_ast) = (tau(current_ast, next_ast)^alpha_aco) * (eta^beta);
                end
            end
            
            if sum(probabilities) == 0, break; end 
            
            prob_norm = probabilities / sum(probabilities);
            next_ast = find(cumsum(prob_norm) >= rand(), 1, 'first');
            
            % 3. Extract the chosen leg details
            chosen_epoch_idx = best_epochs_for_targets(next_ast);
            chosen_tof = clean_tof_tensor(current_ast, next_ast, chosen_epoch_idx);
            
            % 4. Update Ant state
            path = [path, next_ast];
            epochs = [epochs, chosen_epoch_idx];
            total_dv = total_dv + dvs_for_targets(next_ast);
            dv_remaining = dv_remaining - dvs_for_targets(next_ast);
            total_profit = total_profit + mining_values(next_ast);
            
            % 5. Advance the clock for the next iteration
            % Arrival Day = Departure Day of the leg + the time spent in transit
            current_arrival_day = departure_epochs(chosen_epoch_idx) + chosen_tof;
            current_ast = next_ast;
        end
        
        ant_fitness(ant) = total_profit;
        ant_paths{ant} = path;
        ant_epochs{ant} = epochs;
        
        if total_profit > global_best_fitness
            global_best_fitness = total_profit;
            global_best_path = path;
            global_best_epochs = epochs;
            global_best_cost = total_dv;
            global_best_profit = total_profit;
        end
    end
    
    tau = (1 - rho) * tau;
    for ant = 1:num_ants
        if ant_fitness(ant) > 0
            p = ant_paths{ant};
            deposit = Q * (ant_fitness(ant) / max(mining_values));
            for i = 1:(length(p)-1)
                tau(p(i), p(i+1)) = tau(p(i), p(i+1)) + deposit;
            end
        end
    end
    fitness_history(iter) = global_best_fitness;
    waitbar(iter/num_iterations, h);
end
close(h);

%% 4. Display Results and Visualization
if isempty(global_best_path)
    error('No valid mission found.');
end
actual_hops = length(global_best_epochs);
fprintf('\n=== ACO MISSION COMPLETE ===\nHops: %d | Profit: %.2f | Cost: %.2f km/s\n', ...
    actual_hops, global_best_profit, global_best_cost);

% Figure 1: Convergence
figure('Color', 'w'); plot(fitness_history, 'LineWidth', 2); grid on;
title('ACO Convergence (Maximize Profit)'); xlabel('Generation'); ylabel('Total Profit');

% Figure 2: Mission Network Graph
figure('Name', 'Optimal Mission Tour', 'Color', 'w', 'Position', [750, 200, 600, 400]);
s_nodes = string(global_best_path(1:end-1));
t_nodes = string(global_best_path(2:end));
G = digraph(s_nodes, t_nodes);
p_graph = plot(G, 'Layout', 'layered', 'MarkerSize', 10, 'LineWidth', 2, 'ArrowSize', 15);
title(sprintf('Variable Hop Tour (%d Hops | Cost: %.1f km/s)', actual_hops, global_best_cost));
node_ids_in_graph = str2double(G.Nodes.Name);
p_graph.NodeCData = mining_values(node_ids_in_graph);
colormap(parula); colorbar;
highlight(p_graph, string(global_best_path(1)), 'Marker', 's', 'MarkerSize', 18);   
highlight(p_graph, string(global_best_path(end)), 'Marker', 'd', 'MarkerSize', 18); 
text(p_graph.XData(1), p_graph.YData(1)+0.1, 'START', 'Color', 'g', 'FontWeight', 'bold');
text(p_graph.XData(end), p_graph.YData(end)-0.1, 'END', 'Color', 'r', 'FontWeight', 'bold');
axis off;

%% 5. Video Generation (Time-Normalized Cinematic View)
fprintf('Generating Cinematic Mission Video...\n');
video_filename = 'asteroid_mission_tour.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 60;
open(v);
fig_vid = figure('Color', 'k', 'Position', [100 100 1280 720]);
hold on; axis equal; 
view(-30, 20); 
camproj('perspective'); 
set(gca, 'Color', 'k', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');

% --- SPEED CONFIGURATION ---
% Adjust this to make the whole video faster or slower.
% 0.5 means 1 frame for every 2 mission days.
video_speed_factor = 0.1; 

% Create Starfield
num_stars = 500;
star_X = (rand(num_stars, 1) - 0.5) * 15;
star_Y = (rand(num_stars, 1) - 0.5) * 15;
star_Z = (rand(num_stars, 1) - 0.5) * 15;
scatter3(star_X, star_Y, star_Z, 2, 'w', 'filled', 'MarkerFaceAlpha', 0.6);

% Sun
plot3(0,0,0, 'yo', 'MarkerSize', 20, 'MarkerFaceColor', '#FDB813'); 
scatter3(0,0,0, 1200, 'y', 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
axis([-3.5 3.5 -3.5 3.5 -1.5 1.5]);
leg_colors = hsv(actual_hops);

% HUD Setup
hud_text = text(0, 0, 0, 'INITIALIZING...', 'Color', '#00FF00', 'FontSize', 14, ...
    'FontName', 'Courier', 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [0.02 0.9 0]);

% Spacecraft marker
sc_marker = plot3(0,0,0, 'w^', 'MarkerSize', 8, 'MarkerFaceColor', 'w');

for leg = 1:actual_hops
    s_ast = global_best_path(leg);
    e_ast = global_best_path(leg+1);
    e_idx = global_best_epochs(leg);
    
    % TIMING CALCS
    t_dep_days = departure_epochs(e_idx);
    leg_tof = clean_tof_tensor(s_ast, e_ast, e_idx);
    t_arr_days = t_dep_days + leg_tof;
    
    t_dep_TU = (t_dep_days * 24 * 3600) / TU;
    t_arr_TU = (t_arr_days * 24 * 3600) / TU;
    
    % CALCULATE EXACT 3D POINTS
    n_start = sqrt(1 / targets(s_ast,1)^3); 
    x_start = targets(s_ast, 1:6)'; 
    x_start(6) = x_start(6) + n_start * t_dep_TU;
    [X1, Y1, Z1] = elementsToCartesian(x_start);
    R1 = [X1; Y1; Z1];
    
    n_end = sqrt(1 / targets(e_ast,1)^3);
    x_end = targets(e_ast, 1:6)';
    x_end(6) = x_end(6) + n_end * t_arr_TU;
    [X2, Y2, Z2] = elementsToCartesian(x_end);
    R2 = [X2; Y2; Z2];
    
    % SLERP GEOMETRY
    mag1 = norm(R1); mag2 = norm(R2);
    dot_val = max(-1, min(1, dot(R1, R2) / (mag1 * mag2)));
    theta = acos(dot_val);
    n1 = R1 / mag1;
    v_perp = (R2/mag2) - dot_val * n1;
    n2_ortho = v_perp / norm(v_perp);
    
    % DYNAMIC STEP CALCULATION (TRANSIT)
    transit_steps = max(10, round(leg_tof * video_speed_factor));
    dt_days = leg_tof / transit_steps;
    trail = animatedline('Color', leg_colors(leg,:), 'LineWidth', 2.5);
    
    % --- 1. TRANSIT PHASE ---
    for s = 1:transit_steps
        frac = s / transit_steps;
        r_mag = mag1 + (mag2 - mag1) * frac;
        n_curr = n1 * cos(theta * frac) + n2_ortho * sin(theta * frac);
        R_curr = r_mag * n_curr;
        
        X = R_curr(1); Y = R_curr(2); Z = R_curr(3);
        addpoints(trail, X, Y, Z);
        set(sc_marker, 'XData', X, 'YData', Y, 'ZData', Z);
        
        camtarget([X, Y, Z]);
        camorbit(0.6, 0, 'data', [0 0 1]); 
        
        curr_day = t_dep_days + (s * dt_days);
        hud_text.String = sprintf('PHASE: TRANSIT\nMISSION DAY: %d\nTARGET: AST %d', round(curr_day), e_ast);
        hud_text.Color = '#00FF00';
        
        drawnow;
        writeVideo(v, getframe(fig_vid));
    end
    
    % Mark the arrival
    plot3(X, Y, Z, 'wo', 'MarkerSize', 8, 'MarkerFaceColor', leg_colors(leg,:));
    text(X, Y, Z+0.1, sprintf(' Ast %d', e_ast), 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
    
    % --- 2. MINING / LOITERING PHASE ---
    if leg < actual_hops
        next_dep_days = departure_epochs(global_best_epochs(leg+1));
        wait_days = next_dep_days - t_arr_days;
        
        if wait_days > 0
            % DYNAMIC STEP CALCULATION (MINING)
            % This ensures mining moves at the same "days per second" as transit
            loiter_steps = max(10, round(wait_days * video_speed_factor));
            
            x_ast_loiter = x_end; 
            dt_wait_days = wait_days / loiter_steps;
            dt_wait_TU = (dt_wait_days * 24 * 3600) / TU;
            
            loiter_trail = animatedline('Color', '#FFD700', 'LineStyle', ':', 'LineWidth', 2);
            hud_text.Color = '#FFD700';
            
            for w = 1:loiter_steps
                [X_l, Y_l, Z_l] = elementsToCartesian(x_ast_loiter);
                addpoints(loiter_trail, X_l, Y_l, Z_l);
                set(sc_marker, 'XData', X_l, 'YData', Y_l, 'ZData', Z_l);
                
                camtarget([X_l, Y_l, Z_l]);
                camorbit(0.6, 0, 'data', [0 0 1]);
                
                curr_day = t_arr_days + (w * dt_wait_days);
                hud_text.String = sprintf('PHASE: MINING\nMISSION DAY: %d\nLOC: AST %d\nDEPARTURE IN: %d DAYS', ...
                    round(curr_day), e_ast, round(wait_days - (w * dt_wait_days)));
                
                x_ast_loiter(6) = x_ast_loiter(6) + n_end * dt_wait_TU;
                drawnow;
                writeVideo(v, getframe(fig_vid));
            end
        end
    else
        % --- 3. END OF MISSION ---
        hud_text.String = sprintf('PHASE: MISSION COMPLETE\nTOTAL ELAPSED DAYS: %d\nFINAL LOC: AST %d', round(curr_day), e_ast);
        hud_text.Color = '#00FFFF';
        for p = 1:120 
            camorbit(0.2, 0, 'data', [0 0 1]);
            camzoom(0.995); 
            drawnow;
            writeVideo(v, getframe(fig_vid));
        end
    end
end
close(v);

%% 4.5 Extended Visualizations & Analytics
fprintf('Generating extended mission analytics...\n');
hop_dvs = zeros(1, actual_hops);
hop_profits = zeros(1, actual_hops);
transit_start_days = zeros(1, actual_hops);
arrival_days = zeros(1, actual_hops);

for leg = 1:actual_hops

    s_ast = global_best_path(leg);
    e_ast = global_best_path(leg+1);
    e_idx = global_best_epochs(leg);
    leg_tof = clean_tof_tensor(s_ast, e_ast, e_idx);

    hop_dvs(leg) = clean_tensor(s_ast, e_ast, e_idx);
    hop_profits(leg) = mining_values(e_ast);
    transit_start_days(leg) = departure_epochs(e_idx);
    arrival_days(leg) = transit_start_days(leg) + leg_tof;
end
start_val = mining_values(global_best_path(1));
cum_dv = cumsum([0, hop_dvs]);
cum_profit = cumsum([start_val, hop_profits]); 

% --- Figure 3: Mission Resource Dashboard ---
figure('Name', 'Mission Resources', 'Color', 'w', 'Position', [100, 100, 800, 400]);
yyaxis left
stairs(0:actual_hops, cum_dv, 'LineWidth', 2, 'Color', '#0072BD'); hold on;
plot([0, actual_hops], [mission_dv_budget, mission_dv_budget], 'r--', 'LineWidth', 2);
ylabel('Cumulative \DeltaV (km/s)');
ylim([0, max(mission_dv_budget * 1.2, max(cum_dv) * 1.2)]);
yyaxis right
plot(0:actual_hops, cum_profit, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
ylabel('Cumulative Profit (Units)');
title('Mission Resource Accumulation');
xlabel('Hop Number (0 = Start Location)');
legend('Fuel Used', 'Fuel Budget', 'Profit', 'Location', 'northwest');
grid on;

% --- Figure 4: Mission Timeline (Gantt Chart) ---
figure('Name', 'Mission Timeline', 'Color', 'w', 'Position', [150, 150, 800, 350]);
hold on;
for leg = 1:actual_hops
    plot([transit_start_days(leg), arrival_days(leg)], [leg, leg], '-', 'Color', '#0072BD', 'LineWidth', 12);
    if leg < actual_hops
        next_dep = transit_start_days(leg+1);
        plot([arrival_days(leg), next_dep], [leg, leg], '-', 'Color', '#77AC30', 'LineWidth', 12);
    else
        plot([arrival_days(leg), arrival_days(leg) + min_stay_days], [leg, leg], '-', 'Color', '#77AC30', 'LineWidth', 12);
    end
end
yticks(1:actual_hops);
yticklabels(strcat('Leg ', string(1:actual_hops)));
xlabel('Mission Elapsed Time (Days)');
title('Spacecraft Operations Timeline');
p1 = plot(NaN,NaN,'-', 'Color', '#0072BD', 'LineWidth', 6);
p2 = plot(NaN,NaN,'-', 'Color', '#77AC30', 'LineWidth', 6);
legend([p1, p2], {'Transit Phase', 'Mining/Loitering Phase'}, 'Location', 'northeast');
grid on; box on; set(gca, 'YDir', 'reverse'); 
ylim([0, actual_hops+1]);

% --- Figure 5: Full Asteroid Field Context ---
figure('Name', 'Asteroid Field Context', 'Color', 'k', 'Position', [200, 200, 800, 600]);
hold on; axis equal; view(3);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
plot3(0,0,0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', '#FDB813'); 
all_X = zeros(num_targets, 1); all_Y = zeros(num_targets, 1); all_Z = zeros(num_targets, 1);
for i = 1:num_targets
    elem = targets(i, 1:6);
    [all_X(i), all_Y(i), all_Z(i)] = elementsToCartesian(elem);
end
scatter3(all_X, all_Y, all_Z, 20, mining_values, 'filled', 'MarkerFaceAlpha', 0.4);
colormap(parula); cb = colorbar; cb.Color = 'w'; cb.Label.String = 'Mining Value';
m_X = all_X(global_best_path); m_Y = all_Y(global_best_path); m_Z = all_Z(global_best_path);
plot3(m_X, m_Y, m_Z, 'w-', 'LineWidth', 2);
scatter3(m_X, m_Y, m_Z, 60, 'w', 'filled'); 
title('Full Asteroid Field with Optimal Path Overlay', 'Color', 'w');

%% 4.6 Corrected Flight Data Analytics (State & Control History)
fprintf('Generating Synchronized Flight Data Analytics...\n');
t_all = [];
x_all = [];
u_all = [];
target_cart = []; % Track target position continuously

% 1. INITIALIZE STARTING STATE ONLY ONCE
s_start = global_best_path(1);
e_idx_start = global_best_epochs(1); 
t_start_mission_TU = (departure_epochs(e_idx_start) * 24 * 3600) / TU;

n_start = sqrt(1 / targets(s_start,1)^3);
x_curr = targets(s_start, 1:6)';
x_curr(6) = x_curr(6) + n_start * t_start_mission_TU;

for leg = 1:actual_hops
    s_ast = global_best_path(leg);
    e_ast = global_best_path(leg+1);
    e_idx = global_best_epochs(leg);
    
    % Use actual TOF from the synchronized lookup table
    leg_tof_days = clean_tof_tensor(s_ast, e_ast, e_idx);
    t_dep_days = departure_epochs(e_idx);
    t_arr_days = t_dep_days + leg_tof_days;
    
    % --- TRANSIT PHASE SIMULATION ---
    target_elem_base = targets(e_ast, 1:6)';
    n_target = sqrt(1 / target_elem_base(1)^3);
    
    t_start_TU = (t_dep_days * 24 * 3600) / TU;
    leg_tof_days = clean_tof_tensor(s_ast, e_ast, e_idx);
    t_end_TU = ((t_dep_days + leg_tof_days) * 24 * 3600) / TU;

    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    ode_func = @(t, x) compute_state_dot(t, x, target_elem_base, n_target, umax_canonical);
    
    [T_out_TU, X_out] = ode45(ode_func, [t_start_TU, t_end_TU], x_curr, options);

    % Convert and store Transit Data
    leg_t = (T_out_TU * TU) / (24 * 3600);
    leg_u = zeros(length(T_out_TU), 1);
    leg_target_pos = zeros(length(T_out_TU), 3);
    
    for i = 1:length(T_out_TU)
        % Propagate target to exact time step for validation
        xf_target = target_elem_base;
        xf_target(6) = target_elem_base(6) + n_target * T_out_TU(i);
        [leg_target_pos(i,1), leg_target_pos(i,2), leg_target_pos(i,3)] = elementsToCartesian(xf_target);
        
        % Reconstruct Control
        [~, u_mag_canonical] = get_derivatives(X_out(i,:)', xf_target, umax_canonical);
        leg_u(i) = u_mag_canonical * 1000 * (DU / TU^2);
    end

    % --- TERMINAL SANITY CHECK ---
    % Convert the final spacecraft state of THIS specific leg to Cartesian
    [sc_x_f, sc_y_f, sc_z_f] = elementsToCartesian(X_out(end,:));
    sc_final_pos = [sc_x_f, sc_y_f, sc_z_f]; 
    
    % The target's final position is the last row of the target history we just built
    target_final_pos = leg_target_pos(end, :);
    
    % Calculate the miss distance
    miss_distance_AU = norm(sc_final_pos - target_final_pos);
    miss_distance_km = miss_distance_AU * 149597870.7;
    
    fprintf('Leg %d to Ast %d: Miss Distance = %.6f AU (%.2f km)\n', ...
        leg, e_ast, miss_distance_AU, miss_distance_km);
    
    % Store these for the final Figure 7 geometry plot
    intercept_points(leg).sc = sc_final_pos;
    intercept_points(leg).target = target_final_pos;
    
    t_all = [t_all; leg_t];
    x_all = [x_all; X_out];
    u_all = [u_all; leg_u];
    target_cart = [target_cart; leg_target_pos];
    
    % Update state for arrival
    x_curr = target_elem_base;
    x_curr(6) = target_elem_base(6) + n_target * t_end_TU;
    
    % --- MINING / COASTING PHASE ---
    if leg < actual_hops
        next_dep_days = departure_epochs(global_best_epochs(leg+1));
        wait_days = next_dep_days - t_arr_days;
        
        if wait_days > 0
            % Generate points for the loitering duration
            wait_t = linspace(t_arr_days, next_dep_days, round(wait_days)+1)';
            wait_t = wait_t(2:end); % Avoid duplicate timestamps at the join
            
            wait_x = zeros(length(wait_t), 6);
            wait_target_pos = zeros(length(wait_t), 3);
            n_curr_can = sqrt(1 / x_curr(1)^3); 
            
            for ws = 1:length(wait_t)
                dt_coast_TU = (wait_t(ws) - t_arr_days) * 24 * 3600 / TU;
                wait_x(ws, :) = x_curr';
                wait_x(ws, 6) = x_curr(6) + n_curr_can * dt_coast_TU;
                
                % Target position (Spacecraft is sitting on the asteroid)
                xf_target = target_elem_base;
                xf_target(6) = target_elem_base(6) + n_target * ( (wait_t(ws)*24*3600)/TU );
                [wait_target_pos(ws,1), wait_target_pos(ws,2), wait_target_pos(ws,3)] = elementsToCartesian(xf_target);
            end
            
            t_all = [t_all; wait_t];
            x_all = [x_all; wait_x];
            u_all = [u_all; zeros(length(wait_t), 1)];
            target_cart = [target_cart; wait_target_pos];
            
            % Update for next leg departure
            x_curr = wait_x(end, :)';
        end
    end
end

% --- GENERATE FULL CARTESIAN HISTORY ---
sc_cart = zeros(length(t_all), 3);
for k = 1:length(t_all)
    [sc_cart(k,1), sc_cart(k,2), sc_cart(k,3)] = elementsToCartesian(x_all(k,:));
end

% --- Figure 6: Control Input History ---
figure('Name', 'Control Input History', 'Color', 'w', 'Position', [100, 100, 900, 350]);
hold on; grid on;
y_max = max_thrust * 1.2;
for leg = 1:actual_hops-1
     s_ast = global_best_path(leg);
    e_ast = global_best_path(leg+1);
    e_idx = global_best_epochs(leg);
    leg_tof = clean_tof_tensor(s_ast, e_ast, e_idx);

    t_dep_days = departure_epochs(global_best_epochs(leg));
    t_arr_days = t_dep_days + leg_tof;
    t_arr = t_arr_days;
    t_next = departure_epochs(global_best_epochs(leg+1));
    fill([t_arr t_next t_next t_arr], [0 0 y_max y_max], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
plot(t_all, u_all, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Applied Thrust');
plot([t_all(1), t_all(end)], [max_thrust, max_thrust], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Max Thrust Limit');
xlabel('Mission Elapsed Time (Days)');
ylabel('Thrust Magnitude (m/s^2)');
title('Spacecraft Control Effort Over Time (Green = Mining/Coasting)');
legend('Location', 'northeast');
ylim([0, y_max]);
xlim([t_all(1), t_all(end)]);

% --- CORRECTED Figure 7: Cartesian Intersection Validation ---
% --- CORRECTED Figure 7: Intercept Geometry ---
figure('Name', 'Intercept Geometry', 'Color', 'w'); hold on; axis equal; view(3);
grid on;
for leg = 1:actual_hops
    % Plot Spacecraft (Circle) and Target (X)
    plot3(intercept_points(leg).sc(1), intercept_points(leg).sc(2), intercept_points(leg).sc(3), ...
          'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot3(intercept_points(leg).target(1), intercept_points(leg).target(2), intercept_points(leg).target(3), ...
          'kx', 'MarkerSize', 12, 'LineWidth', 1.5);
    
    % If there's a gap, this line will show it
    line([intercept_points(leg).sc(1) intercept_points(leg).target(1)], ...
         [intercept_points(leg).sc(2) intercept_points(leg).target(2)], ...
         [intercept_points(leg).sc(3) intercept_points(leg).target(3)], 'Color', 'r');
end
title('Spatial Overlap at Intercept Moments');
xlabel('X (AU)'); ylabel('Y (AU)'); zlabel('Z (AU)');
legend('Spacecraft Arrival', 'Target Asteroid', 'Miss Distance Vector');


save('aco_results.mat', 'global_best_path', 'global_best_epochs', 'actual_hops');
%% Helper Functions
function dx_dt = compute_state_dot(t_TU, x_curr, target_elem_base, n_target, u_max)
    % 1. Propagate the target's Mean Anomaly to the solver's current time
    xf_target = target_elem_base;
    xf_target(6) = xf_target(6) + n_target * t_TU; 
    
    % 2. Get the state derivatives from your existing controller
    % (We drop the second output, u_mag, because ode45 only wants dx_dt)
    [dx_dt, ~] = get_derivatives(x_curr, xf_target, u_max);
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

function E = Keplers(M, e)
    E = M; 
    for k = 1:10
        E = E - (E - e*sin(E) - M)/(1 - e*cos(E));
    end
end