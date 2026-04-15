%% Multi-Leg Asteroid Tour Solver
clc; clear;

% 1. Mission Constraints & Parameters
load('asteroid_3D_tensor.mat'); % Loads dv_tensor_kms, departure_epochs, mining_values

num_targets = size(dv_tensor_kms, 1);
num_epochs = length(departure_epochs);
transfer_days = 50; 
min_stay_days = 10; % Minimum time spent at the asteroid before leaving
max_leg_dv = 15;    % km/s LIMIT: Ignore any single hop that costs more than this
target_hops = 3;    % Number of transfers in the mission (e.g., 3 hops visits 4 asteroids)

% 2. Clean the Tensor
% Replace diagonal (0s) and impossible transfers (saturation) with NaN so the solver ignores them
clean_tensor = dv_tensor_kms;
clean_tensor(clean_tensor == 0) = NaN;
clean_tensor(clean_tensor > max_leg_dv) = NaN;

% 3. Initialize Global Search Trackers
global best_mission_cost best_mission_path best_mission_epochs;
best_mission_cost = inf;
best_mission_path = [];
best_mission_epochs = [];

disp(['Searching for optimal ', num2str(target_hops), '-hop mission...']);

% 4. Kick off the search from every asteroid at every valid starting epoch
for start_ast = 1:num_targets
    % Only start in the first half of the timeline to leave room for the mission
    for start_epoch = 1:floor(num_epochs/2) 
        path = start_ast;
        epochs = start_epoch;
        cost = 0;
        
        search_branch(start_ast, start_epoch, path, epochs, cost, target_hops, clean_tensor, departure_epochs, transfer_days, min_stay_days);
    end
end

% 5. Print the Results
if isinf(best_mission_cost)
    fprintf('\nNo valid %d-hop mission found under %d km/s per leg.\n', target_hops, max_leg_dv);
    fprintf('Try increasing max_leg_dv, decreasing min_stay_days, or doing fewer hops.\n');
else
    fprintf('\n=== OPTIMAL MISSION ITINERARY ===\n');
    fprintf('Total Delta-V Cost: %.2f km/s\n\n', best_mission_cost);
    
    for i = 1:target_hops
        a_dep = best_mission_path(i);
        a_arr = best_mission_path(i+1);
        e_dep = best_mission_epochs(i);
        dv = clean_tensor(a_dep, a_arr, e_dep);
        
        day_dep = departure_epochs(e_dep);
        day_arr = day_dep + transfer_days;
        
        fprintf('LEG %d: Asteroid %d -> Asteroid %d\n', i, a_dep, a_arr);
        fprintf('   Depart: Day %d (Epoch %d)\n', day_dep, e_dep);
        fprintf('   Arrive: Day %d\n', day_arr);
        fprintf('   Cost:   %.2f km/s\n', dv);
        if i < target_hops
            wait_time = departure_epochs(best_mission_epochs(i+1)) - day_arr;
            fprintf('   [ Stay at Asteroid %d for %d days ]\n\n', a_arr, wait_time);
        end
    end
    
    % Optional: Calculate total profit of visited asteroids
    total_value = sum(mining_values(unique(best_mission_path)));
    fprintf('\nTotal Asteroid Value Mined: %.2f units\n', total_value);
    fprintf('=================================\n');
end


%% Recursive Depth-First Search Function
function search_branch(current_ast, current_epoch, current_path, current_epochs, current_cost, hops_left, tensor, epoch_list, t_transfer, min_stay)
    global best_mission_cost best_mission_path best_mission_epochs;
    
    % Base Case: Mission Complete
    if hops_left == 0
        if current_cost < best_mission_cost
            best_mission_cost = current_cost;
            best_mission_path = current_path;
            best_mission_epochs = current_epochs;
        end
        return;
    end
    
    % Pruning: Stop searching this branch if it's already worse than our best found so far
    if current_cost >= best_mission_cost
        return; 
    end
    
    % Calculate earliest allowed departure for the NEXT leg
    arrival_day = epoch_list(current_epoch) + t_transfer;
    earliest_next_departure_day = arrival_day + min_stay;
    
    % Find the epoch index that matches or exceeds this day
    next_valid_epoch = find(epoch_list >= earliest_next_departure_day, 1, 'first');
    
    % If there's no time left in the simulation, this branch dies
    if isempty(next_valid_epoch)
        return;
    end
    
    num_targets = size(tensor, 1);
    num_epochs = length(epoch_list);
    
    % Try jumping to every possible next asteroid...
    for next_ast = 1:num_targets
        if next_ast == current_ast || ismember(next_ast, current_path)
            continue; % Don't go to the same asteroid, and don't revisit old ones
        end
        
        % ...at every valid future epoch
        for next_epoch = next_valid_epoch:num_epochs
            
            hop_dv = tensor(current_ast, next_ast, next_epoch);
            
            % If it's a valid transfer (not NaN), take the jump!
            if ~isnan(hop_dv)
                new_cost = current_cost + hop_dv;
                new_path = [current_path, next_ast];
                new_epochs = [current_epochs, next_epoch];
                
                % Recurse deeper into the timeline
                search_branch(next_ast, next_epoch, new_path, new_epochs, new_cost, hops_left - 1, tensor, epoch_list, t_transfer, min_stay);
            end
        end
    end
end