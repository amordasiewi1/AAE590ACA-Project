load('asteroid_3D_tensor.mat');

% 1. Check for NaNs (Did the GVE singularities break the math?)
num_nans = sum(isnan(dv_tensor_kms), 'all');
fprintf('Number of NaNs: %d\n', num_nans);

% 2. Check for controller saturation
% 0.1 m/s^2 for 50 days is exactly 432 km/s. 
max_dv = max(dv_tensor_kms, [], 'all');
fprintf('Maximum Delta-V: %.2f km/s\n', max_dv);

% 3. Look at the distribution
% Exclude the 0s on the diagonal (start == end)
valid_dvs = dv_tensor_kms(dv_tensor_kms > 0);
fprintf('Mean Delta-V: %.2f km/s\n', mean(valid_dvs));
fprintf('Median Delta-V: %.2f km/s\n', median(valid_dvs));

% Find the minimum cost across all epochs for each asteroid pair
min_dv_across_time = min(dv_tensor_kms, [], 3);

% Set diagonal to NaN so the color scale isn't ruined by 0s
min_dv_across_time(logical(eye(size(min_dv_across_time)))) = NaN;

figure('Name', 'Best Transfer Costs', 'Color', 'w');
h = heatmap(min_dv_across_time);
h.Title = 'Minimum \Delta V Across All Epochs (km/s)';
h.XLabel = 'Target Asteroid ID';
h.YLabel = 'Departure Asteroid ID';
colormap(parula); % Parula is good for highlighting extremes    

% Find the global minimum (ignoring the diagonal zeros)
dv_search = dv_tensor_kms;
dv_search(dv_search == 0) = NaN; 
[best_dv, idx] = min(dv_search, [], 'all', 'linear');
[best_start, best_end, best_epoch] = ind2sub(size(dv_search), idx);

fprintf('\n--- BEST TRANSFER ---\n');
fprintf('Departure Asteroid: %d\n', best_start);
fprintf('Target Asteroid: %d\n', best_end);
fprintf('Epoch: %d (Day %d)\n', best_epoch, departure_epochs(best_epoch));
fprintf('Cost: %.2f km/s\n', best_dv);