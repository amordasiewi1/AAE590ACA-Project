%% Load and Setup
load('asteroid_laydown.mat'); % asteroid_params: [a e i RAAN w M_deg value]
N_asteroids = size(asteroid_params, 1);
AU = 149597870.7; 

figure('Color', 'w');
hold on; grid on; axis equal;
view(3);

% Pre-allocate arrays for vectorized plotting (more efficient)
X = zeros(N_asteroids, 1);
Y = zeros(N_asteroids, 1);
Z = zeros(N_asteroids, 1);
Values = asteroid_params(:, 7); % Extract the 1-100 values

%% Math Loop
for j = 1:N_asteroids
    a    = asteroid_params(j, 1) * AU;
    e    = asteroid_params(j, 2);
    inc  = deg2rad(asteroid_params(j, 3));
    raan = deg2rad(asteroid_params(j, 4));
    aop  = deg2rad(asteroid_params(j, 5));
    M    = deg2rad(asteroid_params(j, 6));

    % Solve Kepler's Eq
    E = M; 
    for k = 1:5
        E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
    end
    
    nu = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    r = a * (1 - e*cos(E));
    
    % Perifocal to Inertial Rotation
    r_pqw = [r*cos(nu); r*sin(nu); 0];
    R = [cos(raan)*cos(aop)-sin(raan)*sin(aop)*cos(inc), -cos(raan)*sin(aop)-sin(raan)*cos(aop)*cos(inc),  sin(raan)*sin(inc);
         sin(raan)*cos(aop)+cos(raan)*sin(aop)*cos(inc), -sin(raan)*sin(aop)+cos(raan)*cos(aop)*cos(inc), -cos(raan)*sin(inc);
         sin(aop)*sin(inc),                               cos(aop)*sin(inc),                                cos(inc)];
    
    pos = R * r_pqw;
    X(j) = pos(1); Y(j) = pos(2); Z(j) = pos(3);
end

%% The Plot
% Using 'Values' for both Size and Color
scatter3(X, Y, Z, 40, Values, 'filled', 'MarkerFaceAlpha', 0.8);

% Plot Sun
plot3(0,0,0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');

% Force the Colorbar to the 1-100 scale
colormap(turbo); % 'turbo' or 'jet' are good for value gradients
c = colorbar;
clim([1 100]); % Sets the color scale limits explicitly
ylabel(c, 'Asteroid Economic Value (1-100)');

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Asteroid Laydown: Initial Epoch');