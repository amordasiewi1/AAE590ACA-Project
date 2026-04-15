%% Initialize
clc
clear
close all

%% Declare Variables
a_upper = 2.75;
a_lower = 2.5;
i_upper = asind(0.3);
i_lower = asind(0.1);
e_upper = 0.35;
e_lower = 0.1;
RAAN_upper = 30;
RAAN_lower = 0;
AOP_upper = 60;
AOP_lower = 0;
ta_upper = 90;
ta_lower = 0;
N_asteroids = 50;
value_high = 100;
value_low = 1;
asteroid_params = zeros(N_asteroids, 7);
%% Create Asteroids
for j=1:N_asteroids
    a_ast = createRandom(a_upper, a_lower);
    e_ast = createRandom(e_upper, e_lower);
    i_ast = createRandom(i_upper, i_lower);
    RAAN_ast = createRandom(RAAN_upper, RAAN_lower); % FIXED: angle_lower
    w_ast = createRandom(AOP_upper, AOP_lower);
    ta_deg = createRandom(ta_upper, ta_lower); % Generated in degrees
    
    % --- CONVERT TRUE ANOMALY TO MEAN ANOMALY ---
    % 1. Convert to Radians for math
    ta_rad = deg2rad(ta_deg);
    % 2. Get Eccentric Anomaly (E)
    E_rad = 2 * atan(sqrt((1-e_ast)/(1+e_ast)) * tan(ta_rad/2));
    % 3. Get Mean Anomaly (M)
    M_rad = E_rad - e_ast * sin(E_rad);
    % 4. Convert back to Degrees (since your pipeline expects deg and then converts)
    M_deg = mod(rad2deg(M_rad), 360);
    
    value = randi([value_low, value_high], [1,1]);
    
    % Save with Mean Anomaly as the 6th element
    asteroid = [a_ast e_ast i_ast RAAN_ast w_ast M_deg, value];
    asteroid_params(j, :) = asteroid;
end
%% Save
save("asteroid_laydown.mat", "asteroid_params")

%% Declare Functions
function random = createRandom(ub,lb)
    mu = (ub + lb) / 2;
    sigma = (lb - mu) / 3;
    random = (sigma * randn(1, 1)) + mu;
    if (random > ub) || (random < lb)
        random = createRandom(ub,lb);
    end
end