% general rotation matrices - radians
% Alessandro Mordasiewicz - 10/29/2025
% Source: STM script eqns. 1.2 - 1.4

% about 1st axis (1.2)
R1 = @(x) [ 1  0      0;
             0  cos(x) sin(x);
             0 -sin(x) cos(x)];

% about 2nd axis (1.3)
R2 = @(x) [cos(x) 0 -sin(x);
            0      1  0;
            sin(x) 0  cos(x)];

% about 3rd axis (1.4)
R3 = @(x) [ cos(x) sin(x) 0;
            -sin(x) cos(x) 0;
             0      0      1];