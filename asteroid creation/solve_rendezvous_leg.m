function [dv_total_kms, x_final, tof_days] = solve_rendezvous_leg(x0, x_target_base, t_dep_days, u_max, t_cap, DU, TU, mu_sun)
    %#codegen
    dt = (0.5 * 24 * 3600) / TU; 
    x = x0;
    dv_DU = 0;
    t_elapsed_TU = 0;
    
    n_target = sqrt(1 / x_target_base(1)^3); 
    
    while t_elapsed_TU < t_cap
        t_current_TU = (t_dep_days * 24 * 3600 / TU) + t_elapsed_TU;
        x_target_now = x_target_base;
        x_target_now(6) = x_target_base(6) + n_target * t_current_TU;
        
        a_sc = max(1e-3, x(1));
        e_sc = max(1e-4, min(x(2), 0.99));
        
        E_sc = Keplers(mod(x(6), 2*pi), e_sc);
        ta_sc = 2 * atan(sqrt((1 + e_sc) / (1 - e_sc)) * tan(E_sc / 2));
        
        E_tgt = Keplers(mod(x_target_now(6), 2*pi), x_target_now(2));
        ta_tgt = 2 * atan(sqrt((1 + x_target_now(2)) / (1 - x_target_now(2))) * tan(E_tgt / 2));
        
        [r_sc, v_sc] = KeplerianToCartesian(a_sc, e_sc, x(3), x(4), x(5), ta_sc, 1);
        [r_tgt, v_tgt] = KeplerianToCartesian(x_target_now(1), x_target_now(2), x_target_now(3), x_target_now(4), x_target_now(5), ta_tgt, 1);
        
        dr = norm(r_sc - r_tgt);
        dv = norm(v_sc - v_tgt);
        
        if dr < 0.01 && dv < 0.001
            break;
        end
        
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
     
     K = diag([50, 50, 50, 50, 50, 50]);
     
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
    p = a * (1 - e^2);
    r_mag = p / (1 + e * cos(ta));
    
    r_pqw = [r_mag * cos(ta); r_mag * sin(ta); 0];
    v_pqw = sqrt(mu / p) * [-sin(ta); e + cos(ta); 0];
    
    cw = cos(w); sw = sin(w);
    co = cos(omg); so = sin(omg);
    ci = cos(i); si = sin(i);
    
    Q = [ co*cw - so*sw*ci, -co*sw - so*cw*ci,  so*si;
          so*cw + co*sw*ci, -so*sw + co*cw*ci, -co*si;
          sw*si,             cw*si,             ci ];
          
    r_vec = Q * r_pqw;
    v_vec = Q * v_pqw;
end

function E_current = Keplers(M, e)
    tol = 1e-8;
    E_current = M; % Initialize before the loop
    for iter = 1:10000
       E_next = E_current - (E_current - e*sin(E_current) - M) / (1 - e * cos(E_current));
       
       if abs(E_next - E_current) < tol
           E_current = E_next;
           return;
       end
       E_current = E_next;
    end
    error('Kepler solver failed to converge.');
end