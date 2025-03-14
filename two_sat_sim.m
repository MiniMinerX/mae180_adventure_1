function [dydt, d_min_out] = two_sat_sim(t, y, mu, J2, Rearth, A2M, CD)
    global min_distance;
    global min_distance_time;
    global min_distance_r_s1;
    global min_distance_r_s2;

    y_sat_1 = y(1:6);
    y_sat_2 = y(7:12);

    r_s1 = y_sat_1(1:3);
    v_s1 = y_sat_1(4:6);
    r_norm_s1 = norm(r_s1);

    a_gravity_s1 = -mu / r_norm_s1^3 * r_s1;

    % J2 Pertebation
    factor_s1 = -3/2 * J2 * mu * Rearth^2 / r_norm_s1^5;
    z_s1 = r_s1(3);
    aJ2_s1 = factor_s1 * [(1 - 5 * (z_s1^2 / r_norm_s1^2)) * r_s1(1);
                    (1 - 5 * (z_s1^2 / r_norm_s1^2)) * r_s1(2);
                    (3 - 5 * (z_s1^2 / r_norm_s1^2)) * r_s1(3)];

    h_s1 = r_norm_s1 - Rearth;
    % Density (kg/m^3) and scale height
    [rho_s1, H_s1] = atmosphere(h_s1);
    
    %convert rho to kg/km^3
    rho_conv_s1 = rho_s1 * 1e9;

    a_drag_s1 = -0.5 * A2M*CD*rho_conv_s1*norm(v_s1)*v_s1;

    % Compute total acceleration
    a_total_s1 = a_gravity_s1 + a_drag_s1 + aJ2_s1;

    dydt_s1 = [v_s1; a_total_s1;];

    r_s2 = y_sat_2(1:3);
    v_s2 = y_sat_2(4:6);
    r_norm_s2 = norm(r_s2);

    a_gravity_s2 = -mu / r_norm_s2^3 * r_s2;

    % J2 Pertebation
    factor_s2 = -3/2 * J2 * mu * Rearth^2 / r_norm_s2^5;
    z_s2 = r_s2(3);
    aJ2_s2 = factor_s2 * [(1 - 5 * (z_s2^2 / r_norm_s2^2)) * r_s2(1);
                    (1 - 5 * (z_s2^2 / r_norm_s2^2)) * r_s2(2);
                    (3 - 5 * (z_s2^2 / r_norm_s2^2)) * r_s2(3)];

    h_s2 = r_norm_s2 - Rearth;
    % Density (kg/m^3) and scale height
    [rho_s2, H_s2] = atmosphere(h_s2);
    
    %convert rho to kg/km^3
    rho_conv_s2 = rho_s2 * 1e9;

    a_drag_s2 = -0.5 * A2M*CD*rho_conv_s2*norm(v_s2)*v_s2;

    % Compute total acceleration
    a_total_s2 = a_gravity_s2 + a_drag_s2 + aJ2_s2;

    dydt_s2 = [v_s2; a_total_s2;];

    % Compute relative distance between satellites
    d = norm(r_s1 - r_s2);
    if d < min_distance
         min_distance = d;
         min_distance_time = t;
         min_distance_r_s1 = r_s1;
         min_distance_r_s2 = r_s2;


         % Collision avoidance: if within threshold and burn not yet applied, apply Δv.
         %if (d < collision_threshold) && ~burn_executed
              % Impulsive burn on satellite 2 (example: small change in velocity)
              %v2 = v2 + delta_v;
          %    burn_executed = true;
         %     fprintf('Collision avoidance burn executed at t = %.2f s; new sat2 velocity applied.\n', t);
         %end
    end

    % Assemble the derivative of the state vector.
    dydt = [dydt_s1; dydt_s2];

end