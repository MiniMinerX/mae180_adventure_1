function dydt = one_sat_sim(t, y, mu, J2, Rearth, A2M, CD)
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);

    a_gravity = -mu / r_norm^3 * r;

    % J2 Pertebation
    factor = -3/2 * J2 * mu * Rearth^2 / r_norm^5;
    z = r(3);
    aJ2 = factor * [(1 - 5 * (z^2 / r_norm^2)) * r(1);
                    (1 - 5 * (z^2 / r_norm^2)) * r(2);
                    (3 - 5 * (z^2 / r_norm^2)) * r(3)];

    h = r_norm - Rearth;
    % Density (kg/m^3) and scale height
    [rho, H] = atmosphere(h);
    
    %convert rho to kg/km^3
    rho_conv = rho * 1e9;

    a_drag = -0.5 * A2M*CD*rho_conv*norm(v)*v;

    % Compute total acceleration
    a_total = a_gravity + a_drag + aJ2 ;

    dydt = [v; a_total;];

end