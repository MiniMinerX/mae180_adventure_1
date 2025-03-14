function [r0, v0, oe0, r, v, oe] = OrbitComp(lat, lst, alt, ra, dec, JD, JD_prop)
    
    %CONSTANTS
    mu= 398600.4354; %km^3/s^2
    Rearth= 6378.1366; %km
    smallomegaearth= 7.2921159E-5; %rad/s
    Areatomassratio= 0.0123; %m^2/kg
    JDday2secondconversionunit = 86400; 
    J2 = 0.0010826267; %unitless
    CD = 1.28; % drag coefficient corrected from annocement 
    we = 7.2921159E-5;% Earth's inertial rotation rate

    % UCSD coords
    lot=deg2rad(-117.2336137);
    
    % convert JD to seconds and units to radians
    JD = JD * JDday2secondconversionunit;
    JD_prop = JD_prop * JDday2secondconversionunit;
    lat = deg2rad(lat);
    lst = deg2rad(lst);
    ra = deg2rad(ra);
    dec = deg2rad(dec);

    [r0,v0] = gauss(lat, lst, alt, ra, dec, JD, JD_prop);
    
    % 10. Compute initial orbital elements from r0 and v0
    oe0 = ComputeOrbitalElements(r0, v0, mu);
    
    % Initial state vector
    y0 = [r0';v0'];
    
    % propogation span
    JD_init = JD(2);
    JD_end = JD_prop;
    t_prop = JD_end - JD_init;
    tspan = [0, t_prop];
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    A2M_convert = Areatomassratio / 1000000;  % m^2/kg to km^2/kg
    [t, Y] = ode45(@(t, y) one_sat_sim(t, y, mu, J2, Rearth, A2M_convert,CD), tspan, y0, options);
    
    r_end = Y(end,1:3);
    v_end = Y(end,4:6);

    oe = ComputeOrbitalElements(r_end, v_end, mu);

    figure;
    hold on;
    
     % Downsample for smoother plotting
    step = 500;
    Y_down = Y(1:step:end, :);
    numPoints = size(Y_down, 1);

    % Extract position coordinates
    x = Y_down(:, 1);
    y = Y_down(:, 2);
    z = Y_down(:, 3);

    % Create a color vector from 0 to 1 across the orbit
    c = linspace(0, 1, numPoints);

    % Prepare for the "surface trick" to draw a continuous line
    X = repmat(x', 2, 1);
    Y_mat = repmat(y', 2, 1);
    Z = repmat(z', 2, 1);
    C = repmat(c, 2, 1);

    % Plot
    figure; 
    hold on;

    % --- 1) Draw Earth first with constant color data so it doesn’t affect orbit colormap ---
    [Xs, Ys, Zs] = sphere(50);
    Xs = Xs * Rearth;
    Ys = Ys * Rearth;
    Zs = Zs * Rearth;
    % Assign a solid color (c = [0.2 0.6 1]) instead of colormap-based shading
    surf(Xs, Ys, Zs, 'FaceColor', [0.2 0.6 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    
    % --- 2) Draw the orbit line with its own color data ---
    hOrbit = surface(X, Y_mat, Z, C, ...
        'FaceColor', 'none', ...
        'EdgeColor', 'interp', ...
        'LineWidth', 1.5);

    % Use a rainbow-style colormap
    colormap(jet);

    % Force color limits to match our 0..1 data
    caxis([0 1]);

    % Optionally, display a colorbar to illustrate the gradient
    colorbar;
    
    % Plot initial and final positions
    scatter3(r0(1), r0(2), r0(3), 150, 'r', 'filled', 'o', 'DisplayName', 'Initial Position');
    scatter3(r_end(1), r_end(2), r_end(3), 150, 'g', 'filled', 'o', 'DisplayName', 'Final Position');

    % Plot velocity vectors
    quiver3(r0(1), r0(2), r0(3), v0(1), v0(2), v0(3), 1000, 'r', 'LineWidth', 2, 'DisplayName', 'Initial Velocity');
    quiver3(r_end(1), r_end(2), r_end(3), v_end(1), v_end(2), v_end(3), 1000, 'g', 'LineWidth', 2, 'DisplayName', 'Final Velocity');
    
    % Labels and Formatting
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Rainbow Satellite Orbit with J2 Perturbation and Atmospheric Drag');
    grid on;
    axis equal;
    legend('Location', 'best'); % Only shows r0 and r_end
    view(3);
    hold off;

    r = r_end;
    v = v_end;
end


