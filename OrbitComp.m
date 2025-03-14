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
    
    % Downsample the orbit for better performance
    step = 500; % Adjust to balance performance & smoothness
    Y_downsampled = Y(1:step:end, :);
    numPoints = size(Y_downsampled, 1);
    
    % Generate a color gradient
    cmap = jet(numPoints);
    
    figure;
    hold on;
    
    % Plot the rainbow-colored orbit without legend spam
    for i = 1:numPoints-1
        plot3(Y_downsampled(i:i+1,1), Y_downsampled(i:i+1,2), Y_downsampled(i:i+1,3), ...
            'Color', cmap(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Prevent legend spam
    end
    
    % Plot Earth
    [X, Y2, Z] = sphere(50);
    X = X * Rearth;
    Y2 = Y2 * Rearth;
    Z = Z * Rearth;
    surf(X, Y2, Z, 'FaceColor', 'c', 'EdgeColor', 'none'); 
    alpha(0.6);
    
    % Plot r0 and r_end with distinct markers (only these go in legend)
    scatter3(r0(1), r0(2), r0(3), 150, 'r', 'filled', 'o', 'DisplayName', 'Initial Position');
    scatter3(r_end(1), r_end(2), r_end(3), 150, 'g', 'filled', 'o', 'DisplayName', 'Final Position');
    
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


