function oe = ComputeOrbitalElements(r, v, mu)
    norm_r = norm(r);
    norm_v = norm(v);
    
    % Specific mechanical energy
    epsilon = 0.5 * norm_v^2 - mu / norm_r;
    
    % Semimajor axis
    a = -mu / (2 * epsilon);
    
    % Angular momentum vector h
    h = cross(r, v);
    h_mag = norm(h);
    h_hat = h/h_mag;
    
    % Inclination
    i = acosd(h(3) / h_mag);
    
    % eccentricity
    e_vec = (1/mu) * cross(v, h) - (r / norm_r);
    e_mag = norm(e_vec);
    
    % Node Vector (n)
    z_hat = [0; 0; 1];
    n = cross(z_hat, h);
    n_mag = norm(n);
    
    % Longitude of ascending node (omega)
    if i == 0 || i == 180
       omega = []
    else
       % We check sign of Node vector y value
       % to adjust for acosd quadrant
       % note tan(omega) = n_hat_y / n_hat_x, so tan=sin/cos
       % lets us use the cosine argument
       if n(2) >= 0
           omega = acosd(n(1) / n_mag);
       else
           omega = 360 - acosd(n(1) / n_mag);
       end
    end
    
    % Argument of perigee (w) and related angles
    if e_mag > 1e-10 % elliptical orbit
      if i==0 || i==180 % Equatorial
        % Longitude of perigee (ϖ)
        if e_vec(2) >= 0
            w = acosd(e_vec(1)/e_mag);
        else
            w = 360 - acosd(e_vec(1)/e_mag);
        end
      else % inclined elliptical
        % tan_w = dot(h_hat, cross(n_hat,e_vec))/dot(n_hat,e_vec)
        % so taking the cosine will use the denomenator dot(n_hat,e_vec);
        % We can get n_hat by dividing n / (n_mag) and note
        % e_mag should =1 so it doesnt really matter
        if e_vec(3) >= 0
            w = acosd(dot(n,e_vec) / (n_mag * e_mag));
        else
            w = 360 - acosd(dot(n,e_vec) / (n_mag * e_mag));
        end
      end
    else % circular orbit
       w = [] % angle of perigee undefined
    end
    
    % True anomaly (theta0) and related angles
    if e_mag > 1e-10 % elliptical orbit
        if dot(r, v) >= 0
            f = acosd(dot(e_vec, r)/ (e_mag * norm_r));
        else
            f = 360 - acosd(dot(e_vec, r)/ (e_mag * norm_r));
        end
    elseif i==0 || i==180 % equatorial circular orbit (true longitude at epoch)
        if r(2) >= 0
            f = acosd(r(1)/norm_r);
        else
            f = 360 - acosd(r(1)/norm_r);
        end
    else % circular inclined orbit (argument of latitude u)
        if r(3) >= 0
            f = acosd(dot(n,r)/(n_mag*norm_r));
        else
            f = 360 - acosd(dot(n,r)/(n_mag*norm_r));
        end
    end

    E = 2 * atan2( sqrt(1-e_mag)*sind(f/2), sqrt(1+e_mag)*cosd(f/2) );
    M = mod(rad2deg(E - e_mag*sin(E)), 360);

    %theta0 = f + w;
    
    oe = [a; e_mag; i; omega; w; M];
    
end