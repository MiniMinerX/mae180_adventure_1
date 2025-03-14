function [r0,v0] = gauss(lat, lst, alt, ra, dec, JD, JD_prop)
    % units are in seconds and degrees
    
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

    Re = Rearth + alt; % from the annotation reference
    
    D_conv = [
        sin(lat)*cos(lot), -sin(lot), cos(lat)*cos(lot);
        sin(lat)*sin(lot), cos(lot), cos(lat)*sin(lot);
        -cos(lat), 0, sin(lat);
    ];
    
    Rsite = zeros(3,3);
    for j = 1 : 3
            Rsite(:,j) = Re * [ cos(lat) * cos(lst(j)) ; cos(lat) * sin(lst(j)) ; sin(lat)];
    end
    Rsite = Rsite';
    
    %compute the direction cosine vectors
    
    phat = zeros(3,3);
    
    for i = 1 : 3
        phat(i,:) = [cos(dec(i)) * cos(ra(i)) cos(dec(i)) * sin(ra(i)) sin(dec(i))];
    end
    %phat = phat';
    
    %computation of tau
    tau1 = JD(1) - JD(2);
    tau3 = JD(3) - JD(2);
    tau  = tau3  - tau1;
    
    
    p1 = cross(phat(2,:), phat(3,:));
    p2 = cross(phat(1,:), phat(3,:));
    p3 = cross(phat(1,:), phat(2,:));
    p= [p1 ; p2 ; p3];
    
    %p = D_conv*p_sez% ECI
    
    Do = dot(phat(1,:),p1);
    D  = zeros(3,3);
    for i = 1:3
        for j = 1:3
            D(i,j) = [dot( Rsite(i,:) , p(j,:) )];
        end
    end
    
    A = 1 / Do * (-D(1,2) * tau3 / tau + D(2,2) + D(3,2) * tau1 / tau); % A in here have problem
    B = 1 / 6 / Do * (D(1,2) * (tau3^2 - tau^2) * tau3 / tau + D(3,2) * (tau^2 - tau1^2) * tau1 / tau);
    
    %compute E
    E = dot(Rsite(2,:), phat(2,:));
    R2_2 = dot(Rsite(2,:), Rsite(2,:));
    
    %compute a,b,c
    a1 = -(A^2 + 2*A*E + R2_2);
    b = -2 * mu * B * ( A + E );
    c = -mu^2 * B^2;
    
    %According generate diagram
    r2 = 6000;
    %compute converge R2
    tol = 1e-6; % Convergence tolerance
    max_iter = 100;
    for i = 1:max_iter
        r2_new = r2 - ( (r2^8 + a1 * r2^6 + b * r2^3 + c) / (8 * r2^7 + a1 * 6 * r2^5 + b * 3 * r2^2 ) );
        if abs(r2_new - r2) < tol
            break;
        end
        r2 = r2_new;
    end
    
    %compute the slant range
    slant2 = A + (mu * B) / r2^3;
    
    slant1 = 1 / Do * ( ( (6 * (D(3,1) * tau1 / tau3 + D(2,1) * tau / tau3 ) * r2^3 + mu * D(3,1) * ...
             (tau^2 - tau1^2) * tau1 / tau3 ) / (6 * r2^3 + mu * (tau^2 - tau3^2)) )- D(1,1));
    
    slant3 = 1 / Do * ( ( (6 * (D(1,3) * tau3 / tau1 - D(2,3) * tau / tau1 ) * r2^3 + mu * D(1,3) * ...
             (tau^2 - tau3^2) * tau3 / tau1 ) / (6 * r2^3 + mu * (tau^2 - tau1^2)) )- D(3,3));  
    
    
    %compute the radius position
    %r_vec  = zeros(3,3);
    r1_vec = Rsite(1,:) + slant1 .* phat(1,:);
    r2_vec = Rsite(2,:) + slant2 .* phat(2,:);
    r3_vec = Rsite(3,:) + slant3 .* phat(3,:);
    r_vec  = [r1_vec ; r2_vec ; r3_vec];
    
    %compute the velocity 
    f1 = 1 - 0.5 * mu / r2^3 * tau1^2;
    f3 = 1 - 0.5 * mu / r2^3 * tau3^2;
    g1 = tau1 - 1 / 6 * mu / r2^3 * tau1^3;
    g3 = tau3 - 1 / 6 * mu / r2^3 * tau3^3;
    
    %compute the v2 
    v2 = 1 / (f1 * g3 - f3 * g1) * (-f3 .* r1_vec + f1 .* r3_vec);
    
    t31 = JD(3)-JD(1);
    t32 = JD(3)-JD(2);
    t21 = JD(2)-JD(1);
    
    summer1 = -t32*(1/(t21*t31) + mu/(12*(norm(r1_vec)^3)))*r1_vec;
    summer2 = (t32-t21)*(1/(t21*t32)+mu/(12*(norm(r2_vec)^3)))*r2_vec;
    summer3 = t21*(1/(t32*t31)+mu/(12*(norm(r3_vec)^3)))*r3_vec;
    
    v2_vec = summer1+summer2+summer3; %[output:4d775596]

    r0 = r2_vec;
    v0 = v2_vec;
end