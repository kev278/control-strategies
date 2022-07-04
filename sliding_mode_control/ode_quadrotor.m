% In association with Lauren Stanley

%ODE Function for quadrotor
function dxdt=ode_quadrotor(t,X)
    % Print Time
    %t
    % Variables in SI Units
    m = 27 * 10^-3;
    l = 46 * 10^-3;
    Ix = 16.571710 * 10^-6;
    Iy = 16.571710 * 10^-6;
    Iz = 29.261652 * 10^-6;
    Ip = 12.65625 * 10^-8;
    kF = 1.28192 * 10^-8;
    kM = 5.964552 * 10^-3;
    wMax = 2618;  % in rad/s
    wMin = 0;     % in rad/s
    g = 9.81;

    % Set omega as a global variable so it can be used between interations
    global omega
    if t == 0
        omega = 0;
    end

    % States
    x = X(1);
    y = X(2);
    z = X(3);
    phi = X(4);
    theta = X(5);
    psi = X(6);

    x_dot = X(7);
    y_dot = X(8);
    z_dot = X(9);
    phi_dot = X(10);
    theta_dot = X(11);
    psi_dot = X(12);
    
    
    % Desired Trajectory values based on time
    [x_desired, y_desired, z_desired, x_dot_desired, y_dot_desired, z_dot_desired, x_double_dot_desired, y_double_dot_desired, z_double_dot_desired] = TrajectoryGenerator(t);
    
    % Relationship used to calculate desired roll and pitch from x and y
    % desired trajectories   
    kp = .01;
    kd = 25;
    Fx = m * (-kp * (x - x_desired) - kd * (x_dot - x_dot_desired) + x_double_dot_desired);
    Fy = m * (-kp * (y - y_desired) - kd * (y_dot - y_dot_desired) + y_double_dot_desired);
    
    if t<=5
        K1 = 0.1;
        K2 = 0.1;
        K3 = 0.1;
        K4 = 0.1;
    elseif 5<t && t<=20
        K1 = 0.1;
        K2 = 0.001;
        K3 = 0.001;
        K4 = .1;
     elseif 20<t && t<=35
        K1 = 0.1;
        K2 = .001;
        K3 = .1;
        K4 = .1;
    elseif 35<t && t<=50
        K1 = 0.1;
        K2 = 0.1;
        K3 = 0.1;
        K4 = 0.1;
    elseif 50<t && t<=65
        K1 = 0.1;
        K2 = 0.1;
        K3 = 0.1;
        K4 = 0.1;
    end
    
    % Saturation Limits
    p=1;
    st_z = 0.5;
    st_phi = 0.5;
    st_theta = 0.5;
    st_psi = 0.5;
    
    % Calculate u1
    e_z = z - z_desired;
    e_dot_z = z_dot - z_dot_desired;
    lambda = 1;

    s_z = e_dot_z+lambda*e_z;
    sat_z = min(max(s_z/p, -st_z),st_z);
    %ur1 = -K1*sign(s_z);
    ur1 = -K1*sat_z;    
    u1 = (((g + z_double_dot_desired -lambda*e_dot_z)/((cos(phi)*cos(theta))/m))  + ur1);
    
    % Calculate desired roll, pitch, and yaw with their associated 
    % velocities and accelerations
    phi_desired = asin(-Fy/u1);
    theta_desired = asin(Fx/u1);
    psi_desired = 0;
    
    % Assumeing zero velocity for roll pitch and yaw
    phi_dot_desired = 0;
    theta_dot_desired = 0;
    psi_dot_desired = 0;
    
    theta_double_dot_desired = 0;
    phi_double_dot_desired = 0;
    psi_double_dot_desired = 0;
    
    % Error Dynamics
    e_theta = wrapToPi(theta - theta_desired);
    e_phi = wrapToPi(phi - phi_desired);
    e_psi = wrapToPi(psi - psi_desired);
    
    e_dot_theta = theta_dot - theta_dot_desired;
    e_dot_phi = phi_dot - phi_dot_desired;
    e_dot_psi = psi_dot - psi_dot_desired;
  
    % Sliding Surface for phi, theta and psi
    s_phi = e_dot_phi+lambda*e_phi;
    sat_phi = min(max(s_phi/p, -st_phi),st_phi);
    s_theta = e_dot_theta+lambda*e_theta;
    sat_theta = min(max(s_theta/p, -st_theta),st_theta);
    s_psi = e_dot_psi+lambda*e_psi;
    sat_psi = min(max(s_psi/p, -st_psi),st_psi);
    
    % Calculate inputs 2, 3, and 4
    %ur2 = -K2*sign(s_phi);
    %ur3 = -K3*sign(s_theta);
    %ur4 = -K4*sign(s_psi);
    
    ur2 = -K2*sat_phi;
    ur3 = -K3*sat_theta;
    ur4 = -K4*sat_psi;
    
    u2 = (((-psi_dot*theta_dot*(Iy - Iz)+Ip*omega*theta_dot)/Ix + phi_double_dot_desired-lambda*e_dot_phi)/(1/Ix)) + ur2;
    u3 = (((-phi_dot*psi_dot*(Iz - Ix)-Ip*omega*phi_dot)/Iy + theta_double_dot_desired-lambda*e_dot_theta)/(1/Iy)) + ur3;
    u4 = (((-phi_dot*theta_dot*(Ix - Iy))/Iz + psi_double_dot_desired - lambda*e_dot_psi)/(1/Iz)) + ur4;
    
    % Use allocation matrix to determine rotor speeds
    AM = [1 / (4 * kF), -sqrt(2) / (4 * kF * l), -sqrt(2) / (4 * kF * l), -1 / (4 * kM * kF);
      1 / (4 * kF), -sqrt(2) / (4 * kF * l),  sqrt(2) / (4 * kF * l),  1 / (4 * kM * kF);
      1 / (4 * kF),  sqrt(2) / (4 * kF * l),  sqrt(2) / (4 * kF * l), -1 / (4 * kM * kF);
      1 / (4 * kF),  sqrt(2) / (4 * kF * l), -sqrt(2) / (4 * kF * l),  1 / (4 * kM * kF)]; 
    u = [u1; u2; u3; u4];
    w = AM * u;
    
    % Enforce rotor speed limitations
    w1 = real(sqrt(w(1)));
    if w1 < wMin
        w1 = wMin;
    elseif w1 > wMax
        w1 = wMax;
    end
    w2 = real(sqrt(w(2)));
    if w2 < wMin
        w2 = wMin;
    elseif w2 > wMax
        w2 = wMax;
    end
    w3 = real(sqrt(w(3)));
    if w1 < wMin
        w3 = wMin;
    elseif w3 > wMax
        w3 = wMax;
    end
    w4 = real(sqrt(w(4)));
    if w4 < wMin
        w4 = wMin;
    elseif w4 > wMax
        w4 = wMax;
    end
    omega = w1 - w2 + w3 - w4;

    % Robot Dynamics
    x_double_dot = (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi)) * u1 / m;
    y_double_dot = (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi)) * u1 / m;
    z_double_dot = (cos(phi) * cos(theta)) * u1 / m - g;
    phi_double_dot = theta_dot * psi_dot * (Iy - Iz) / Ix - Ip * omega * theta_dot / Ix + u2 / Ix;
    theta_double_dot = phi_dot * psi_dot * (Iz - Ix) / Iy + Ip * omega * phi_dot / Iy + u3 / Iy;
    psi_double_dot = phi_dot * theta_dot * (Ix - Iy) / Iz + u4 / Iz;
    
    % Pass variables for next ODE iteration
    dxdt_1 = x_dot;
    dxdt_2 = y_dot;
    dxdt_3 = z_dot;
    dxdt_4 = phi_dot;
    dxdt_5 = theta_dot;  
    dxdt_6 = psi_dot;
    
    dxdt_7 = x_double_dot;
    dxdt_8 = y_double_dot;
    dxdt_9 = z_double_dot;
    dxdt_10 = phi_double_dot;
    dxdt_11 = theta_double_dot;
    dxdt_12 = psi_double_dot;
    
    dxdt=[dxdt_1;dxdt_2;dxdt_3;dxdt_4;dxdt_5;dxdt_6;...
        dxdt_7;dxdt_8;dxdt_9;dxdt_10;dxdt_11;dxdt_12];
end
