% Code in association with Lauren Stanley

%% Initialization
clf 
% Time Step Settings
t0 = 0;
tf = 65;
t = 0:0.01:tf;
% Initial Conditions
initial_x=0;
initial_y=0;
initial_z=0;
initial_phi=0;
initial_theta=0;
initial_psi=0;
initial_x_dot=0;
initial_y_dot=0;
initial_z_dot=0;
initial_phi_dot=0;
initial_theta_dot=0;
initial_psi_dot=0;

%% Computation
% Reduce ODE Overhead
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

% Call ODE
[t,X]=ode45(@ode_quadrotor,t,[initial_x;initial_y;initial_z;initial_phi;initial_theta;initial_psi;...
    initial_x_dot;initial_y_dot;initial_z_dot;initial_phi_dot;initial_theta_dot;initial_psi_dot],opts);

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
g=9.81;
omega = 0;

% Get desired trajectories for plotting
[x_desiredplot,y_desiredplot,z_desiredplot]=TrajectoryGeneratorv0();

% Compute Inputs 
u1list = [];
u2list = [];
u3list = [];
u4list = [];

w1list = [];
w2list = [];
w3list = [];
w4list = [];

phi_desired_list = [];
theta_desired_list = [];
psi_desired_list = [];
% This for-loop recalculates the inputs are other variables from the ODE
% for use in analysis and plotting
for i=1:1:length(t)
    t_current = t(i);
    [x_desired, y_desired, z_desired, x_dot_desired, y_dot_desired, z_dot_desired, x_double_dot_desired, y_double_dot_desired, z_double_dot_desired] = TrajectoryGenerator(t_current);       
    
    % Pull states and their velocities from ODE results
    x = X(i,1);
    y = X(i,2);
    z = X(i,3);
    phi = wrapToPi(X(i,4));
    theta = wrapToPi(X(i,5));
    psi = wrapToPi(X(i,6));

    x_dot = X(i,7);
    y_dot = X(i,8);
    z_dot = X(i,9);
    theta_dot = X(i,10);
    phi_dot = X(i,11);
    psi_dot = X(i,12);
    
    % Relationship used to calculate desired roll and pitch from x and y
    % desired trajectories
    kp = .01;
    kd = 25;
    Fx = m * (-kp * (x - x_desired) - kd * (x_dot - x_dot_desired) + x_double_dot_desired);
    Fy = m * (-kp * (y - y_desired) - kd * (y_dot - y_dot_desired) + y_double_dot_desired);
    
    % Gains for each stage of the robots motion
    if t_current<=5
        K1 = 0.1;
        K2 = 0.1;
        K3 = 0.1;
        K4 = 0.1;
    elseif 5<t_current && t_current<=20
        K1 = 0.1;
        K2 = 0.001;
        K3 = .001;
        K4 = .1;
     elseif 20<t_current && t_current<=35
        K1 = 0.1;
        K2 = .001;
        K3 = .1;
        K4 = .1;
    elseif 35<t_current && t_current<=50
        K1 = 0.1;
        K2 = 0.1;
        K3 = 0.1;
        K4 = 0.1;
    elseif 50<t_current && t_current<=65
        K1 = 0.1;
        K2 = 0.1;
        K3 = 0.1;
        K4 = 0.1;
    end

    % Saturation Limits
    st_z = 0.5;
    st_phi = 0.5;
    st_theta = 0.5;
    st_psi = 0.5;

    % Calculate u1
    e_z = z - z_desired;
    e_dot_z = z_dot - z_dot_desired;
    lambda = 1;

    s_z = e_dot_z+lambda*e_z;
    sat_z = min(max(s_z, -st_z),st_z);
    %ur1 = -K1*sign(s_z);
    ur1 = -K1*sat_z;
    u1 = ((g + z_double_dot_desired -lambda*e_dot_z)/((cos(phi)*cos(theta))/m))*1  + ur1;
    
    % Calculate desired roll, pitch, and yaw with their associated
    % velocities and accelerations
    phi_desired = asin(-Fy/u1);
    theta_desired = asin(Fx/u1);
    psi_desired = 0;
    
    phi_dot_desired = 0;
    theta_dot_desired = 0;
    psi_dot_desired = 0;
        
    phi_double_dot_desired = 0;
    theta_double_dot_desired = 0;
    psi_double_dot_desired = 0;
    
    % Error Dynamics
    e_phi = wrapToPi(phi - phi_desired);
    e_theta = wrapToPi(theta - theta_desired);
    e_psi = wrapToPi(psi - psi_desired);
    
    e_dot_phi = phi_dot - phi_dot_desired;
    e_dot_theta = theta_dot - theta_dot_desired;
    e_dot_psi = psi_dot - psi_dot_desired;
  
    % Sliding Surfaces for phi, theta and psi
    s_phi = e_dot_phi+lambda*e_phi;
    sat_phi = min(max(s_phi, -st_phi),st_phi);
    s_theta = e_dot_theta+lambda*e_theta;
    sat_theta = min(max(s_theta, -st_theta),st_theta);
    s_psi = e_dot_psi+lambda*e_psi;
    sat_psi = min(max(s_psi, -st_psi),st_psi);

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

    % Append Desired Angles
    phi_desired_list(end+1) = phi_desired;
    theta_desired_list(end+1) = theta_desired;
    psi_desired_list(end+1) = psi_desired;

    % Append inputs to lists for plotting
    u1list(end+1) = u1;
    u2list(end+1) = u2;
    u3list(end+1) = u3;
    u4list(end+1) = u4;
    
    % Append rotor speeds to lists for plotting
    w1list(end+1) = w1;
    w2list(end+1) = w2;
    w3list(end+1) = w3;
    w4list(end+1) = w4;
end

%% Visualization
% Position Plots (x,y,z)
subplot(3,1,1)
fplot(x_desiredplot,[t0,tf],'Color',[0 0 1],'LineWidth',4.0);
hold on
plot(t,X(:,1),'Color',[1 0 0],'LineWidth',1.5);
xlabel('t (s)')
ylabel('m')
title('X Position');

subplot(3,1,2)
fplot(y_desiredplot,[t0,tf],'Color',[0 0 1],'LineWidth',4.0);
hold on
plot(t,X(:,2),'Color',[1 0 0],'LineWidth',1.5);
xlabel('t (s)');
ylabel('m')
title('Y Position');

subplot(3,1,3)
fplot(z_desiredplot,[t0,tf],'Color',[0 0 1],'LineWidth',4.0);
hold on
plot(t,X(:,3),'Color',[1 0 0],'LineWidth',1.5);
xlabel('t (s)');
ylabel('m')
title('Z Position');

pause();
clf;

% Orientation Plots (roll, pitch, yaw)
subplot(3,1,1)
plot(t,phi_desired_list,'Color',[0 0 1],'LineWidth',4.0)
hold on
plot(t,X(:,4),'Color',[1 0 0],'LineWidth',1.5);
axis([0 65 -5*10^-3 5*10^-3])
xlabel('t (s)');
ylabel('rad')
title('Phi');

subplot(3,1,2)
plot(t,theta_desired_list,'Color',[0 0 1],'LineWidth',4.0)
hold on
plot(t,X(:,5),'Color',[1 0 0],'LineWidth',1.5);
axis([0 65 -5*10^-3 5*10^-3])
xlabel('t (s)');
ylabel('rad')
title('Theta');

subplot(3,1,3)
plot(t,psi_desired_list,'Color',[0 0 1],'LineWidth',4.0)
hold on
plot(t,X(:,6),'Color',[1 0 0],'LineWidth',1.5);
axis([0 65 -1 1])
xlabel('t (s)');
ylabel('rad')
title('Psi');

pause();
clf;

% Input Plots
subplot(2,2,1)
plot(t,u1list,'Color',[0    0    0],'LineWidth',2.0);
axis([0 65 .225 .275])
xlabel('t (s)');
ylabel('N')
title('Input 1');

subplot(2,2,2)
plot(t,u2list,'Color',[0    0    0],'LineWidth',2.0);
axis([0 65 -5*10^-4 5*10^-4])
xlabel('t (s)');
ylabel('Nm')
title('Input 2');

subplot(2,2,3)
plot(t,u3list,'Color',[0    0    0],'LineWidth',2.0);
axis([0 65 -5*10^-4 5*10^-4])
xlabel('t (s)');
ylabel('Nm')
title('Input 3');

subplot(2,2,4)
plot(t,u4list,'Color',[0    0    0],'LineWidth',2.0);
axis([0 65 -1 1])
xlabel('t (s)');
ylabel('Nm')
title('Input 4');

pause();
clf;

% Rotor Speed Plots
subplot(2,2,1)
plot(t,w1list,'Color',[0    0    0],'LineWidth',2.0);
xlabel('t (s)');
ylabel('Rad/s')
title('Rotor 1');

subplot(2,2,2)
plot(t,w2list,'Color',[0    0    0],'LineWidth',2.0);
xlabel('t (s)');
ylabel('Rad/s')
title('Rotor 2');

subplot(2,2,3)
plot(t,w3list,'Color',[0    0    0],'LineWidth',2.0);
xlabel('t (s)');
ylabel('Rad/s')
title('Rotor 3');

subplot(2,2,4)
plot(t,w4list,'Color',[0    0    0],'LineWidth',2.0);
xlabel('t (s)');
ylabel('Rad/s')
title('Rotor 4');
