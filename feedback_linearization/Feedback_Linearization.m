% Time
syms t 'real'

% Joint Trajectories as found from calculations
joint1Traj = (180 - 5.4 * t^2 + 18 / 50 * t^3) * pi / 180;
joint2Traj = (90 - 2.7 * t^2 + 9 / 50 * t^3) * pi / 180;

% States
syms theta1 theta1_dot theta1_double_dot theta2 theta2_dot theta2_double_dot 'real'
% Inputs
syms tau1 tau2 'real'
% Constants
syms g m1 m2 l1 l2 r1 r2 I1 I2 'real'

% Equations
tau1_prime = I1*theta1_double_dot + 4*I2*theta1_double_dot + 2*I2*theta2_double_dot + l1^2*m2*theta1_double_dot + m1*r1^2*theta1_double_dot + m2*r2^2*theta1_double_dot + m2*r2^2*theta2_double_dot - g*m1*r1*sin(theta1) - g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - l1*m2*r2*sin(theta2)*theta2_dot^2 + 2*l1*m2*r2*cos(theta2)*theta1_double_dot + l1*m2*r2*cos(theta2)*theta2_double_dot - 2*l1*m2*r2*sin(theta2)*theta1_dot*theta2_dot;
tau2_prime = 2*I2*theta1_double_dot + I2*theta2_double_dot + m2*r2^2*theta1_double_dot + m2*r2^2*theta2_double_dot - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*sin(theta2)*theta1_dot^2 + l1*m2*r2*cos(theta2)*theta1_double_dot;

% Manipulator Form
% Get g
g1 = subs(tau1_prime, [theta1_dot theta1_double_dot theta2_dot theta2_double_dot], [0 0 0 0]);
g2 = subs(tau2_prime, [theta1_dot theta1_double_dot theta2_dot theta2_double_dot], [0 0 0 0]);

% Get M
M1 = subs(tau1 - g1,[theta1_dot theta2_dot], [0 0]);
M2 = subs(tau2 - g2,[theta1_dot theta2_dot], [0 0]);

% Get C
C1 = tau1 - M1 - g1;
C2 = tau2 - M2 - g2;

% ODE
[t,y] = ode45(@myode, [0, 10], [200 * pi / 180, 125 * pi / 180, 0, 0]);

joint1Traj = 180 - 5.4 * t.^2 + 18 / 50 * t.^3;
joint2Traj = 90 - 2.7 * t.^2 + 9 / 50 * t.^3;
theta1_dot_traj = (-10.8 * t + 54 / 50 * t.^2) * pi / 180;
theta2_dot_traj = (-5.4 * t + 27 / 50 * t.^2) * pi / 180;

for i=1:size(y,1)
    theta1 = y(i, 1);
    theta2 = y(i, 2);
    theta1_dot = y(i, 3);
    theta2_dot = y(i, 4);

    % Variables
    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 1; % not needed
    r1 = 0.45;
    r2 = 0.45;
    I1 = 0.084;
    I2 = 0.084;
    g = 9.81;
    
    G = [-g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1); -g*m2*r2*sin(theta1 + theta2)];
    M11 = I1 + 4*I2 + l1^2*m2 + m1*r1^2 + m2*r2^2 + 2*l1*m2*r2*cos(theta2);
    M12 = 2*I2 + m2*r2^2 + l1*m2*r2*cos(theta2);
    M21 = 2*I2 + m2*r2^2 + l1*m2*r2*cos(theta2);
    M22 = I2 + m2*r2^2;
    M = [M11, M12; M21, M22];
    C = [0 -l1*m2*r2*sin(theta2)*theta2_dot - 2*l1*m2*r2*theta1_dot*sin(theta2); l1*m2*r2*theta1_dot*sin(theta2) 0];
    K = [2 0 3 0; 0 6 0 5];

    % Desired Trajectories
    theta1_desired = (180 - 5.4 * t(i)^2 + 18 / 50 * t(i)^3) * pi / 180;
    theta1_dot_desired = (-10.8 * t(i) + 54 / 50 * t(i)^2) * pi / 180;
    theta1_double_dot_desired = (-10.8 + 10.8 / 5 * t(i)) * pi / 180;

    theta2_desired = (90 - 2.7 * t(i)^2 + 9 / 50 * t(i)^3) * pi / 180;
    theta2_dot_desired = (-5.4 * t(i) + 27 / 50 * t(i)^2) * pi / 180;
    theta2_double_dot_desired = (-5.4 + 27 / 25 * t(i)) * pi / 180;
    
    % Feedback Linearization
    tau = M * (-K * ([theta1; theta2; theta1_dot; theta2_dot] - [theta1_desired; theta2_desired; theta1_dot_desired; theta2_dot_desired]) + [theta1_double_dot_desired; theta2_double_dot_desired]) + C * [theta1_dot; theta2_dot] + G;
    tau1(i) = tau(1, :);
    tau2(i) = tau(2, :);

end

% plot(t, (joint1Traj))
% hold on
% plot(t, (y(:, 1) * 180 / pi))
% xlabel('Time')
% ylabel('theta1')
% legend('desired','output')
% figure;
% plot(t, (joint2Traj))
% hold on
% plot(t, (y(:, 2) * 180 / pi))
% xlabel('Time')
% ylabel('theta2')
% legend('desired','output')
% 
% figure;
% plot(t, (theta1_dot_traj))
% hold on
% plot(t, (y(:, 3) ))
% xlabel('Time')
% ylabel('theta1dot')
% legend('desired','output')
% figure;
% plot(t, (theta2_dot_traj))
% hold on
% plot(t, (y(:, 4)) )
% xlabel('Time')
% ylabel('theta2dot')
% legend('desired','output')

figure;
plot(t, (tau1))
xlabel('Time')
ylabel('Tau1')
figure;
plot(t, (tau2))
xlabel('Time')
ylabel('Tau2')


% Plotting
% figure;
% plot(t, y(:, 1));
% xlabel('Time')
% ylabel('theta1')
% figure;
% plot(t, y(:, 2));
% xlabel('Time')
% ylabel('theta2')
% figure;
% plot(t, y(:, 3));
% xlabel('Time')
% ylabel('theta1d')
% figure;
% plot(t, y(:, 4));
% xlabel('Time')
% ylabel('theta2d')
% figure;
% plot(t, tau1);
% xlabel('Time')
% ylabel('Torque1')
% figure;
% plot(t, tau2);
% xlabel('Time')
% ylabel('Torque2')

%% ODE Function
function dX = myode(t, X)

    % Variables
    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 1; % not needed
    r1 = 0.45;
    r2 = 0.45;
    I1 = 0.084;
    I2 = 0.084;
    g = 9.81;

    dX = zeros(4, 1);
    X = num2cell(X);
    [theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});

    % I get an error with wrapTo2pi
    if abs(theta1) > 2 * pi
        theta1 = mod(theta1, 2 * pi);
    end

    if abs(theta2) > 2 * pi
        theta2 = mod(theta2, 2 * pi);
    end

    G = [-g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1); -g*m2*r2*sin(theta1 + theta2)];
    M11 = I1 + 4*I2 + l1^2*m2 + m1*r1^2 + m2*r2^2 + 2*l1*m2*r2*cos(theta2);
    M12 = 2*I2 + m2*r2^2 + l1*m2*r2*cos(theta2);
    M21 = 2*I2 + m2*r2^2 + l1*m2*r2*cos(theta2);
    M22 = I2 + m2*r2^2;
    M = [M11, M12; M21, M22];
    C = [0 -l1*m2*r2*sin(theta2)*theta2_dot - 2*l1*m2*r2*theta1_dot*sin(theta2); l1*m2*r2*theta1_dot*sin(theta2) 0];
    K = [2 0 3 0; 0 6 0 5];

    % Desired Trajectories
    theta1_desired = (180 - 5.4 * t^2 + 18 / 50 * t^3) * pi / 180;
    theta1_dot_desired = (-10.8 * t + 54 / 50 * t^2) * pi / 180;
    theta1_double_dot_desired = (-10.8 + 10.8 / 5 * t) * pi / 180;

    theta2_desired = (90 - 2.7 * t^2 + 9 / 50 * t^3) * pi / 180;
    theta2_dot_desired = (-5.4 * t + 27 / 50 * t^2) * pi / 180;
    theta2_double_dot_desired = (-5.4 + 27 / 25 * t) * pi / 180;
    
    % Feedback Linearization
    tau = M * (-K * ([theta1; theta2; theta1_dot; theta2_dot] - [theta1_desired; theta2_desired; theta1_dot_desired; theta2_dot_desired]) + [theta1_double_dot_desired; theta2_double_dot_desired]) + C * [theta1_dot; theta2_dot] + G;
    tau1 = tau(1, :);
    tau2 = tau(2, :);

    % State vector
    dX(1) = theta1_dot;
    dX(2) = theta2_dot;
    dX(3) = -(382000*tau1 - 382000*tau2 + 5433759*sin(theta1) + 171900*theta1_dot^2*sin(theta2) + 171900*theta2_dot^2*sin(theta2) - 600000*tau2*cos(theta2) - 2648700*sin(theta1 + theta2)*cos(theta2) + 343800*theta1_dot*theta2_dot*sin(theta2) + 270000*theta1_dot^2*cos(theta2)*sin(theta2))/(270000*cos(theta2)^2 - 491443);
    dX(4) = (1146000*tau1 - 6292000*tau2 - 22717017*sin(theta1 + theta2) + 16301277*sin(theta1) + 2831400*theta1_dot^2*sin(theta2) + 515700*theta2_dot^2*sin(theta2) + 25604100*cos(theta2)*sin(theta1) + 1800000*tau1*cos(theta2) - 3600000*tau2*cos(theta2) - 7946100*sin(theta1 + theta2)*cos(theta2) + 1031400*theta1_dot*theta2_dot*sin(theta2) + 1620000*theta1_dot^2*cos(theta2)*sin(theta2) + 810000*theta2_dot^2*cos(theta2)*sin(theta2) + 1620000*theta1_dot*theta2_dot*cos(theta2)*sin(theta2))/(3*(270000*cos(theta2)^2 - 491443));

end