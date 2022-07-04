clear all
% Time
syms t 'real'
% States
syms theta1 theta1_dot theta1_double_dot theta2 theta2_dot theta2_double_dot 'real'
% Inputs
syms tau1 tau2 'real'
% Constants
syms g m1 m2 l1 l2 r1 r2 I1 I2 'real'

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

% Desired Trajectories
theta1_desired = (180 - 5.4 * t^2 + 18 / 50 * t^3) * pi / 180;
theta1_dot_desired = (-10.8 * t + 54 / 50 * t^2) * pi / 180;
theta1_double_dot_desired = (-10.8 + 10.8 / 5 * t) * pi / 180;

theta2_desired = (90 - 2.7 * t^2 + 9 / 50 * t^3) * pi / 180;
theta2_dot_desired = (-5.4 * t + 27 / 50 * t^2) * pi / 180;
theta2_double_dot_desired = (-5.4 + 27 / 25 * t) * pi / 180;

% Dynamics
a = I1 + I2 + m1 * r1^2 + m2 * (l1^2 + r2^2);
b = m2 * l1 * r2;
d = I2 + m2 * r2^2;

M = [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
C = [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
G = [-m1*g*r1*sin(theta1)-m2*g*(l1*sin(theta1)+r2*sin(theta1+theta2)); -m2*g*r2*sin(theta1+theta2)];

% State Feedback for Gains
eigen_values = [-3, -3, -4, -4];
A = [0, 0, 1, 0;
     0, 0, 0, 1; 
     0, 0, 0, 0; 
     0, 0, 0, 0];

B = [0, 0;
     0, 0;
     1, 0;
     0, 1];

K = place(A, B, eigen_values);
Kp = K(1:2, 1:2);
Kd = K(1:2, 3:4);

%Find Matrix P
Q = eye(4) * 10;
Acl = A - B * K;
P = lyap(Acl', Q);
% P = 0;

% Alpha
alpha = [m2*l1^2 + m1*r1^2 + m2*r2^2 + I1 + I2
        m2*l1*r2
        m2*r2^2 + I2
        m1*r1 + m2*l1
        m2*r2];
alpha = 0.75 * alpha;

% ODE Output
[t, y] = ode45(@(t,y) myode(t, y, K, B, P), [0, 10], [200 * pi / 180; 125 * pi / 180; 0; 0; alpha]);

% Adaptive Control
for i = 1:size(y,1)

    theta1 = y(i, 1);
    theta2 = y(i, 2);
    theta1_dot = y(i, 3);
    theta2_dot = y(i, 4);

    a = y(i, 5);
    b = y(i, 6);
    c = y(i, 7);
    g1 = y(i, 8);
    g2 = y(i, 9);

    % Estimates
    M_hat = [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
    C_hat = [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
    G_hat = [-g1*g*sin(theta1)-g2*g*sin((theta1+theta2)); -g2*g*sin(theta1+theta2)];

    % Desired Trajectories
    theta1_desired = (180 - 5.4 * t(i)^2 + 18 / 50 * t(i)^3) * pi / 180;
    theta1_dot_desired = (-10.8 * t(i) + 54 / 50 * t(i)^2) * pi / 180;
    theta1_double_dot_desired = (-10.8 + 10.8 / 5 * t(i)) * pi / 180;

    theta2_desired = (90 - 2.7 * t(i)^2 + 9 / 50 * t(i)^3) * pi / 180;
    theta2_dot_desired = (-5.4 * t(i) + 27 / 50 * t(i)^2) * pi / 180;
    theta2_double_dot_desired = (-5.4 + 27 / 25 * t(i)) * pi / 180;
    
    % Input to cancel uncertainities
    x = [theta1 - theta1_desired; theta2 - theta2_desired; theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];

    % Robust Control
    v = [theta1_double_dot_desired; theta2_double_dot_desired] - K * x;
    tau = M_hat * v + C_hat * [theta1_dot; theta2_dot] + G_hat;
    tau1(i) = tau(1, :);
    tau2(i) = tau(2, :);

end

% Desired Trajectories
theta1_desired = 180 - 5.4 * t.^2 + 18 / 50 * t.^3;
theta2_desired = 90 - 2.7 * t.^2 + 9 / 50 * t.^3;
theta1_dot_desired = (-10.8 * t + 54 / 50 * t.^2) * pi / 180;
theta2_dot_desired = (-5.4 * t + 27 / 50 * t.^2) * pi / 180;

% Plotting
plot(t, (theta1_desired))
hold on
plot(t, (y(:, 1) * 180 / pi))
xlabel('Time')
ylabel('theta1')
legend('desired','output')
figure;
plot(t, (theta2_desired))
hold on
plot(t, (y(:, 2) * 180 / pi))
xlabel('Time')
ylabel('theta2')
legend('desired','output')

figure;
plot(t, (theta1_dot_desired))
hold on
plot(t, (y(:, 3) ))
xlabel('Time')
ylabel('theta1dot')
legend('desired','output')
figure;
plot(t, (theta2_dot_desired))
hold on
plot(t, (y(:, 4)) )
xlabel('Time')
ylabel('theta2dot')
legend('desired','output')

figure;
plot(t, y(:, 5))
hold on
plot(t, y(:, 6))
plot(t, y(:, 7))
plot(t, y(:, 8))
plot(t, y(:, 9))
xlabel('Time')
ylabel('Alpha')
legend('alpha1', 'alpha2', 'alpha3', 'alpha4', 'alpha5');

figure;
plot(t, (tau1))
xlabel('Time')
ylabel('Tau1')
figure;
plot(t, (tau2))
xlabel('Time')
ylabel('Tau2')




%% ODE Function
function dX = myode(t, X, K, B, P)

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

    dX = zeros(9, 1);
    X = num2cell(X);
    [theta1, theta2, theta1_dot, theta2_dot, a, b, d, g1, g2] = deal(X{:});

    % I get an error with wrapTo2pi
    if abs(theta1) > 2 * pi
        theta1 = mod(theta1, 2 * pi);
    end

    if abs(theta2) > 2 * pi
        theta2 = mod(theta2, 2 * pi);
    end

    % Estimates
    M_hat = [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
    C_hat = [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
    G_hat = [-g1*g*sin(theta1)-g2*g*sin((theta1+theta2)); -g2*g*sin(theta1+theta2)];

    % Desired Trajectories
    theta1_desired = (180 - 5.4 * t^2 + 18 / 50 * t^3) * pi / 180;
    theta1_dot_desired = (-10.8 * t + 54 / 50 * t^2) * pi / 180;
    theta1_double_dot_desired = (-10.8 + 10.8 / 5 * t) * pi / 180;

    theta2_desired = (90 - 2.7 * t^2 + 9 / 50 * t^3) * pi / 180;
    theta2_dot_desired = (-5.4 * t + 27 / 50 * t^2) * pi / 180;
    theta2_double_dot_desired = (-5.4 + 27 / 25 * t) * pi / 180;

    % Input to cancel uncertainities
    x = [theta1 - theta1_desired; theta2 - theta2_desired; theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];

    % Robust Control
    v = [theta1_double_dot_desired; theta2_double_dot_desired] - K * x;
    tau = M_hat * v + C_hat * [theta1_dot; theta2_dot] + G_hat;
    tau1 = tau(1, :);
    tau2 = tau(2, :);

    % State vector
    dX(1) = theta1_dot;
    dX(2) = theta2_dot;
    dX(3) = -(382000*tau1 - 382000*tau2 + 5433759*sin(theta1) + 171900*theta1_dot^2*sin(theta2) + 171900*theta2_dot^2*sin(theta2) - 600000*tau2*cos(theta2) - 2648700*sin(theta1 + theta2)*cos(theta2) + 343800*theta1_dot*theta2_dot*sin(theta2) + 270000*theta1_dot^2*cos(theta2)*sin(theta2))/(270000*cos(theta2)^2 - 491443);
    dX(4) = (1146000*tau1 - 6292000*tau2 - 22717017*sin(theta1 + theta2) + 16301277*sin(theta1) + 2831400*theta1_dot^2*sin(theta2) + 515700*theta2_dot^2*sin(theta2) + 25604100*cos(theta2)*sin(theta1) + 1800000*tau1*cos(theta2) - 3600000*tau2*cos(theta2) - 7946100*sin(theta1 + theta2)*cos(theta2) + 1031400*theta1_dot*theta2_dot*sin(theta2) + 1620000*theta1_dot^2*cos(theta2)*sin(theta2) + 810000*theta2_dot^2*cos(theta2)*sin(theta2) + 1620000*theta1_dot*theta2_dot*cos(theta2)*sin(theta2))/(3*(270000*cos(theta2)^2 - 491443));

    % Parametric Equations
    q1d = dX(1);
    q2d = dX(2);
    ddq1 = dX(3);
    ddq2 = dX(4);
    q1 = theta1;
    q2 = theta2;

    Y = [ddq1, ...
        cos(q2)*(2*ddq1 + ddq2) - 2*sin(q2)*q1d*q2d - sin(q2)*q2d^2, ...
        ddq2, ...
        -sin(q1)*g, ...
        -sin(q1 + q2)*g; ...
        0, ...
        sin(q2)*q1d^2 + cos(q2)*ddq1, ...
        ddq1 + ddq2, ...
        0, ...
        -sin(q1+q2)*g];

    % Adaptation Law
    phi = M_hat \ Y;
    gamma = eye(5);
    alpha_hat = -gamma \ (phi' * B' * P * x);
    
    dX(5) = alpha_hat(1);
    dX(6) = alpha_hat(2);
    dX(7) = alpha_hat(3);
    dX(8) = alpha_hat(4);
    dX(9) = alpha_hat(5);

end



