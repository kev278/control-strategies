clear; close; clc;
% ROS Setup
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);
tic;
sample = 1;

%Variables
m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1; % not needed
r1 = 0.45;
r2 = 0.45;
I1 = 0.084;
I2 = 0.084;
g = 9.81;

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

% Alpha
alpha = [m2*l1^2 + m1*r1^2 + m2*r2^2 + I1 + I2
        m2*l1*r2
        m2*r2^2 + I2
        m1*r1 + m2*l1
        m2*r2];
alpha = 0.75 * alpha;

a = alpha(1);
b = alpha(2);
d = alpha(3);
g1 = alpha(4);
g2 = alpha(5);

while(toc < 10)

   t = toc;

    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    % design your state feedback controller in the following
    
    X = [wrapTo2Pi(jointData.Position(1)); wrapTo2Pi(jointData.Position(2)) ;jointData.Velocity(1); jointData.Velocity(2); a; b; d; g1; g2];
    
    theta1 = jointData.Position(1);
    theta2 = jointData.Position(2);
    theta1_dot = jointData.Velocity(1); 
    theta2_dot = jointData.Velocity(2);
   
    % Desired Traj
    theta1_desired = (180 - 5.4 * t^2 + 18 / 50 * t^3) * pi / 180;
    theta1_dot_desired = (-10.8 * t + 54 / 50 * t^2) * pi / 180;
    theta1_double_dot_desired = (-10.8 + 10.8 / 5 * t) * pi / 180;

    theta2_desired = (90 - 2.7 * t^2 + 9 / 50 * t^3) * pi / 180;
    theta2_dot_desired = (-5.4 * t + 27 / 50 * t^2) * pi / 180;
    theta2_double_dot_desired = (-5.4 + 27 / 25 * t) * pi / 180;

    % Estimates
    M_hat = [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
    C_hat = [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
    G_hat = [-g1*g*sin(theta1)-g2*g*sin((theta1+theta2)); -g2*g*sin(theta1+theta2)];

    % Input to cancel uncertainities
    x = [theta1 - theta1_desired; theta2 - theta2_desired; theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];

    % Robust Control
    v = [theta1_double_dot_desired; theta2_double_dot_desired] - K * x;
    tau = M_hat * v + C_hat * [theta1_dot; theta2_dot] + G_hat;

   tau1.Data = tau(1, :);
   tau2.Data = tau(2, :);

   send(j1_effort,tau1);
   send(j2_effort,tau2);

   
   X1(sample) = jointData.Position(1);
   X2(sample) = jointData.Position(2);
   X3(sample) = jointData.Velocity(1);
   X4(sample) = jointData.Velocity(2);

   % Parametric Equations
    q1d = X1(sample);
    q2d = X2(sample);
    ddq1 = X3(sample);
    ddq2 = X4(sample);
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
    
    X5(sample) = alpha_hat(1);
    X6(sample) = alpha_hat(2);
    X7(sample) = alpha_hat(3);
    X8(sample) = alpha_hat(4);
    X9(sample) = alpha_hat(5);

   Tau1(sample) = tau1.Data;
   Tau2(sample) = tau2.Data;  

   time(sample) = toc;

   sample = sample + 1;
   
% you can sample data here to be plotted at the end
end

figure;
plot(time, X1);
xlabel('Time')
ylabel('theta1')
figure;
plot(time, X2);
xlabel('Time')
ylabel('theta2')
figure;
plot(time, X3);
xlabel('Time')
ylabel('theta1dot')
figure;
plot(time, X4);
xlabel('Time')
ylabel('theta2dot')
figure;
plot(time, Tau1);
xlabel('Time')
ylabel('Torque1')
figure;
plot(time, Tau2);
xlabel('Time')
ylabel('Torque2')

figure;
plot(time, X5)
hold on
plot(time, X6)
plot(time, X7)
plot(time, X8)
plot(time, X9)
xlabel('Time')
ylabel('Alpha')
legend('alpha1', 'alpha2', 'alpha3', 'alpha4', 'alpha5');


tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;