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

% Nominal values
m1_nominal = 0.75;
m2_nominal = 0.75;
I1_nominal = 0.063;
I2_nominal = 0.063;

I1 = I1_nominal;
I2 = I2_nominal;
m1 = m1_nominal;
m2 = m2_nominal;

% Dynamics
a = I1 + I2 + m1 * r1^2 + m2 * (l1^2 + r2^2);
b = m2 * l1 * r2;
d = I2 + m2 * r2^2;

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

% Upper Bound
rho = 8;
phi = 0.1;

while(toc < 10)

   t = toc;

    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    % design your state feedback controller in the following
    
    X = [wrapTo2Pi(jointData.Position(1)); wrapTo2Pi(jointData.Position(2)) ;jointData.Velocity(1); jointData.Velocity(2)];
    
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

    % Dynamics
    M = [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
    C = [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
    G = [-m1*g*r1*sin(theta1)-m2*g*(l1*sin(theta1)+r2*sin(theta1+theta2)); -m2*g*r2*sin(theta1+theta2)];
    
   % Input to cancel uncertainities
    x = [theta1 - theta1_desired; theta2 - theta2_desired; theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];

    % Boundary Layer
    
    if phi > 0
        if norm(x' * P * B) > phi
            vr = -((x' * P * B) / norm(x' * P * B)) * rho;
        else
            vr = -(x' * P * B) / phi * rho;
        end
    else
        if norm(x' * P * B) == 0
            vr = 0; 
        else
            vr = -((x' * P * B) / norm(x' * P * B)) * rho;
        end
    end

    % Robust Control
    v = [theta1_double_dot_desired; theta2_double_dot_desired] - K * x + vr';
    tau = M * v + C * [theta1_dot; theta2_dot] + G;

   tau1.Data = tau(1, :);
   tau2.Data = tau(2, :);

   send(j1_effort,tau1);
   send(j2_effort,tau2);

   
   X1(sample) = jointData.Position(1);
   X2(sample) = jointData.Position(2);
   X3(sample) = jointData.Velocity(1);
   X4(sample) = jointData.Velocity(2);

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


tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;