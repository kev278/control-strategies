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

while(toc < 10)
    
    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    % design your state feedback controller in the following
    
    X = [jointData.Position(1);jointData.Position(2);jointData.Velocity(1); jointData.Velocity(2)];
    K = [23.5850, 5.8875, 5.1470, 2.6104;
         5.8875,  4.9875, 1.5543, 0.9970];

    U = -K * X;
   
   tau1.Data = U(1);
   tau2.Data = U(2);
   send(j1_effort,tau1);
   send(j2_effort,tau2);
   
   X1(sample) = jointData.Position(1);
   X2(sample) = jointData.Position(2);
   X3(sample) = jointData.Velocity(1);
   X4(sample) = jointData.Velocity(2);

   Tau1(sample) = tau1.Data;
   Tau2(sample) = tau2.Data;  

   t(sample) = toc;

   sample = sample + 1;
   
% you can sample data here to be plotted at the end
end

figure;
plot(t, X1);
xlabel('Time')
ylabel('theta1')
figure;
plot(t, X2);
xlabel('Time')
ylabel('theta2')
figure;
plot(t, X3);
xlabel('Time')
ylabel('theta1d')
figure;
plot(t, X4);
xlabel('Time')
ylabel('theta2d')
figure;
plot(t, Tau1);
xlabel('Time')
ylabel('Torque1')
figure;
plot(t, Tau2);
xlabel('Time')
ylabel('Torque2')


tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;