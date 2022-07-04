%% Two DOF Manipulator, Lagrangian Formulation

function [tau1, tau2] = TwoDOF() 
    %% Define variables
 
    % Time
    syms t 'real'
    % States
    syms x1 y1 x2 y2 theta1(t) theta2(t) v1 v2 'real'
    % Input
    syms tau1 tau2 'real'
    % Constants
    syms g m1 m2 l1 l2 r1 r2 I1 I2 'real'
    % Energies
    syms K P 'real'
    
    %% Define states
    
    % Position
    x1 = r1 * sin(theta1); 
    y1 = r1 * cos(theta1);
    
    x2 = l1 * sin(theta1) + r2 * sin(theta1 + theta2);
    y2 = l1 * cos(theta1) + r2 * cos(theta1 + theta2);
    
    % Velocity
    x1_dot = diff(x1, t);
    y1_dot = diff(y1, t);
    x2_dot = diff(x2, t);
    y2_dot = diff(y2, t);
    
    v1 = sqrt(x1_dot^2 + y1_dot^2);
    v2 = sqrt((x2_dot^2 + y2_dot^2));
    
    w1 = diff(theta1);
    w2 = diff(theta2);
    
    %% Define Energies
    
    % KE
    K1 = 0.5 * m1 * v1^2 + 0.5 * I1 * w1^2;
    K2 = 0.5 * m2 * v2^2 + 0.5 * I2 * w2^2;
    K =K1 + K2;
    
    % PE
    P1 = m1 * g * y1;
    P2 = m2 * g * y2;
    P = P1 + P2;
    
    %% Lagrangian
    
    L = K - P;
    
    % Tau1
    partial_theta1 = diff(L, theta1);
    partial_theta1_dot = diff(L, diff(theta1, t));
    dT = diff(partial_theta1_dot, t);
    tau1 = dT - partial_theta1;
    tau1 = simplify(tau1);
    
    % Tau2
    partial_theta2 = diff(L, theta2);
    partial_theta2_dot = diff(L, diff(theta2, t));
    dT = diff(partial_theta2_dot, t);
    tau2 = dT - partial_theta2;
    tau2 = simplify(tau2);

end

