%clear all
% Time
syms t 'real'
% States
syms theta1 theta2 theta1_dot theta2_dot theta1_double_dot theta2_double_dot 'real'
% Input
syms tau1 tau2 'real'
% Constants
syms g m1 m2 l1 l2 r1 r2 I1 I2 'real'

%% a)
% tau1 = I1*theta1_double_dot + I2*theta1_double_dot + I2*theta2_double_dot + l1^2*m2*theta1_double_dot + m1*r1^2*theta1_double_dot + m2*r2^2*theta1_double_dot + m2*r2^2*theta2_double_dot - g*m1*r1*sin(theta1) - g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - l1*m2*r2*sin(theta2)*theta2_dot^2 + 2*l1*m2*r2*cos(theta2)*theta1_double_dot + l1*m2*r2*cos(theta2)*theta2_double_dot - 2*l1*m2*r2*sin(theta2)*theta1_dot*theta2_dot;
% tau2 = I2*theta1_double_dot + I2*theta2_double_dot + m2*r2^2*theta1_double_dot + m2*r2^2*theta2_double_dot - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*sin(theta2)*theta1_dot^2 + l1*m2*r2*cos(theta2)*theta1_double_dot;
% 
% EOM = [tau1; tau2];
% EOMEquilibrium = subs(EOM, [theta1_double_dot, theta2_double_dot, theta1_dot, theta2_dot, tau1, tau2], [0, 0, 0, 0, 0, 0]);
% 
% solution = solve(EOMEquilibrium == 0, [theta1, theta2]);
% display(solution.theta1);
% display(solution.theta2);

%% b)

% Comment a) to run b) and further steps

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

X = sym('X', [4, 1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_dot;
X(4) = theta2_dot;

eq1 = -tau1 + I1*theta1_double_dot + I2*theta1_double_dot + I2*theta2_double_dot + l1^2*m2*theta1_double_dot + m1*r1^2*theta1_double_dot + m2*r2^2*theta1_double_dot + m2*r2^2*theta2_double_dot - g*m1*r1*sin(theta1) - g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - l1*m2*r2*sin(theta2)*theta2_dot^2 + 2*l1*m2*r2*cos(theta2)*theta1_double_dot + l1*m2*r2*cos(theta2)*theta2_double_dot - 2*l1*m2*r2*sin(theta2)*theta1_dot*theta2_dot;
eq2 = -tau2 + I2*theta1_double_dot + I2*theta2_double_dot + m2*r2^2*theta1_double_dot + m2*r2^2*theta2_double_dot - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*sin(theta2)*theta1_dot^2 + l1*m2*r2*cos(theta2)*theta1_double_dot;

sol = solve([eq1 == 0, eq2 == 0], [theta1_double_dot, theta2_double_dot]);

dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = sol.theta1_double_dot;
dX(4) = sol.theta2_double_dot;

X_dot = [dX(1); dX(2); dX(3); dX(4)];

A = jacobian(X_dot, X);
B = jacobian(X_dot, [tau1; tau2]);
A_lin = A;
B_lin = B;

%% c)

% About (0, 0)
A1 = A;
B1 = B;
A1 = subs(A1,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);
B1 = subs(B1,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);

A1 = double(A1);
B1 = double(B1);

eigenValues1 = eig(A1);

% About (pi, 0)
A2 = A;
B2 = B;
A2 = subs(A2,[theta1,theta2,theta1_dot,theta2_dot],[pi,0,0,0]);
B2 = subs(B2,[theta1,theta2,theta1_dot,theta2_dot],[pi,0,0,0]);

A2 = double(A2);
B2 = double(B2);

eigenValues2 = eig(A2);

% About (0, pi)
A3 = A;
B3 = B;
A3 = subs(A3,[theta1,theta2,theta1_dot,theta2_dot],[0,pi,0,0]);
B3 = subs(B3,[theta1,theta2,theta1_dot,theta2_dot],[0,pi,0,0]);

A3 = double(A3);
B3 = double(B3);

eigenValues3 = eig(A3);

%% d)
% MATLAB documentation syntax used
Co = ctrb(A1,B1);
unco = length(A1) - rank(Co);
% unco is 0 so system in Controllble

%% e)
syms K11 K12 K13 K14 K21 K22 K23 K24 'real'
syms lambda

K = [K11, K12, K13, K14;
     K21, K22, K23, K24];

Acl = A1 - (B1 * K);
% det of A - BK
characteristicPoly = det(Acl - lambda * eye(4));
characteristicPoly = simplify(characteristicPoly);

lambda = [-0.5, -1, -1.5, -2];
K = place(A1,B1,lambda);

%% f)

% ODE
[t, y] = ode45(@myode, [0, 10], [30 * pi / 180, 45 * pi / 180, 0, 0]);

for i = 1:size(y,1)
    tau1(i)= -(K(1,1) * y(i, 1) + K(1, 2) * y(i, 2)+ K(1, 3) * y(i, 3) + K(1, 4) *y (i, 4));
    tau2(i)= -(K(2,1) * y(i, 1) + K(2, 2) * y(i, 2)+ K(2, 3) * y(i, 3) + K(2, 4) *y (i, 4));
end

figure;
plot(t, y(:, 1));
xlabel('Time')
ylabel('theta1')
figure;
plot(t, y(:, 2));
xlabel('Time')
ylabel('theta2')
figure;
plot(t, y(:, 3));
xlabel('Time')
ylabel('theta1d')
figure;
plot(t, y(:, 4));
xlabel('Time')
ylabel('theta2d')
figure;
plot(t, tau1);
xlabel('Time')
ylabel('Torque1')
figure;
plot(t, tau2);
xlabel('Time')
ylabel('Torque2')


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

    K = [20.5200, 4.8846, 4.9672, 1.4070;
         4.9505,  4.6809, 1.4598, 0.6141];

    % X is a cell so we need to convert it to integer
    X = [theta1, theta2, theta1_dot, theta2_dot]';
    U = -K * X;
    
    tau1 = U(1, :);
    tau2 = U(2, :);

    dX(1) = theta1_dot;
    dX(2) = theta2_dot;
    dX(3) = -(382000*tau1 - 382000*tau2 + 5433759*sin(theta1) + 171900*theta1_dot^2*sin(theta2) + 171900*theta2_dot^2*sin(theta2) - 600000*tau2*cos(theta2) - 2648700*sin(theta1 + theta2)*cos(theta2) + 343800*theta1_dot*theta2_dot*sin(theta2) + 270000*theta1_dot^2*cos(theta2)*sin(theta2))/(270000*cos(theta2)^2 - 491443);
    dX(4) = (1146000*tau1 - 6292000*tau2 - 22717017*sin(theta1 + theta2) + 16301277*sin(theta1) + 2831400*theta1_dot^2*sin(theta2) + 515700*theta2_dot^2*sin(theta2) + 25604100*cos(theta2)*sin(theta1) + 1800000*tau1*cos(theta2) - 3600000*tau2*cos(theta2) - 7946100*sin(theta1 + theta2)*cos(theta2) + 1031400*theta1_dot*theta2_dot*sin(theta2) + 1620000*theta1_dot^2*cos(theta2)*sin(theta2) + 810000*theta2_dot^2*cos(theta2)*sin(theta2) + 1620000*theta1_dot*theta2_dot*cos(theta2)*sin(theta2))/(3*(270000*cos(theta2)^2 - 491443));
  
end