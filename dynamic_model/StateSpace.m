%% Sate Sapce Representation

clear all

%% Define variable

syms T1 T2 'real'
syms theta1 theta2 theta1_d theta2_d theta1_dd theta2_dd 'real'
syms g m1 m2 l1 l2 r1 r2 I1 I2 'real'

%% Compute EOM

%[tau1, tau2] = TwoDOF();

%% Define State-Vector

X = sym('X', [4, 1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_d;
X(4) = theta2_d;

eq1 =  I1*theta1_dd + l1^2*m2*theta1_dd + m1*r1^2*theta1_dd + m2*r2^2*theta1_dd + m2*r2^2*theta2_dd - g*m1*r1*sin(theta1) - g*m2*r2*sin(theta1 + theta2) - g*l1*m2*sin(theta1) - l1*m2*r2*sin(theta2)*theta2_d^2 + 2*l1*m2*r2*cos(theta2)*theta1_dd + l1*m2*r2*cos(theta2)*theta2_dd - 2*l1*m2*r2*sin(theta2)*theta1_d*theta2_d - T1;
eq2 = I2*theta2_dd + m2*r2^2*theta1_dd + m2*r2^2*theta2_dd - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*sin(theta2)*theta1_d^2 + l1*m2*r2*cos(theta2)*theta1_dd - T2;

sol = solve([eq1 == 0, eq2 == 0], [theta1_dd, theta2_dd]);

%% ODE

[t, y] = ode45(@myode, [0, 10], [200 * pi / 180, 125 * pi / 180, 0, 0]);

%% Plotting

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
    [theta1, theta2, theta1_d, theta2_d] = deal(X{:});
    
    % I get an error with wrapTo2pi
    if abs(theta1) > 2 * pi
        theta1 = mod(theta1, 2 * pi);
    end

    if abs(theta2) > 2 * pi
        theta2 = mod(theta2, 2 * pi);
    end

    T1 = 0;
    T2 = 0;

    dX(1) = theta1_d;
    dX(2) = theta2_d;
    dX(3) = (I2*T1 + T1*m2*r2^2 - T2*m2*r2^2 + l1*m2^2*r2^3*theta1_d^2*sin(theta2) + l1*m2^2*r2^3*theta2_d^2*sin(theta2) + I2*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - T2*l1*m2*r2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_d*theta2_d*sin(theta2) + l1^2*m2^2*r2^2*theta1_d^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta2_d^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_d*theta2_d*sin(theta2))/(I1*I2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 + I2*m2*r2^2 - l1^2*m2^2*r2^2*cos(theta2)^2 + m1*m2*r1^2*r2^2 + 2*I2*l1*m2*r2*cos(theta2));
    dX(4) = -(T1*m2*r2^2 - T2*l1^2*m2 - T2*m1*r1^2 - I1*T2 - T2*m2*r2^2 + l1*m2^2*r2^3*theta1_d^2*sin(theta2) + l1^3*m2^2*r2*theta1_d^2*sin(theta2) + l1*m2^2*r2^3*theta2_d^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + T1*l1*m2*r2*cos(theta2) - 2*T2*l1*m2*r2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_d*theta2_d*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_d^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_d^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_d^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_d*theta2_d*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta1_d^2*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(I1*I2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 + I2*m2*r2^2 - l1^2*m2^2*r2^2*cos(theta2)^2 + m1*m2*r1^2*r2^2 + 2*I2*l1*m2*r2*cos(theta2));

end