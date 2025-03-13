clear all;
close all;
clc;
addpath("C:\Users\botto\OneDrive\Backup\Bertram lenovo\DTU\34746 Robust and Fault-tolerant Control\Exercise\SA Tool\SaTool_3_0100\sa_tool_3_0100\");
addpath("C:\Users\nicol\OneDrive\Documents\34746 - Robust & Fault Tolerent Control\Exercises\SA Tool\SaTool_3_0100\sa_tool_3_0100");


% The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;                    % Sampling period
sigma_meas = 0.0093*eye(3);     % Measurements covariance matrix

%% SA TOOL
% Inputs
syms u_1(t) u_2(t) % Symbolic input declearation
ST_input(1,:) = {u_1,u_2}; % symbolic variable
ST_input(2,:) = {'u_1(t)','u_2(t)'}; % LaTeX expression


% Measurements
syms y_1(t) y_2(t) y_3(t) % Symbolic measurement declearation
ST_meas(1,:) = {y_1,y_2,y_3}; % symbolic variable
ST_meas(2,:) = {'y_1(t)','y_2(t)','y_3(t)'}; % LaTeX expression


% Unknowns
syms Theta_1(t) dTheta_1(t) omega_1(t) domega_1(t) Theta_2(t) dTheta_2(t) omega_2(t) domega_2(t) Theta_3(t) dTheta_3(t) omega_3(t) domega_3(t) d % Symbolic unknowns declearation
ST_unknowns(1,:) = {Theta_1, dTheta_1, omega_1, domega_1, Theta_2, dTheta_2, omega_2, domega_2, Theta_3, dTheta_3, omega_3, domega_3, d}; % symbolic variable
ST_unknowns(2,:) = {'theta_1', 'dtheta_1', 'omega_1', 'domega_1', 'theta_2', 'dtheta_2', 'omega_2', 'domega_2', 'theta_3', 'dtheta_3', 'omega_3', 'domega_3', 'd'}; % LaTeX expression


% Parameters
syms J_1 J_2 J_3 k_1 k_2 b_1 b_2 b_3 % Symbolic parameter declearation
ST_parameters(1,:) = {J_1, J_2, J_3, k_1, k_2, b_1, b_2, b_3}; % symbolic variable
ST_parameters(2,:) = {'J_1', 'J_2', 'J_3', 'k_1', 'k_2', 'b_1', 'b_2', 'b_3'};% LaTeX expression


ST_cons(1,:) = {'c1','c2', 'c3','c4','c5','c6','m13','m14','m15','d7','d8', 'd9','d10','d11', 'd12'};
ST_cons(2,:) = {'c_1','c_2', 'c_3','c_4','c_5','c_6','m_{13}','m_{14}','m_{15}','d_7','d_8', 'd_9','d_{10}','d_{11}', 'd_{12}'};     % LaTeX   
ST_cons(3,:) = {...
    0 == dTheta_1 - omega_1, ...
    0 == J_1*domega_1 - u_1 + b_1*omega_1 + k_1*(Theta_1 - Theta_2) + d, ...
    0 == dTheta_2 - omega_2, ...
    0 == J_2 * domega_2 - u_2 + b_2*omega_2 + k_1*(Theta_2 - Theta_1) + k_2*(Theta_2 - Theta_3), ...
    0 == dTheta_3 - omega_3, ...
    0 == J_3*domega_3 + b_3*omega_3 + k_2*(Theta_3 - Theta_2), ...
    0 == y_1 - Theta_1, ...
    0 == y_2 - Theta_2, ...
    0 == y_3 - Theta_3, ...
    0 == dTheta_1 - diff(Theta_1,t), ...
    0 == domega_1 - diff(omega_1,t), ...
    0 == dTheta_2 - diff(Theta_2,t), ...
    0 == domega_2 - diff(omega_2,t), ...
    0 == dTheta_3 - diff(Theta_3,t), ...
    0 == domega_3 - diff(omega_3,t)};

    ST_canfail = [1:9];

    %husk {} omkring causale parameter, (differential)
    cons_oneway = {[],[],[],[],[],[],[],[],[],{Theta_1},{omega_1},{Theta_2},{omega_2},{Theta_3},{omega_3}};
    ST_domains = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

    ST_IncMat = sa_incidencematrix(ST_cons,ST_input,ST_meas,ST_unknowns,cons_oneway)

    % create system
ST_sys =...
    sa_create(ST_IncMat,ST_cons,...
    ST_input, ST_meas,ST_unknowns,...
    ST_domains, ST_canfail,ST_parameters)

% Find matchings
ST_sys=sa_match(ST_sys,'rank')

% Create report
%sa_report(ST_sys,'Assignment_A','pdf','analytic_expressions',false); 

% Create  matching
disp('Obtained matching:');
sa_disp(ST_sys, 't');

% Create Symbolic form
disp('Parity relation symbolic form:')
sa_disp(ST_sys, 's');

% Create analytical form
disp('Parity relation analytic form:')
sa_disp(ST_sys, 'a');

% Analystisk form in Laplase
syms s
syms y_1 y_2 y_3
syms u_1 u_2
a1 = (u_2 - b_2*y_2*s + k_1*y_1 - k_1*y_2 - k_2*y_2 + k_2*y_3)/J_2 - y_2*s^2;
a2 = (k_2*(y_2 - y_3) - b_3*y_3*s)/J_3 - y_3*s^2;

syms lambda

a1-lambda*a2
% Define lambda
solve(0==k_2/J_2+(lambda*k_2)/J_3+(lambda*b_3*s)/J_3-lambda*s^(2),lambda)
lambda = - (J_3*k_2) / (J_2*J_3*s^2 - J_2*b_3*s + J_2*k_2);

% Compute a3 = a1 - lambda * a2
a3 = simplify(a1 - lambda * a2);

% Solve for y3 in terms of other variables
y3_expr = solve(a2, y_3);

% Substitute y3 into a3 to eliminate it
a3_no_y3 = simplify(subs(a3, y_3, y3_expr));
a3 = a3_no_y3
disp('Simplified a3 without y3:');
disp(a3);
subs(subs(a3,y_1,0),y_2,0);
%% State space representation
load('ECP_values.mat');
load('ECP502Data.mat');
% Physical system parameters
J_1 = ECP_values(1);            % Disk 1 inertia kgm^2
J_2 = ECP_values(2);            % Disk 2 inertia kgm^2
J_3 = ECP_values(3);            % Disk 3 inertia kgm^2
k_1 = ECP_values(4);            % Shaft 1-2 stiffness Nm/rad
k_2 = ECP_values(5);            % Shaft 2-3 stiffness Nm/rad
b_1 = mean(ECP_values([6 7]));  % Disk 1 damping and friction Nms/rad
b_2 = mean(ECP_values([8 9]));  % Disk 2 damping and friction Nms/rad
b_3 = mean(ECP_values([10 11]));% Disk 3 damping and friction Nms/rad
T_Cp = ECP_values(12);          % Disk 1 Coulomb friction in positive direction
T_Cm = ECP_values(13);          % Disk 1 Coulomb friction in negative direction
atan_scale = 100;               % Sign approximation factor
w_th = 0.75;                    % Threshold angular velocity rad/s

A = [0 1 0 0 0 0
    -k_1/J_1 -b_1/J_1 k_1/J_1 0 0 0
    0 0 0 1 0 0
    k_1/J_2 0 -(k_1+k_2)/J_2 -b_2/J_2 k_2/J_2 0
    0 0 0 0 0 1
    0 0 k_2/J_3 0 -k_2/J_3 -b_3/J_3];
B = [0 0
    1/J_1 0
    0 0
    0 1/J_2
    0 0
    0 0];
C = [  1 0 0 0 0 0
       0 0 1 0 0 0
       0 0 0 0 1 0 ];
D = [0 0
     0 0
     0 0];

% Disturbance matrices
E_x = [ 0 
        1/J_1 
        0 
        0 
        0 
        0 ];
E_y = [ 0 
        0 
        0 ];

% Actuator and system fault matrix (6x5) %f ={row = size of input output, cols = size of state] 
F_x = [0 0 0 0 0  
       0 0 0 1/J_1 0   
       0 0 0 0 0   
       0 1/J_2 0 0 1/J_2   
       0 0 0 0 0  
       0 0 0 0 0 ];


% Sensor fault matrices Skal være 3x5
F_y = [ 1 0 0 0 0
        0 1 0 0 0
        0 0 1 0 0 ];
%%%%%%%% SKIP THAT WHEN SIMULATING IN OPEN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete time
sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, T_s, 'zoh');
F_d = sys_d.A;
G_d = sys_d.B;

load('ECP502Data.mat');
ymeas_ts = timeseries(y_meas, t);
u1_ts = timeseries(u_1, t);
u2_ts = timeseries(u_2, t);

% State-feedback LQR design
Q_c = diag([2 0 2 0 2.5 0.0024]);
R_c = diag([10 10]);
K_c = dlqr(F_d, G_d, Q_c, R_c);

% Scaling of reference
C_3 = [0 0 0 0 1 0];
C_ref = pinv(C_3 * ((eye(size(F_d,1)) - F_d + G_d * K_c) \ G_d) * K_c);

% Kalman filter with friction estimation - DO NOT MODIFY
F_aug = [F_d G_d(:,1);zeros(1,6) 1];
G_aug = [G_d;0 0];
C_aug = [C zeros(3,1)];
% Kalman gain
L_aug = dlqe(F_aug, eye(7), C_aug, 1e-3 * eye(7), sigma_meas(1,1).^2 * eye(3));
L_o = L_aug(1:6,:);
L_d = L_aug(7,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the one we run with
wn = 15;  % Cutoff frequency (rad/s), chosen below the system's resonant frequencies
zeta = 0.707;  % Damping ratio (Butterworth design)

% Define second-order low-pass filter
s = tf('s');
H = (wn^2) / (s^2 + 2*zeta*wn*s + wn^2);
% Display the transfer function
disp('Adjusted Second-Order Low-Pass Filter Transfer Function:');
H

H_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2]
H_ru = [0 tf(1); 0 0]

sys_y = H*H_ry
sys_u = H*H_ru

sys = sys_u+sys_u
T_s = 4e-3; % Sampling period (4 ms)
dsys = c2d(sys,T_s, 'zoh')


% Display the discretized model
disp(dsys);

figure;
step(dsys);
title('Step Response of Discretized Residual System');
xlabel('Time (seconds)');
ylabel('Response');
grid on;



%% Simulation Opgave 3
% Extract variables
ECP_t = t;       % Time vector
ECP_u_1 = u_1;   % First control input
ECP_u_2 = u_2;   % Second control input
ECP_y_meas = y_meas; % Measured output

% Create figure with subplots
%figure;
%simIn = Simulink.SimulationInput("ActuatorSystem");
%simIn = simIn.setExternalInput([ECP_t, ECP_u_1, ECP_u_2]); % Time + Inputs
simTime = 50;                   % Simulation duration in seconds
f_m_time = 8.5;                 % Sensor fault occurence time
u_1 = [t u_1]
u_2 = [t u_2]
y_meas =[t y_meas]
simout = sim("threeDiskOscillatorRig_v2");
r1 = simout.residual_out.signals.values(:,1); 
r2 = simout.residual_out.signals.values(:,2);

% First subplot: r_1 and r_2 vs time
figure;
subplot(2,1,1);
plot(t, r1, 'b', 'LineWidth', 1.5); hold on;
plot(t, r2, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('residual magnitude');
legend('r_1', 'r_2');
title('residuals vs Time');
grid on;


% First subplot: u_1 and u_2 vs time
figure;
subplot(2,1,1);
plot(ECP_t, ECP_u_1, 'b', 'LineWidth', 1.5); hold on;
plot(ECP_t, ECP_u_2, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Inputs');
legend('u_1', 'u_2');
title('ECP Control Inputs vs Time');
grid on;

% Second subplot: y_meas vs time
subplot(2,1,2);
plot(ECP_t, ECP_y_meas, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Measured Output');
title('ECP Measured Output vs Time');
grid on;
% Display the figure
hold off;

%%
%%Opgave4  
%from slide 11,12,27, lecture 5:
syms s
H_yu = C*inv((s*eye(size(A))-A))*B+D;
H_yd = C*(inv(s*eye(size(A))-A))*E_x+E_y;
H_yf = C*(inv(s*eye(size(A))-A))*F_x+F_y;

% V_ry = H_ry as we defined earlier however this time se use the symbolic s. 
V_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2]; 
H_rf= V_ry*H_yf;

%lecture 5 slide 23
Rank([H_yd H_yf] >rank(H_yd))

%% Strong and weak detectability
H_rf = tf(0);

%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
h = 0;                  % Put the threshold from GLR here

%% Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
va_eig_d = [];  % Discrete time eigenvalues
va_eig = log(va_eig_d)/T_s;     % Continuous time eigenvalues
% Then discretise your VA

B_change = [1 0;0 0];

%% Simulation for sensor fault (f_u = 0)
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m_time = 8.5;                 % Sensor fault occurence time
sim('threeDiskOscillatorRig');

%% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
sim('threeDiskOscillatorRig');

%% Plot settings
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);

