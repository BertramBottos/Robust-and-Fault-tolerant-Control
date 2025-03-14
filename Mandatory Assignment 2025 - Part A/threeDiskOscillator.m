clear all;
close all;
clc;
addpath("C:\Users\botto\OneDrive\Backup\Bertram lenovo\DTU\34746 Robust and Fault-tolerant Control\Exercise\SA Tool\SaTool_3_0100\sa_tool_3_0100\");
addpath("C:\Users\nicol\OneDrive\Documents\34746 - Robust & Fault Tolerent Control\Exercises\SA Tool\SaTool_3_0100\sa_tool_3_0100");


% The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;                    % Sampling period
sigma_meas = (0.0093^2)*eye(3);     % Measurements covariance matrix

%% SA TOOL
% Inputus
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

%% SA Tool done
clear all;
close all;
clc;

% Analystisk form
syms y_1 y_2 y_3 u_1 u_2 s b_1 b_2 J_1 J_2 J_3 k_1 k_2 b_3
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

partfrac(subs(subs(a3,y_1,0),y_2,0))  
partfrac(subs(subs(a3,y_2,0),y_1,0))


J_1 = 0.0025; % kgm2 Bottom disk moment of inertia
J_2 = 0.0018; %kgm2 Middle disk moment of inertia
J_3 = 0.0018; %kgm2 Top disk moment of inertia
k_1 = 2.7;% Nmrad−1 Stiffness of the bottom shaft
k_2 = 2.6;% Nmrad−1 Stiffness of the middle shaft
b_1 = 0.0029;% Nmsrad−1 Damping/friction on the bottom disk
b_2 = 0.0002;% Nmsrad−1 Damping/friction on the middle disk
b_3 =  0.00015;% Nmsrad−1 Damping/friction on the top disk
T_s = 0.004;                    % Sampling period


% Laplace 
H1_s = (1/J_2) * (1 / (s^2 + (b_2/J_2) * s + (k_1 + k_2)/J_2));
H2_s = (1/J_3) * (k_2 / (s^2 + (b_3/J_3) * s + k_2/J_3));
%H3_s

[H1_s_num, H1_s_den] = numden(H1_s);
[H2_s_num, H2_s_den] = numden(H2_s);

% Transfer Function
H1_tf = tf(sym2poly(H1_s_num), sym2poly(H1_s_den));
H2_tf = tf(sym2poly(H2_s_num), sym2poly(H2_s_den));

% Discreticed with a sampling period of T_s = 4 ms.
H1_d = c2d(H1_tf, T_s, 'zoh');  % Discretizing
H2_d = c2d(H2_tf, T_s, 'zoh');  % Discretizing


disp('Discrete-time Residual Generator H1(z):');
H1_d

disp('Discrete-time Residual Generator H2(z):');
H2_d



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

% Disturbance matrices               %why is it not just 1??? 
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
sigma_meas = (0.0093^2)*eye(3);     % Measurements covariance matrix
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

%H_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
%   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2]
%H_ru = [0 tf(1); 0 0]

H_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2;
   (k_1)/J_2 -((k_1*k_2+b_2*k_2*s+b_3*k_1*s+b_3*k_2*s+J_2*J_3*s^4+J_2*b_3*s^3+J_3*b_2*s^3+J_2*k_2*s^2+J_3*k_1*s^2+J_3*k_2*s^2+b_2*b_3*s^2))/(J_2*(J_3*s^2+b_3*s+k_2)) 0]
H_ru = [0 tf(1); 0 0; 0 1/J_2]

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
u_1_in = [t u_1];
u_2_in = [t u_2];
y_meas_in =[t y_meas];
simout = sim("threeDiskOscillatorRig_v2");
r1 = simout.r.signals.values(:,1); 
r2 = simout.r.signals.values(:,2);
r3 = simout.r.signals.values(:,3);


%%
%creating plots for assignmnet 3
% First subplot: r_1 and r_2 vs time
tsim_nested = simout.r.time; %dont know dont care why t = simulated time
figure;
subplot(2,1,1);
plot(tsim_nested,r1, 'b', 'LineWidth', 1.5); hold on;
plot(tsim_nested,r2, 'r', 'LineWidth', 1.5); hold on; 
plot(tsim_nested,r3, 'r', 'LineWidth', 1.5,'Color',"g");
xlabel('Time (s)');
ylabel('residual magnitude');
legend('r_1', 'r_2','r3');
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
%%Opgave4  - strong and weak fault detectablity  
% PLEASE NOTE we are
%%working with the unfiltered residuals. (because that is whats
%%written in the slides) & it dosent affect detectablity.

%from slide 11,12,27, lecture 5:
syms y_1 y_2 y_3 u_1 u_2 s b_1 b_2 J_1 J_2 J_3 k_1 k_2 b_3 s %can uncomment if you want to see numerical solution
InputSyms()

H_yu = C*inv((s*eye(size(A))-A))*B+D;
H_yd = C*(inv(s*eye(size(A))-A))*E_x+E_y;
H_yf = C*(inv(s*eye(size(A))-A))*F_x+F_y;



% V_ry = H_ry as we defined earlier however this time se use the symbolic s. 
V_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2]; 
H_rf= simplify(V_ry*H_yf);

H_faulty = [H_yu H_yd;eye(size(H_yu,2),size(H_yu,2)) zeros(size(H_yu,2),size(H_yd,2))]; %see slide 13 lecture 5 
F=simplify(null(H_faulty')');  %by design which means r(s) only depend on Vry*Hyf*f(s) 
%check if F is designed correct
simplify(F*H_faulty)

%lecture 5 slide 23/25
for i=1:size(H_yf,2)
if(rank([H_yd H_yf(:,i)]) >rank(H_yd) )  
    fprintf("faut %d is weakly detectable\n", i)
else 
   fprintf("faut %d is NOT weakly detectable\n", i)
end
end



% is it good enough that we detect it in one residual and
%that means the entire i'th fault is strongly detectible or de we need to
%check for r1 and r2 seperately???
for j=1:size(H_yf,2)
 strong_detectability = subs(simplify( F * [H_yf(:,j); zeros(2,1)]), s, 0);

 %checking if fault is strongly detetible on any residual
if (strong_detectability(1) ~= 0) || (strong_detectability(2) ~= 0)
    fprintf("Fault %d is strongly detectable\n", i);
else
    fprintf("Fault %d is NOT strongly detectable\n", i);
end

end


%% Opgave 5 
%r1_sym = V_ry(1,:);
%r2_sym = V_ry(2,:);

%after introducting F we have: (non sensitive to y) 
r_f = V_ry*H_yf;
%okay so this is not the correct residual generator they are asking for.

%earlier we had a very nice one "sys_y" which is sensitive to y. (should i
%take the filtered one) ??? ) 

%sys_u_padded = [sys_u, zeros(2,1)];  % Adding an extra zero column
H_ru_padded = [H_ru,zeros(3,1)]

% Add the filtered transfer functions
%sys = sys_y + sys_u_padded;
systest = H_ry+H_ru_padded

ss(sys)
A1 = ss(systest).A;
B1 =ss(systest).B
C1 =ss(systest).C
D1 =  ss(systest).D
Qx = lyap(A1,B1*sigma_meas*B1');
Qy = C1*Qx*C1'+D1*sigma_meas*D1';

sigma_r1 = Qy(1,1)
sigma_r2 = Qy(2,2)


%GLR 
M=20; %initial guess of window size
mu_0 = 0;
sigma = sqrt(var(r1));
[g,mu_1,idx] = GLR(r1,M,mu_0,sigma)



%GLR tuning 
P_F=0.0001;
P_M=0.01;
%P_M = probability of missed detection 
%P_D = probabilty of detection (im not sure if its the same). 
P_D=1-P_M
%finding threshold h
k = 1/2;
theta = 2;
% Density function expression
syms zz gg xx gg1;
pd_zz = 1/(theta^k*gamma(k))*zz^(k-1)*exp(-zz/theta);
p_zz = int(pd_zz,zz,2*gg,Inf);  % Integrate over the probability space
eq_1 = P_F - p_zz == 0;  % Equation to be solved
h = double(vpasolve(eq_1,gg));  % % Put the threshold from GLR here
pf_calc = double(int(pd_zz,zz,2*h,Inf)); %calculated probability of false alarm slide 24 lecture 7
% check
disp(P_F - pf_calc); %check error of vpa solve
disp(h - chi2inv(1 - P_F,1)/2); % Compare with the other method
%h = chi2inv(1 - P_F,1)/2;
fault = (g>h)*40;

%glr plot
figure;
aa = subplot(2,1,1);
plot(tsim_nested, r1, 'DisplayName', 'r1'); % Plot r1
hold on;
plot(tsim_nested, mu_1, 'DisplayName', 'mu_1'); % Plot mu_1
hold off;
grid on;
title("r1 with mean")
ylabel('$r(t)$', 'Interpreter', 'latex'); % Use latex for labels
legend('show'); % Display legend

% Plotting the second subplot (GLR plot g and h)
bb = subplot(2,1,2);
plot(tsim_nested, g, 'DisplayName', 'g(t)'); % Plot GLR statistic
hold on;
plot(tsim_nested, ones(size(tsim_nested,1),1)*h, 'DisplayName', 'h'); % Plot h (ensure h is the same size as g)
hold on; 
plot(tsim_nested,fault,'DisplayName','fault')
hold off;
title('GLR PLOT');
xlabel('$t$ in s', 'Interpreter', 'latex'); % Latex formatting
ylabel('$g(t)$', 'Interpreter', 'latex'); % Latex formatting
legend('show'); % Display legend
grid on;


%%
%finding window size slide 27 week 7
mu_1 = 75; %% looked at r plot - we need to look into this !!!!!!!!!!!!!!!!!!!!!!!!!!!
mu_0 = 0;
lambda_1 =(gg1*(mu_1-mu_0)^2)/sigma^2;
bes = besseli(-1/2, sqrt(lambda_1 * xx));
pd_xx = (1/2) * (xx/lambda_1)^(-1/4) * exp(-(xx+lambda_1)/2)*bes; 
p_xx = int(pd_xx,xx,2*h,Inf)
eq_2 = P_D-p_xx == 0
M = double(vpasolve(eq_2,gg1))
%pd_calc = double(int(pd_xx,xx,2*h,Inf)); %calculated probability of false alarm slide 24 lecture 7
%check
%disp(P_D - pd_calc)
%disp(P_d-ncx2cdf(2*h,1,lambda_1))

f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
                 


%%
Atest= {1, s, 3; 1, 1, 1; 1, 1, 1};
Btest = eye(3)
lyap(Atest,Btest*sigma_meas*Btest')

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

