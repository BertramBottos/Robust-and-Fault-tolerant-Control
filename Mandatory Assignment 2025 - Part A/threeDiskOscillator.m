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

% Analystisk form (la place transformed)
syms y_1 y_2 y_3 u_1 u_2 s b_1 b_2 J_1 J_2 J_3 k_1 k_2 b_3
a1 = (u_2 - b_2*y_2*s + k_1*y_1 - k_1*y_2 - k_2*y_2 + k_2*y_3)/J_2 - y_2*s^2;
a2 = (k_2*(y_2 - y_3) - b_3*y_3*s)/J_3 - y_3*s^2;

syms lambda

a1-lambda*a2;
% Define lambda
solve(0==k_2/J_2+(lambda*k_2)/J_3+(lambda*b_3*s)/J_3-lambda*s^(2),lambda)
lambda = - (J_3*k_2) / (J_2*J_3*s^2 - J_2*b_3*s + J_2*k_2);

% Compute a3 = a1 - lambda * a2
a3 = simplify(a1 - lambda * a2);

% Solve for y3 in terms of other variables
y3_expr = solve(a2, y_3);

% Substitute y3 into a3 to eliminate it
a3_no_y3 = simplify(subs(a3, y_3, y3_expr));
a3 = a3_no_y3;
disp('Simplified a3 without y3:');
disp(a3);
subs(subs(a3,y_1,0),y_2,0);

partfrac(subs(subs(a3,y_1,0),y_2,0));  
partfrac(subs(subs(a3,y_2,0),y_1,0));


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
%% This is the one we run with
%Plot settings
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);

wn = 15;  % Cutoff frequency (rad/s), chosen below the system's resonant frequencies
zeta = 0.707;  % Damping ratio (Butterworth design)

% Define second-order low-pass filter
s = tf('s');
%H = (wn^2) / (s^2 + 2*zeta*wn*s + wn^2)* (wn^2) / (s^2 + 2*zeta*wn*s + wn^2);
H = 1 / (s^2 + 2*zeta*wn*s + wn^2)* (wn^2) / (s^2 + 2*zeta*wn*s + wn^2)
% Display the transfer function
disp('Adjusted Second-Order Low-Pass Filter Transfer Function:');
H


H_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2;
   (k_1)/J_2 -((k_1*k_2+b_2*k_2*s+b_3*k_1*s+b_3*k_2*s+J_2*J_3*s^4+J_2*b_3*s^3+J_3*b_2*s^3+J_2*k_2*s^2+J_3*k_1*s^2+J_3*k_2*s^2+b_2*b_3*s^2))/(J_2*(J_3*s^2+b_3*s+k_2)) 0];
H_ru = [0 tf(1/J_2); 0 0; 0 tf(1/J_2)];
H_ru_padded = [H_ru,zeros(3,1)];

sys_y = H*H_ry;
sys_u = H*H_ru;
sys_u_normal = H*H_ru; %need thsi for running our simulink model 

sys = [sys_y,sys_u]
T_s = 4e-3; % Sampling period (4 ms)
dsys = c2d(sys,T_s, 'zoh')
dsys_y = c2d(sys_y,T_s,'zoh');


% Display the discretized model
disp(dsys);

figure;
step(dsys);
title('Step Response of Discretized System');
xlabel('Time (seconds)');
ylabel('Response');
grid on;



%% Simulation Opgave 3
x_0 = [0;0;0;0;0;0];            % Initial conditions

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
u1_test = [zeros(1000,1); ones((size(u_1,1)-1000),1)]; %after 4 sec
u2_test = [zeros(2000,1); ones((size(u_1,1)-2000),1)]; %after 8 sec
u_1_in = [t u1_test];
u_2_in = [t u2_test];
%y_meas_in =[t y_meas];
y_meas_in = [t,zeros(size(y_meas,1),3)];
simout = sim("threeDiskOscillatorRig_v2");
r1 = simout.r.signals.values(:,1); 
r2 = simout.r.signals.values(:,2);
r3 = simout.r.signals.values(:,3);



%creating plots for assignmnet 3
% First subplot: r_1 and r_2 vs time
tsim_nested = simout.r.time; %dont know dont care why t = simulated time
figure;
subplot(2,1,1);
plot(tsim_nested,r1, 'b', 'LineWidth', 1.5); hold on;
plot(tsim_nested,r2, 'r', 'LineWidth', 1.5); hold on; 
plot(tsim_nested,r3, 'r', 'LineWidth', 1.5,'Color',"g");
xlabel('Time (s)');
ylabel('residual ');
legend('r_1', 'r_2','r3');
title('residuals');
grid on;
subplot(2,1,2);
plot(ECP_t, u1_test, 'b', 'LineWidth', 1.5); hold on;
plot(ECP_t, u2_test, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Inputs');
legend('u_1', 'u_2');
title('ECP Control Inputs');
grid on;
saveas(gcf, fullfile('fig', 'r_insensitive_to_U.png'));

u_1_in = [t u_1];
u_2_in = [t u_2];
y_meas_in =[t y_meas];
%y_meas_in = [t,zeros(size(y_meas,1),3)];
simout = sim("threeDiskOscillatorRig_v2");
r1 = simout.r.signals.values(:,1); 
r2 = simout.r.signals.values(:,2);
r3 = simout.r.signals.values(:,3);

% First subplot: u_1 and u_2 vs time
figure;
subplot(2,1,1);
plot(tsim_nested,r1, 'b', 'LineWidth', 1.5); hold on;
plot(tsim_nested,r2, 'r', 'LineWidth', 1.5); hold on; 
plot(tsim_nested,r3, 'r', 'LineWidth', 1.5,'Color',"g");
xlabel('Time (s)');
ylabel('residuals');
legend('r_1', 'r_2','r3');
title('residuals');
grid on;

% Second subplot: y_meas vs time
subplot(2,1,2);
%plot(ECP_t, zeros(size(y_meas,1),3), 'k', 'LineWidth', 1.5);
plot(ECP_t, y_meas, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Output y');
title('ECP Measured Output');
grid on;
% Display the figure
hold off;
saveas(gcf, fullfile('fig', 'R_sensitive_to_Y.png'));

%%
%%Opgave4  - strong and weak fault detectablity  
% PLEASE NOTE we are
%%working with the unfiltered residuals. (because that is whats
%%written in the slides) & it dosent affect detectablity.

%from slide 11,12,27, lecture 5:
%redefine for symbolic values
syms y_1 y_2 y_3 u_1 u_2 s b_1 b_2 J_1 J_2 J_3 k_1 k_2 b_3 s %can uncomment if you want to see numerical solution
H_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2;
   (k_1)/J_2 -((k_1*k_2+b_2*k_2*s+b_3*k_1*s+b_3*k_2*s+J_2*J_3*s^4+J_2*b_3*s^3+J_3*b_2*s^3+J_2*k_2*s^2+J_3*k_1*s^2+J_3*k_2*s^2+b_2*b_3*s^2))/(J_2*(J_3*s^2+b_3*s+k_2)) 0];
H_ru = [0 1/J_2; 0 0; 0 1/J_2];

% V_ry = H_ry as we defined earlier however this time se use the symbolic s. 
H_yu = C*inv((s*eye(size(A))-A))*B+D;
H_yd = C*(inv(s*eye(size(A))-A))*E_x+E_y;
H_yf = C*(inv(s*eye(size(A))-A))*F_x+F_y;

V_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2]; 
H_rf= simplify(V_ry*H_yf);

%lecture 5 slide 23/25
for i=1:size(H_yf,2)
if(rank([H_yd H_yf(:,i)]) >rank(H_yd) )  
    fprintf("faut %d is weakly detectable\n", i)
else 
   fprintf("faut %d is NOT weakly detectable\n", i)
end
end

%s = tf('s');
H_yu = C*inv((s*eye(size(A))-A))*B+D;
H_yd = C*(inv(s*eye(size(A))-A))*E_x+E_y;
H_yf = C*(inv(s*eye(size(A))-A))*F_x+F_y;





%H_faulty = [H_yu H_yd;eye(size(H_yu,2),size(H_yu,2)) zeros(size(H_yu,2),size(H_yd,2))]; %see slide 13 lecture 5 
%F=simplify(null(H_faulty')');  %by design which means r(s) only depend on Vry*Hyf*f(s) 
F = [H_ry H_ru];
%check if F is designed correct
%simplify(F*H_faulty)



% is it good enough that we detect it in one residual and
%that means the entire i'th fault is strongly detectible or de we need to
%check for r1 and r2 seperately???

for j=1:size(H_yf,2)
 strong_detectability = subs(simplify(F * [H_yf(:,j); zeros(2,1)]), s, 0);

 %checking if fault is strongly detetible on any residual
if (strong_detectability(1) ~= 0) || (strong_detectability(2) ~= 0)
    fprintf("Fault %d is strongly detectable\n", j);
else
    fprintf("Fault %d is NOT strongly detectable\n", j);
end

end

%% Opgave 5 
[J_1, J_2, J_3, k_1, k_2, b_1, b_2, b_3, T_Cp, T_Cm] = InputValues();
sigma_meas = (0.0093^2)*eye(3);     % Measurements covariance matrix
%r1_sym = V_ry(1,:);
%r2_sym = V_ry(2,:);

%after introducting F we have: (non sensitive to y) 
r_f = V_ry*H_yf;
%okay so this is not the correct residual generator they are asking for.

%earlier we had a very nice one "sys_y" which is sensitive to y. (should i
%take the filtered one) ??? ) 

%sys_u_padded = [sys_u, zeros(2,1)];  % Adding an extra zero column
H_ru_padded = [H_ru,zeros(3,1)];

% Add the filtered transfer functions
%sys = sys_y + sys_u_padded;
sys_y_2 = H_ry;
sys_u_2 = H_ru_padded;
%sys_u_normal = H*H_ru; %need thsi for running our simulink model 

sys_2 = sys_y_2+sys_u_2
T_s = 4e-3; % Sampling period (4 ms)
%dsys_2 = c2d(dsys,T_s, 'zoh')
systest = H_ry+H_ru_padded;
%systest = dsys


F1 = ss(dsys_y(2,:)).A;
G1 =ss(dsys_y(2,:)).B;
C1 =ss(dsys_y(2,:)).C;
D1 =  ss(dsys_y(2,:)).D;
Qx = dlyap(F1,G1*sigma_meas*G1');
Qy = C1*Qx*C1'+D1*sigma_meas*D1';

%sigma_r1 = sqrt(Qy(1,1));
sigma_r2 = sqrt(Qy);
var_diff_precent = ((var(r2)-Qy)/var(r2))*100
%sigma_r3 = sqrt(Qy(3,3));


%GLR 
M_window=20; %initial guess of window size
mu_0 = 0;
%sigma = sqrt(var(r2)); %residuals from the given data. 
[g,mu_1,idx] = GLR(r2,M_window,mu_0,sigma_r2);



%GLR tuning 
P_F=0.0001;
P_M=0.01;
%P_M = probability of missed detection 
%P_D = probabilty of detection (im not sure if its the same). 
P_D=1-P_M;
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
fault = (g>h)*400;

%glr plot
figure;
aa = subplot(2,1,1);
plot(tsim_nested, r2, 'DisplayName', 'r1'); % Plot r1
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
saveas(gcf, fullfile('fig', 'GLR_theory.png'));

%finding window size
mu_1_ex = mu_1(10000,1); % when fault have occured
M= 1:70;
lambda_1 = (M * (mu_1_ex - mu_0)^2) / (var(r2));
PD_calc = 1 - ncx2cdf(2*h, 1, lambda_1);
M_window = find(PD_calc > P_D, 1);
%% GLR simulation Question 5.3) 
B_change = [1 0;0 0];
% Simulation for no actuator or sensor fault (f_u = 0, fm=u =0)
x_0 = [0;0;0;0;0;0];            % Initial conditions
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m = [0;0;0];     % Sensor fault vector (added to [y1;y2;y3])
f_m_time = 8.5;                 % Sensor fault occurence time
r_in0 = [zeros(2000,1); -25*ones(size(r2,1)-2000,1)];
r_in = [tsim_nested r_in0];
glr_sim = sim('threeDiskOscillatorRig_v3_1');

figure
subplot(2,1,1);  % First subplot in a 2x1 grid
plot(glr_sim.tout,glr_sim.GLR.Data, 'b');  % Plot GLR data in blue
hold on;  % Hold the plot to add another plot
yline(h, 'r--', 'LineWidth', 2);  % Add a horizontal line at y = h (red dashed line)
xlabel('time');
ylabel('GLR Data');
title('GLR Data with threshold h');
grid on;

% Second subplot: Plot the signal r_in
subplot(2,1,2);  % Second subplot in a 2x1 grid
plot(tsim_nested,r_in(:,2), 'g');  % Plot r_in in green
xlabel('time');
ylabel('r_{in}');
title('residual input');
grid on;
saveas(gcf, fullfile('fig', 'GLR_sim.png'));

%% GLR implementation question 6 
B_change = [1 0;0 0];
% Simulation for no actuator or sensor fault (f_u = 0, fm=u =0)
x_0 = [0;0;0;0;0;0];            % Initial conditions
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m = [0;0;0];     % Sensor fault vector (added to [y1;y2;y3])
f_m_time = 8.5;                 % Sensor fault occurence time
glr_imp = sim('threeDiskOscillatorRig_v3');
r2 = glr_imp.r.signals.values(:,2);

figure
subplot(2,1,1);  % First subplot in a 2x1 grid
plot(glr_imp.GLR.Data, 'b');  % Plot GLR data in blue
hold on;  % Hold the plot to add another plot
yline(h, 'r--', 'LineWidth', 2);  % Add a horizontal line at y = h (red dashed line)
xlabel('Index');
ylabel('GLR Data');
title('GLR Data with threshold h');
grid on;

% Second subplot: Plot the signal r_in
subplot(2,1,2);  % Second subplot in a 2x1 grid
plot(r2, 'g');  % Plot r_in in green
xlabel('Index');
ylabel('r2');
title('residual input');
grid on;
saveas(gcf, fullfile('fig', 'GLR_implementation.png'));

%%
%finding window size slide 27 week 7
%%mu_1 = 75; %% looked at r plot - we need to look into this !!!!!!!!!!!!!!!!!!!!!!!!!!!
%%mu_0 = 0;
%lambda_1 =(gg1*(mu_1-mu_0)^2)/sigma^2;
%bes = besseli(-1/2, sqrt(lambda_1 * xx));
%pd_xx = (1/2) * (xx/lambda_1)^(-1/4) * exp(-(xx+lambda_1)/2)*bes; 
%p_xx = int(pd_xx,xx,2*h,Inf)
%eq_2 = P_D-p_xx == 0
%M = double(vpasolve(eq_2,gg1))
%%pd_calc = double(int(pd_xx,xx,2*h,Inf)); %calculated probability of false alarm slide 24 lecture 7
%check
%disp(P_D - pd_calc)
%disp(P_d-ncx2cdf(2*h,1,lambda_1))


%%
%%%%%%%% SKIP THAT WHEN SIMULATING IN OPEN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-feedback LQR design
% Discrete time
sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, T_s, 'zoh');
F_d = sys_d.A;
G_d = sys_d.B;
C_d = sys_d.C;

load('ECP502Data.mat');
ymeas_ts = timeseries(y_meas, t);
u1_ts = timeseries(u_1, t);
u2_ts = timeseries(u_2, t);

% State-feedback dLQR design
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
L_aug = dlqe(F_aug,eye(7),C_aug,1e-3*eye(7),sigma_meas(1,1).^2*eye(3));
L_o = L_aug(1:6,:);
L_d = L_aug(7,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% LQR simulation
B_change = [1 0;0 0];
% Simulation for no actuator or sensor fault (f_u = 0, fm=u =0)
x_0 = [0;0;0;0;0;0];            % Initial conditions
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m = [0;0;0];     % Sensor fault vector (added to [y1;y2;y3])
f_m_time = 8.5;                 % Sensor fault occurence time
testout = sim('threeDiskOscillatorRig_v4');

theta_ref = testout.theta_ref.signals.values;
r1 = testout.r1.signals.values(:,1);
r2 = testout.r1.signals.values(:,2);
r3 = testout.r1.signals.values(:,3);
y1_nofault = testout.y_real.Data(:,1);
y2_nofault = testout.y_real.Data(:,2);
y3_nofault = testout.y_real.Data(:,3);


 %Time Vector (if available)
t = testout.r1.time;  % Assuming all signals share the same time vector

%Plot
figure;
plot(t, r1, 'r', 'LineWidth', 1.5); hold on;
plot(t, r2, 'g', 'LineWidth', 1.5);
plot(t, r3, 'b', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Signal Value');
title('Reidual Signals no fault present');
legend('r1', 'r2', 'r3');
grid on;

% Simulation for actuator fault fault (f_u = -0.1) no sensor fault
x_0 = [0;0;0;0;0;0];            % Initial conditions
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;-0.35];                 % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % enable/disable VA meachanism (0 = disable) 
f_m = [0;0;0];     % Sensor fault vector (added to [y1;y2;y3])
f_m_time = 8.5;                 % Sensor fault occurence time
sim_fault = sim('threeDiskOscillatorRig_v4');

y1_wfault = sim_fault.y_real.Data(:,1);
y2_wfault = sim_fault.y_real.Data(:,2);
y3_wfault = sim_fault.y_real.Data(:,3);


%Plot
figure;
plot(t, y1_nofault, 'r--', 'LineWidth', 3.5);  hold on;
plot(t, y2_nofault, 'g--', 'LineWidth', 3.5); 
plot(t, y3_nofault, 'b--', 'LineWidth', 3.5);  
plot(t, theta_ref, 'black', 'LineWidth', 3.5); 
hold off;
xlabel('Time (s)');
ylabel('Signal Value');
title('Output response wiht no fault');
legend('y1', 'y2', 'y3','theta ref');
grid on;
set(gcf, 'Position', get(0, 'Screensize'))
 saveas(gcf, fullfile('fig', 'LQR_nofault.png'));

figure;
plot(t, y1_wfault, 'r', 'LineWidth', 3.5); hold on;
plot(t, y2_wfault, 'g', 'LineWidth', 3.5);
plot(t, y3_wfault, 'b', 'LineWidth', 3.5);
%plot(t, y1_nofault, 'r--', 'LineWidth', 3.5);  % Dashed line for nofault
%plot(t, y2_nofault, 'g--', 'LineWidth', 3.5);  % Dashed line for nofault
%plot(t, y3_nofault, 'b--', 'LineWidth', 3.5);  % Dashed line for nofault
plot(t, theta_ref, 'black', 'LineWidth', 3.5)
hold off;
xlabel('Time (s)');
ylabel('Signal Value');
title('Faulty Output Response');
%legend('y1_{f}', 'y2_{f}', 'y3_{f}', 'y1', 'y2', 'y3','theta ref');
legend('y1', 'y2', 'y3','theta ref');
grid on;
set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, fullfile('fig', 'LQR_fault.png'));



y1_diff = y1_wfault - y1_nofault;
y2_diff = y2_wfault - y2_nofault;
y3_diff = y3_wfault - y3_nofault;

figure;
plot(t, y1_diff, 'r', 'LineWidth', 1.5); hold on;
plot(t, y2_diff, 'g', 'LineWidth', 1.5);
plot(t, y3_diff, 'b', 'LineWidth', 1.5);
hold off;

xlabel('Time (s)');
ylabel('Difference in Signal Value');
title('Difference between Fault and No Fault Signals');
legend('y1_{fault difference}', 'y2_{fault difference}', 'y3_{fault difference}');
grid on;

%% Virtual actuator design
Mc = ctrb(F_d,G_d);

if (rank(Mc) == size(F_d,2))
    disp('System is controllable');
else
    disp('System is not controllable');
end

B_f = [B(:,1) zeros(size(B(:,2)))];
G_f = [G_d(:,1) zeros(size(G_d(:,2)))];
if (rank(B_f) == rank([B B_f]))
    disp('Perfect static matching for actuator fault');
else
    disp('Imperfect static matching for actuator fault');
end

if (rank(G_f) == rank([G_d G_f]))
    disp('Perfect static matching for actuator fault');
else
    disp('Imperfect static matching for actuator fault');
end

if (rank(ctrb(A,B_f)) == size(A,2))
    disp('Faulty system is controllable');
else
    disp('Faulty system is not controllable');
end

% Continuous time
va_eig = log(eig(F_d - G_d*K_c))/T_s;
M = place(A, B_f, va_eig);
A_D = A - B_f*M;
N_D = pinv(B_f)*B;
B_D = B - B_f*N_D;
C_D = C;

% Discrete time
va_eig_d = exp(va_eig*T_s);
M_d = place(F_d, G_f, va_eig_d);
F_D = F_d - G_f*M_d;
N_D_d = pinv(G_f)*G_d;
G_D = G_d - G_f*N_D_d;
C_D = C;
% Simulation for actuator fault fault (f_u = -0.1) no sensor fault
x_0 = [0;0;0;0;0;0];            % Initial conditions
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;-0.35];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % enable/disable VA meachanism (0 = disable) 
f_m = [0;0;0];     % Sensor fault vector (added to [y1;y2;y3])
f_m_time = 8.5;                 % Sensor fault occurence time
%%
sim_done = sim('threeDiskOscillatorRig')
ymeas = sim_done.y_meas.Data;
theta_ref = sim_done.theta_ref.Data;

r = sim_done.r.Data;
time = sim_done.r.Time; % Adjust according to your time step if available

% Plot ymeas and theta_ref on the same plot
figure;

% Plot ymeas (3D) and theta_ref (1D) on the same figure
subplot(2, 1, 1); % Two plots: one for ymeas and theta_ref, the other for r
hold on;
plot(time, ymeas(:, 1), 'r', 'LineWidth', 3.5 ,'DisplayName', 'y1 meas ');
plot(time, ymeas(:, 2), 'g','LineWidth', 3.5 ,'DisplayName', 'y2 meas ');
plot(time, ymeas(:, 3), 'b', 'LineWidth', 3.5,'DisplayName', 'y3 meas ');
plot(time, theta_ref, 'k', 'LineWidth', 3.5, 'DisplayName', 'theta ref');
xline(detect_time, '--m', 'LineWidth', 3.5, 'DisplayName', 'detect\_time');  % vertical line
xlabel('Time');
ylabel('Values');
legend('show');
title('measured output y and theta_ref');

% Plot r in a separate plot

subplot(2, 1, 2);
hold on
plot(time, r(:, 1), 'b', 'LineWidth', 3.5, 'DisplayName', 'r1 ');
plot(time, r(:, 2), 'c', 'LineWidth', 3.5, 'DisplayName', 'r2 ');
plot(time, r(:, 3), 'r', 'LineWidth', 3.5, 'DisplayName', 'r3 ');
xline(detect_time, '--m', 'LineWidth', 3.5, 'DisplayName', 'detect\_time');  % vertical line
xlabel('Time');
ylabel('r');
legend('r1','r2','r3','detect\_time');
title('Plot of r');
set(gcf, 'Position', get(0, 'Screensize'))
saveas(gcf, fullfile('fig', 'Virtual_Actuator.png'));


