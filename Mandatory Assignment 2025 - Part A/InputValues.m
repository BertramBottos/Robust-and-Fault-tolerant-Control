function  InputValues()
%UNTITLED Summary of this function goes here
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

H_ry = [k_1/J_2 (-k_1-k_2-b_2*s)/J_2-s^2 k_2/J_2; ...
   0 k_2/J_3 (-k_2-b_3*s)/J_3-s^2]

end

