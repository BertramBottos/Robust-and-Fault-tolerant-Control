function   InputSyms()
%INPUTSYMS Summary of this function goes here
syms y_1 y_2 y_3 u_1 u_2 s b_1 b_2 J_1 J_2 J_3 k_1 k_2 b_3 s 
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


end

