function [T1, T2] = prop_loads(X, U)

global aircraft
run("create_aircraft.m");

h = X(5);
rho = ISA(h);

throttle_1 = U(1);
throttle_2 = U(2);
n_rho = aircraft.n_rho;
Tmax = aircraft.Tmax;

rho0 = 1.225;

T1 = throttle_1*Tmax*(rho/rho0)^n_rho;
T2 = throttle_2*Tmax*(rho/rho0)^n_rho;

% 1ª rotação: ângulo tau_p em torno do eixo z
C_tau_p_1 = [ cos(aircraft.tau_1_deg)  sin(aircraft.tau_1_deg)  0;
           -sin(aircraft.tau_1_deg)  cos(aircraft.tau_1_deg)  0;
            0           0           1 ];

C_tau_p_2 = [ cos(aircraft.tau_2_deg)  sin(aircraft.tau_2_deg)  0;
           -sin(aircraft.tau_2_deg)  cos(aircraft.tau_2_deg)  0;
            0           0           1 ];


% 2ª rotação: ângulo i_p em torno do eixo y do sistema intermediário
C_i_p_1 = [  cos(aircraft.iota_1_deg)   0   -sin(aircraft.iota_1_deg);
           0          1    0;
           sin(aircraft.iota_1_deg)   0    cos(aircraft.iota_1_deg) ];

C_i_p_2 = [  cos(aircraft.iota_2_deg)   0   -sin(aircraft.iota_2_deg);
           0          1    0;
           sin(aircraft.iota_2_deg)   0    cos(aircraft.iota_2_deg) ];


C_p_b_1 = C_i_p_1 * C_tau_p_1;
C_b_p_1 = C_p_b_1';

C_p_b_2 = C_i_p_2 * C_tau_p_2;
C_b_p_2 = C_p_b_2';

F_prop_p_1 = [T1 
    0 
    0];

F_prop_p_2 = [T2 
    0 
    0];

F_prop_b_1 = C_b_p_1 * F_prop_p_1;
F_prop_b_2 = C_b_p_2 * F_prop_p_2;
F_prop_b = F_prop_b_1 + F_prop_b_2;

x_1 = aircraft.x_1;
y_1 = aircraft.y_1;
z_1 = aircraft.z_1;

R_1 = [0 -z_1 y_1
    z_1 0 -x_1
    -y_1 x_1 0];

x_2 = aircraft.x_2;
y_2 = aircraft.y_2;
z_2 = aircraft.z_2;

R_2 = [0 -z_2 y_2
    z_2 0 -x_2
    -y_2 x_2 0];

M_prop_C_b_1 = R_1 * F_prop_b_1;
M_prop_C_b_2 = R_2 * F_prop_b_2;
M_prop_C_b = M_prop_C_b_1 + M_prop_C_b_2;

