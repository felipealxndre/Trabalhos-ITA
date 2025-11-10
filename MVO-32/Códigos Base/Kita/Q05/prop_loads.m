function [F_prop_b, M_prop_C_b, T1, T2] = prop_loads(X,U)

% MVO-32 ESTABILIDADE E CONTROLE DE AERONAVES
% Dupla: Thiago Akira Missato e Eduardo Pinto Ferreira.
% Código adaptado para separar os dois motores.

global aircraft

h = X(5);
rho = ISA(h);
n_rho = aircraft.n_rho;
Tmax = aircraft.Tmax;

delta_1 = U(1); 
delta_2 = U(2);

T1 = delta_1 * Tmax * (rho/1.225)^n_rho;
T2 = delta_2 * Tmax * (rho/1.225)^n_rho;

% --- Motor 1 (Esquerdo) ---
iota_1_rad = aircraft.iota_1_deg * pi/180;
tau_1_rad = aircraft.tau_1_deg * pi/180;
r_1_b = [aircraft.x_1; aircraft.y_1; aircraft.z_1];

C_tau_1 = DCM(3, tau_1_rad);
C_iota_1 = DCM(2, iota_1_rad);
C_p1_b = C_iota_1 * C_tau_1;

F_prop_p1 = [T1; 0; 0]; 
F_prop_b1 = C_p1_b.' * F_prop_p1; 
M_prop_C_b1 = skew(r_1_b) * F_prop_b1;

% --- Motor 2 (Direito) ---
iota_2_rad = aircraft.iota_2_deg * pi/180;
tau_2_rad = aircraft.tau_2_deg * pi/180;
r_2_b = [aircraft.x_2; aircraft.y_2; aircraft.z_2];

C_tau_2 = DCM(3, tau_2_rad);
C_iota_2 = DCM(2, iota_2_rad);
C_p2_b = C_iota_2 * C_tau_2; 

F_prop_p2 = [T2; 0; 0]; 
F_prop_b2 = C_p2_b.' * F_prop_p2; 
M_prop_C_b2 = skew(r_2_b) * F_prop_b2;

% --- Total ---
F_prop_b = F_prop_b1 + F_prop_b2;
M_prop_C_b = M_prop_C_b1 + M_prop_C_b2;

end