function [F_prop_b, M_prop_C_b, T1, T2] = prop_loads(X,U)
% Modificado para Lista 3: lida com 2 motores independentes.

global aircraft

% --- Parâmetros comuns ---
h = X(5);
rho = ISA(h);
n_rho = aircraft.n_rho;
Tmax = aircraft.Tmax;

% --- Controles ---
delta_1 = U(1); % Controle motor esquerdo (1)
delta_2 = U(2); % Controle motor direito (2)

% --- Cálculo de Empuxo (Eq 1 e 2 do PDF) ---
T1 = delta_1 * Tmax * (rho/1.225)^n_rho;
T2 = delta_2 * Tmax * (rho/1.225)^n_rho;

% --- Motor 1 (Esquerdo) ---
iota_1_rad = aircraft.iota_1_deg * pi/180;
tau_1_rad = aircraft.tau_1_deg * pi/180;
r_1_b = [aircraft.x_1; aircraft.y_1; aircraft.z_1];

% Matrizes de transformação (Slide Cap 4 / PDF Tabela 2)
C_tau_1 = DCM(3, tau_1_rad);
C_iota_1 = DCM(2, iota_1_rad);
C_p1_b = C_iota_1 * C_tau_1; % C_p/b = C_iota * C_tau

% Força no sistema do corpo (Slide Cap 8)
F_prop_p1 = [T1; 0; 0]; % Força no sistema do motor
F_prop_b1 = C_p1_b.' * F_prop_p1; % Força no sistema do corpo

% Momento no sistema do corpo (Slide Cap 8)
M_prop_C_b1 = skew(r_1_b) * F_prop_b1;

% --- Motor 2 (Direito) ---
iota_2_rad = aircraft.iota_2_deg * pi/180;
tau_2_rad = aircraft.tau_2_deg * pi/180;
r_2_b = [aircraft.x_2; aircraft.y_2; aircraft.z_2];

% Matrizes de transformação
C_tau_2 = DCM(3, tau_2_rad);
C_iota_2 = DCM(2, iota_2_rad);
C_p2_b = C_iota_2 * C_tau_2; % C_p/b = C_iota * C_tau

% Força no sistema do corpo
F_prop_p2 = [T2; 0; 0]; % Força no sistema do motor
F_prop_b2 = C_p2_b.' * F_prop_p2; % Força no sistema do corpo

% Momento no sistema do corpo
M_prop_C_b2 = skew(r_2_b) * F_prop_b2;

% --- Totais ---
F_prop_b = F_prop_b1 + F_prop_b2;
M_prop_C_b = M_prop_C_b1 + M_prop_C_b2;

end