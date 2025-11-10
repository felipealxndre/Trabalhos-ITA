function [CD,CY,CL,Cl,Cm,Cn] = aero_databank(X,U,i_cond)
% Modificado para Lista 3: 
% 1. Aceita i_cond (para compatibilidade com aero_loads)
% 2. Lê coeficientes da 'global aircraft'
% 3. Implementa equações completas (7-19) do PDF
% 4. Usa novos índices de controle U (6 elementos)

% NOTA: O 'i_cond' é recebido, mas não é usado neste databank,
% pois os dados são lidos diretamente da 'global aircraft'.
% Mantido para compatibilidade com a cadeia de chamadas.

global aircraft

% --- Carrega todos os coeficientes do 'aircraft' ---
% Geometria
c = aircraft.c;
b = aircraft.b;

% Coeficiente de Sustentação (CL)
CL_0 = aircraft.CL_0;
CL_alpha = aircraft.CL_alpha;
CL_q = aircraft.CL_q;
CL_i_t = aircraft.CL_i_t;
CL_delta_e = aircraft.CL_delta_e;

% Coeficiente de Momento de Arfagem (Cm)
Cm_0 = aircraft.Cm_0;
Cm_alpha = aircraft.Cm_alpha;
Cm_q = aircraft.Cm_q;
Cm_i_t = aircraft.Cm_i_t;
Cm_delta_e = aircraft.Cm_delta_e;

% Coeficiente de Arrasto (CD)
CD_0 = aircraft.CD_0;
CD_alpha = aircraft.CD_alpha;
CD_alpha2 = aircraft.CD_alpha2;
CD_q = aircraft.CD_q;
CD_i_t_0 = aircraft.CD_i_t_0;
CD_i_t_alpha = aircraft.CD_i_t_alpha;
CD_delta_e_0 = aircraft.CD_delta_e_0;
CD_delta_e_alpha = aircraft.CD_delta_e_alpha;
CD_beta2 = aircraft.CD_beta2;
CD_p2 = aircraft.CD_p2;
CD_r2 = aircraft.CD_r2;
CD_delta_a2 = aircraft.CD_delta_a2;
CD_delta_r2 = aircraft.CD_delta_r2;

% Coeficiente de Força Lateral (CY)
CY_beta = aircraft.CY_beta;
CY_p_0 = aircraft.CY_p_0;
CY_p_alpha = aircraft.CY_p_alpha;
CY_r = aircraft.CY_r;
CY_delta_a = aircraft.CY_delta_a;
CY_delta_r = aircraft.CY_delta_r;

% Coeficiente de Momento de Rolamento (Cl)
Cl_beta_0 = aircraft.Cl_beta_0;
Cl_beta_alpha = aircraft.Cl_beta_alpha;
Cl_p = aircraft.Cl_p;
Cl_r_0 = aircraft.Cl_r_0;
Cl_r_alpha = aircraft.Cl_r_alpha;
Cl_delta_a = aircraft.Cl_delta_a;
Cl_delta_r = aircraft.Cl_delta_r;

% Coeficiente de Momento de Guinada (Cn)
Cn_beta = aircraft.Cn_beta;
Cn_p_0 = aircraft.Cn_p_0;
Cn_p_alpha = aircraft.Cn_p_alpha;
Cn_r = aircraft.Cn_r;
Cn_delta_a_0 = aircraft.Cn_delta_a_0;
Cn_delta_a_alpha = aircraft.Cn_delta_a_alpha;
Cn_delta_r = aircraft.Cn_delta_r;

% --- Extrai variáveis de estado (X) ---
V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
beta_deg = X(7);
p_deg_s = X(9);
r_deg_s = X(10);

% --- Extrai controles (U) ---
% U = [delta_1, delta_2, i_t, delta_e, delta_a, delta_r]
i_t_deg = U(3);
delta_e_deg = U(4);
delta_a_deg = U(5);
delta_r_deg = U(6);

% --- Conversões para Radianos ---
alpha_rad = alpha_deg*pi/180;
beta_rad = beta_deg*pi/180;

p_rad_s = p_deg_s*pi/180;
q_rad_s = q_deg_s*pi/180;
r_rad_s = r_deg_s*pi/180;

i_t_rad = i_t_deg*pi/180;
delta_e_rad = delta_e_deg*pi/180;
delta_a_rad = delta_a_deg*pi/180;
delta_r_rad = delta_r_deg*pi/180;

% --- Adimensionalização das velocidades angulares ---
p_adim = p_rad_s * b / (2*V);
q_adim = q_rad_s * c / (2*V);
r_adim = r_rad_s * b / (2*V);

% --- Implementação das Equações (7-19) do PDF ---

% Coeficientes dependentes de alpha
CD_i_t = CD_i_t_0 + CD_i_t_alpha * alpha_rad; % Eq (10)
CD_delta_e = CD_delta_e_0 + CD_delta_e_alpha * alpha_rad; % Eq (11)
CY_p = CY_p_0 + CY_p_alpha * alpha_rad; % Eq (13)
Cl_beta = Cl_beta_0 + Cl_beta_alpha * alpha_rad; % Eq (15)
Cl_r = Cl_r_0 + Cl_r_alpha * alpha_rad; % Eq (16)
Cn_p = Cn_p_0 + Cn_p_alpha * alpha_rad; % Eq (18)
Cn_delta_a = Cn_delta_a_0 + Cn_delta_a_alpha * alpha_rad; % Eq (19)

% Equações principais
CL = CL_0 + ...
    CL_alpha * alpha_rad + ...
    CL_q * q_adim + ...
    CL_i_t * i_t_rad + ...
    CL_delta_e * delta_e_rad; % Eq (7)

Cm = Cm_0 + ...
    Cm_alpha * alpha_rad + ...
    Cm_q * q_adim + ...
    Cm_i_t * i_t_rad + ...
    Cm_delta_e * delta_e_rad; % Eq (8)

CD = CD_0 + ...
    CD_alpha * alpha_rad + ...
    CD_alpha2 * alpha_rad^2 + ...
    CD_q * q_adim + ...
    CD_i_t * i_t_rad + ...
    CD_delta_e * delta_e_rad + ...
    CD_beta2 * beta_rad^2 + ...
    CD_p2 * p_adim^2 + ...
    CD_r2 * r_adim^2 + ...
    CD_delta_a2 * delta_a_rad^2 + ...
    CD_delta_r2 * delta_r_rad^2; % Eq (9)

CY = CY_beta * beta_rad + ...
    CY_p * p_adim + ...
    CY_r * r_adim + ...
    CY_delta_a * delta_a_rad + ...
    CY_delta_r * delta_r_rad; % Eq (12)

Cl = Cl_beta * beta_rad + ...
    Cl_p * p_adim + ...
    Cl_r * r_adim + ...
    Cl_delta_a * delta_a_rad + ...
    Cl_delta_r * delta_r_rad; % Eq (14)

Cn = Cn_beta * beta_rad + ...
    Cn_p * p_adim + ...
    Cn_r * r_adim + ...
    Cn_delta_a * delta_a_rad + ...
    Cn_delta_r * delta_r_rad; % Eq (17)

end