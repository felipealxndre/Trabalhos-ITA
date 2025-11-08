function [CD,CL,Cm,CY,Cl,Cn] = aero_databank(X,U)

global aircraft
run("create_aircraft.m");

CL_0 = aircraft.CL_0;
CL_alpha = aircraft.CL_alpha;
CL_q = aircraft.CL_q;
CL_it = aircraft.CL_it;
CL_delta_e = aircraft.CL_delta_e; 

Cm_0 = aircraft.Cm_0;
Cm_alpha = aircraft.Cm_alpha;
Cm_q = aircraft.Cm_q;
Cm_it = aircraft.Cm_it;
Cm_delta_e = aircraft.Cm_delta_e;

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

CY_beta = aircraft.CY_beta;
CY_p_0 = aircraft.CY_p_0;
CY_p_alpha = aircraft.CY_p_alpha;
CY_r = aircraft.CY_r;
CY_delta_a = aircraft.CY_delta_a;
CY_delta_r = aircraft.CY_delta_r;

Cl_beta_0 = aircraft.Cl_beta_0;
Cl_beta_alpha = aircraft.Cl_beta_alpha;
Cl_p = aircraft.Cl_p;
Cl_r_0 = aircraft.Cl_r_0;
Cl_r_alpha = aircraft.Cl_r_alpha;
Cl_delta_a = aircraft.Cl_delta_a;
Cl_delta_r = aircraft.Cl_delta_r;

Cn_beta = aircraft.Cn_beta;
Cn_p_0 = aircraft.Cn_p_0;
Cn_p_alpha = aircraft.Cn_p_alpha;
Cn_r = aircraft.Cn_r;
Cn_delta_a_0 = aircraft.Cn_delta_a_0;
Cn_delta_a_alpha = aircraft.Cn_delta_a_alpha;
Cn_delta_r = aircraft.Cn_delta_r;

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);

i_t_deg = U(3);
delta_e_deg = U(4);

alpha_rad = alpha_deg*pi/180;
q_rad_s = q_deg_s*pi/180;

i_t_rad = i_t_deg*pi/180;
delta_e_rad = delta_e_deg*pi/180;

c = aircraft.c;

CL = CL_0 + ...
    CL_alpha*alpha_rad + ...
    CL_q*(q_rad_s*c/(2*V))+ ...
    CL_it*i_t_rad+ ...
    CL_delta_e*delta_e_rad;


CD_it = CD_i_t_0 + CD_i_t_alpha*alpha_rad;
CD_delta_e = CD_delta_e_0 + CD_delta_e_alpha*alpha_rad;

CD = CD_0 + ...
    CD_alpha*alpha_rad + ...
    CD_alpha2*alpha_rad^2 + ...
    CD_q*q_adim + ...
    CD_it*i_t_rad + ...
    CD_delta_e*delta_e_rad + ...
    CD_beta2*beta_rad^2 + ...
    CD_p2*p_adim^2 + ...
    CD_r2*r_adim^2 + ...
    CD_delta_a2*delta_a_rad^2 + ...
    CD_delta_r2*delta_r_rad^2;



Cm = Cm_0 + ...
    Cm_alpha*alpha_rad + ...
    Cm_q*(q_rad_s*c/(2*V))+ ...
    Cm_it*i_t_rad+ ...
    Cm_delta_e*delta_e_rad;

% CY (Força Lateral)
CY_p = CY_p_0 + CY_p_alpha*alpha_rad;

CY = CY_beta*beta_rad + ...
    CY_p*p_adim + ...
    CY_r*r_adim + ...
    CY_delta_a*delta_a_rad + ...
    CY_delta_r*delta_r_rad;

% Cl (Momento de Rolamento)
Cl_beta = Cl_beta_0 + Cl_beta_alpha*alpha_rad;
Cl_r = Cl_r_0 + Cl_r_alpha*alpha_rad;

Cl = Cl_beta*beta_rad + ...
    Cl_p*p_adim + ...
    Cl_r*r_adim + ...
    Cl_delta_a*delta_a_rad + ...
    Cl_delta_r*delta_r_rad;

% Cn (Momento de Guinada)
Cn_p = Cn_p_0 + Cn_p_alpha*alpha_rad;
Cn_delta_a = Cn_delta_a_0 + Cn_delta_a_alpha*alpha_rad;

Cn = Cn_beta*beta_rad + ...
    Cn_p*p_adim + ...
    Cn_r*r_adim + ...
    Cn_delta_a*delta_a_rad + ...
    Cn_delta_r*delta_r_rad;