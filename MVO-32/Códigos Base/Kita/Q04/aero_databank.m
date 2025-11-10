function [CD,CY,CL,Cl,Cm,Cn] = aero_databank(X,U,i_cond)

global aircraft

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
beta_deg = X(7);
p_deg_s = X(9);
r_deg_s = X(10);

i_t_deg = U(3);
delta_e_deg = U(4);
delta_a_deg = U(5);
delta_r_deg = U(6);

alpha_rad = alpha_deg*pi/180;
beta_rad = beta_deg*pi/180;

p_rad_s = p_deg_s*pi/180;
q_rad_s = q_deg_s*pi/180;
r_rad_s = r_deg_s*pi/180;

i_t_rad = i_t_deg*pi/180;
delta_e_rad = delta_e_deg*pi/180;
delta_a_rad = delta_a_deg*pi/180;
delta_r_rad = delta_r_deg*pi/180;

c = aircraft.c;
b = aircraft.b;

switch i_cond
    case 1
        % Cruise Flight Condition - Lista 03
        CL_0 = aircraft.CL_0;
        CL_alpha = aircraft.CL_alpha;
        CL_q = aircraft.CL_q;
        CL_it = aircraft.CL_i_t;
        CL_delta_e = aircraft.CL_delta_e; 
        CD_0 = aircraft.CD_0;
        CD_alpha = aircraft.CD_alpha;
        CD_alpha2 = aircraft.CD_alpha2;
        CD_q = aircraft.CD_q;
        CD_it = aircraft.CD_i_t_0 + aircraft.CD_i_t_alpha * alpha_rad;
        CD_delta_e = aircraft.CD_delta_e_0 + aircraft.CD_delta_e_alpha * alpha_rad;
        Cm_0 = aircraft.Cm_0;
        Cm_alpha = aircraft.Cm_alpha;
        Cm_q = aircraft.Cm_q;
        Cm_it = aircraft.Cm_i_t;
        Cm_delta_e = aircraft.Cm_delta_e;
        CY_beta = aircraft.CY_beta;
        CY_p = aircraft.CY_p_0 + aircraft.CY_p_alpha * alpha_rad;
        CY_r = aircraft.CY_r;
        CY_delta_r = aircraft.CY_delta_r;
        Cl_beta = aircraft.Cl_beta_0 + aircraft.Cl_beta_alpha * alpha_rad;
        Cl_p = aircraft.Cl_p;
        Cl_r = aircraft.Cl_r_0 + aircraft.Cl_r_alpha * alpha_rad;
        Cl_delta_a = aircraft.Cl_delta_a;
        Cl_delta_r = aircraft.Cl_delta_r;
        Cn_beta = aircraft.Cn_beta;
        Cn_p = aircraft.Cn_p_0 + aircraft.Cn_p_alpha * alpha_rad;
        Cn_r = aircraft.Cn_r;
        Cn_delta_a = aircraft.Cn_delta_a_0 + aircraft.Cn_delta_a_alpha * alpha_rad;
        Cn_delta_r = aircraft.Cn_delta_r;
end

CL = CL_0 + ...
    CL_alpha*alpha_rad + ...
    CL_q*(q_rad_s*c/(2*V))+ ...
    CL_it*i_t_rad+ ...
    CL_delta_e*delta_e_rad;

CD = CD_0 + ...
    CD_alpha*alpha_rad + ...
    CD_alpha2*alpha_rad^2 + ...
    CD_q*(q_rad_s*c/(2*V))+ ...
    CD_it*i_t_rad+ ...
    CD_delta_e*delta_e_rad;

Cm = Cm_0 + ...
    Cm_alpha*alpha_rad + ...
    Cm_q*(q_rad_s*c/(2*V))+ ...
    Cm_it*i_t_rad+ ...
    Cm_delta_e*delta_e_rad;

CY = CY_beta*beta_rad+ ...
    CY_p*(p_rad_s*b/(2*V))+ ...
    CY_r*(r_rad_s*b/(2*V))+ ...
    CY_delta_r*delta_r_rad;

Cl = Cl_beta*beta_rad+ ...
    Cl_p*(p_rad_s*b/(2*V))+ ...
    Cl_r*(r_rad_s*b/(2*V))+ ...
    Cl_delta_a*delta_a_rad+ ...
    Cl_delta_r*delta_r_rad;

Cn = Cn_beta*beta_rad+ ...
    Cn_p*(p_rad_s*b/(2*V))+ ...
    Cn_r*(r_rad_s*b/(2*V))+ ...
    Cn_delta_a*delta_a_rad+ ...
    Cn_delta_r*delta_r_rad;
