function [CD,CY,CL,Cl,Cm,Cn] = aero_databank(X,U,i_cond)

global aircraft

CL_0 = 0.221;
CL_alpha = 6.00;
CL_q = 8.78;
CL_delta_e = 0.349;

Cm_0 = 0.127;
Cm_alpha = -1.678;
Cm_q = -52.0;
Cm_delta_e = -1.93;

CD_0 = 0.027;
CD_alpha = 0.025;
CD_delta_e = 0;

CY_beta = 1.195;
CY_p = -0.0214;
CY_r = -0.113;
CY_delta_r = -0.389;

Cl_beta = -0.0280;
Cl_p = -0.0215;
Cl_r = 0.0057;
Cl_delta_a = -0.0293;
Cl_delta_r = 0.0097;

Cn_beta = 0.0359;
Cn_p = -0.0027;
Cn_r = -0.0119;
Cn_delta_a = -0.0001;
Cn_delta_r = -0.0384;

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
beta_deg = X(7);
p_deg_s = X(9);
r_deg_s = X(10);

delta_e_deg = U(2);
delta_a_deg = U(3);
delta_r_deg = U(4);

alpha_rad = alpha_deg*pi/180;
beta_rad = beta_deg*pi/180;
p_rad_s = p_deg_s*pi/180;
q_rad_s = q_deg_s*pi/180;
r_rad_s = r_deg_s*pi/180;

delta_e_rad = delta_e_deg*pi/180;
delta_a_rad = delta_a_deg*pi/180;
delta_r_rad = delta_r_deg*pi/180;

c = aircraft.c;
b = aircraft.b;

CL = CL_0 + ...
    CL_alpha*alpha_rad + ...
    CL_q*(q_rad_s*c/(2*V)) + ...
    CL_delta_e*delta_e_rad;

CD = CD_0 + ...
    CD_alpha*alpha_rad + ...
    CD_delta_e*delta_e_rad; 

Cm = Cm_0 + ...
    Cm_alpha*alpha_rad + ...
    Cm_q*(q_rad_s*c/(2*V)) + ...
    Cm_delta_e*delta_e_rad;

CY = CY_beta*beta_rad + ...
    CY_p*(p_rad_s*b/(2*V)) + ...
    CY_r*(r_rad_s*b/(2*V)) + ...
    CY_delta_r*delta_r_rad;

Cl = Cl_beta*beta_rad + ...
    Cl_p*(p_rad_s*b/(2*V)) + ...
    Cl_r*(r_rad_s*b/(2*V)) + ...
    Cl_delta_a*delta_a_rad + ...
    Cl_delta_r*delta_r_rad;

Cn = Cn_beta*beta_rad + ...
    Cn_p*(p_rad_s*b/(2*V)) + ...
    Cn_r*(r_rad_s*b/(2*V)) + ...
    Cn_delta_a*delta_a_rad + ...
    Cn_delta_r*delta_r_rad;

end