function [CD, CL, Cm, CY, Cl, Cn] = aero_databank(X, U)

global aircraft

d2r = pi/180;

V         = X(1);
alpha_deg = X(2);
q_deg_s   = X(3);
beta_deg  = X(7);
p_deg_s   = X(9);
r_deg_s   = X(10);

delta_e_deg = U(3);
delta_a_deg = U(4);
delta_r_deg = U(5);

alpha = alpha_deg * d2r;
beta  = beta_deg  * d2r;
p     = p_deg_s   * d2r;
q     = q_deg_s   * d2r;
r     = r_deg_s   * d2r;

delta_e = delta_e_deg * d2r;
delta_a = delta_a_deg * d2r;
delta_r = delta_r_deg * d2r;

S     = aircraft.S;
b     = aircraft.b;
c_bar = aircraft.c;
i_t   = aircraft.i_t; 

if V == 0
    V_norm_p = 0;
    V_norm_q = 0;
    V_norm_r = 0;
else
    V_norm_p = (p * b) / (2 * V);
    V_norm_q = (q * c_bar) / (2 * V);
    V_norm_r = (r * b) / (2 * V);
end

CL = aircraft.CL_0 + aircraft.CL_alpha*alpha + aircraft.CL_q*V_norm_q + aircraft.CL_i_t*i_t + aircraft.CL_delta_e*delta_e;

Cm = aircraft.Cm_0 + aircraft.Cm_alpha*alpha + aircraft.Cm_q*V_norm_q + aircraft.Cm_i_t*i_t + aircraft.Cm_delta_e*delta_e;

CDit = aircraft.CD_i_t_0 + aircraft.CD_i_t_alpha*alpha;
CDde = aircraft.CD_delta_e_0 + aircraft.CD_delta_e_alpha*alpha;

CD = aircraft.CD_0 + aircraft.CD_alpha*alpha + aircraft.CD_alpha2*alpha^2 + ...
     aircraft.CD_q*V_norm_q + CDit*i_t + CDde*delta_e + ...
     aircraft.CD_beta2*beta^2 + aircraft.CD_p2*V_norm_p^2 + aircraft.CD_r2*V_norm_r^2 + ...
     aircraft.CD_delta_a2*delta_a^2 + aircraft.CD_delta_r2*delta_r^2;

CYp = aircraft.CY_p_0 + aircraft.CY_p_alpha*alpha;

CY = aircraft.CY_beta*beta + CYp*V_norm_p + aircraft.CY_r*V_norm_r + aircraft.CY_delta_a*delta_a + aircraft.CY_delta_r*delta_r;

Clbeta = aircraft.Cl_beta_0 + aircraft.Cl_beta_alpha*alpha;
Clr = aircraft.Cl_r_0 + aircraft.Cl_r_alpha*alpha;

Cl = Clbeta*beta + aircraft.Cl_p*V_norm_p + Clr*V_norm_r + aircraft.Cl_delta_a*delta_a + aircraft.Cl_delta_r*delta_r;

Cnp = aircraft.Cn_p_0 + aircraft.Cn_p_alpha*alpha;
Cnda = aircraft.Cn_delta_a_0 + aircraft.Cn_delta_a_alpha*alpha;

Cn = aircraft.Cn_beta*beta + Cnp*V_norm_p + aircraft.Cn_r*V_norm_r + Cnda*delta_a + aircraft.Cn_delta_r*delta_r;

end