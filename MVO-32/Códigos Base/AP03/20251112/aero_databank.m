function [CD,CY,CL,Cl,Cm,Cn] = aero_databank(X,U,i_cond)

global aircraft

switch i_cond
    case 1
        % A1 flight condition:
        CL_0 = 0.738;
        CL_alpha = 5.66;
        CL_q = 6.40;
        CL_it = 0.779;
        CL_delta_e = 0.433; 
        CD_0 = 0.0958;
        CD_alpha = 0.3691;
        CD_alpha2 = 1.6257;
        CD_q = 0;
        CD_it = 0.1121;
        CD_delta_e = 0.0623;
        Cm_0 = -0.2422;
        Cm_alpha = -1.203;
        Cm_q = -27.22;
        Cm_it = -3.038;
        Cm_delta_e = -1.688;
        CY_beta = 1.090;
        CY_p = -0.692;
        CY_r = -0.657;
        CY_delta_r = -0.253;
        Cl_beta = -0.211;
        Cl_p = -0.394;
        Cl_r = 0.408;
        Cl_delta_a = -0.0387;
        Cl_delta_r = 0.0139;
        Cn_beta = 0.157;
        Cn_p = -0.276;
        Cn_r = -0.335;
        Cn_delta_a = -0.0075;
        Cn_delta_r = -0.141;
    case 2
        % A2 flight condition:
        CL_0 = 0.329;
        CL_alpha = 4.72;
        CL_q = 7.00;
        CL_it = 0.711;
        CL_delta_e = 0.395; 
        CD_0 = 0.0238;
        CD_alpha = 0.1418;
        CD_alpha2 = 1.1477;
        CD_q = 0;
        CD_it = 0.0455;
        CD_delta_e = 0.0253;
        Cm_0 = -0.1080;
        Cm_alpha = -0.747;
        Cm_q = -27.06;
        Cm_it = -2.774;
        Cm_delta_e = -1.541;
        CY_beta = 1.034;
        CY_p = -0.184;
        CY_r = -0.657;
        CY_delta_r = -0.230;
        Cl_beta = -0.184;
        Cl_p = -0.366;
        Cl_r = 0.179;
        Cl_delta_a = -0.0343;
        Cl_delta_r = 0.0206;
        Cn_beta = 0.163;
        Cn_p = -0.125;
        Cn_r = -0.276;
        Cn_delta_a = -0.0072;
        Cn_delta_r = -0.130;
    case 3
        % A3 flight condition:
        CL_0 = 0.347;
        CL_alpha = 6.22;
        CL_q = 7.60;
        CL_it = 0.349;
        CL_delta_e = 0.194; 
        CD_0 = 0.0232;
        CD_alpha = 0.2190;
        CD_alpha2 = 1.9973;
        CD_q = 0;
        CD_it = 0.0122;
        CD_delta_e = 0.0068;
        Cm_0 = -0.0334;
        Cm_alpha = -1.081;
        Cm_q = -35.44;
        Cm_it = -1.388;
        Cm_delta_e = -0.771;
        CY_beta = 1.149;
        CY_p = 0.0866;
        CY_r = -0.687;
        CY_delta_r = -0.176;
        Cl_beta = -0.235;
        Cl_p = -0.426;
        Cl_r = 0.148;
        Cl_delta_a = -0.0184;
        Cl_delta_r = 0.0193;
        Cn_beta = 0.238;
        Cn_p = -0.0981;
        Cn_r = -0.294;
        Cn_delta_a = -0.0084;
        Cn_delta_r = -0.101;
end

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
beta_deg = X(7);
p_deg_s = X(9);
r_deg_s = X(10);

i_t_deg = U(2);
delta_e_deg = U(3);
delta_a_deg = U(4);
delta_r_deg = U(5);

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
