function [F_aero_b, M_aero_C_b] = aero_loads(X,U,i_cond)

global aircraft

[CD,CY,CL,Cl,Cm,Cn] = aero_databank(X,U,i_cond);

V = X(1); h = X(5); alpha_deg = X(2) ;beta_deg = X(7);

rho = ISA(h);
q_bar = 0.5*rho*V^2;

S = aircraft.S;
c = aircraft.c;
b = aircraft.b;

F_aero_a = q_bar*S*[-CD; -CY; -CL];

C_mbeta = DCM(3, -deg2rad(beta_deg));
C_alpha = DCM(2, deg2rad(alpha_deg));

C_ba = C_alpha*C_mbeta;

F_aero_b = C_ba*F_aero_a;

M_aero_C_b = q_bar*S*[b*Cl; c*Cm; b*Cn];
