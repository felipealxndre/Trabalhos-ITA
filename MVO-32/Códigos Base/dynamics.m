function [Xdot, Y] = dynamics(t, X, U)
% Dinâmica completa 6 GDL (eixos do corpo, ângulos de Euler em graus)
% Implementa as equações de 6-DOF nos eixos do corpo (Cap. 8)

% -------------------------------------------------------------------------
% Leitura de estados (todos em unidades usadas no curso: deg, deg/s)
global g
global aircraft

V          = X(1);
alpha_deg  = X(2);
q_deg_s    = X(3);
theta_deg  = X(4);
h          = X(5);
% x        = X(6);
beta_deg   = X(7);
phi_deg    = X(8);
p_deg_s    = X(9);
r_deg_s    = X(10);
psi_deg    = X(11);
% y        = X(12); 

% Conversões úteis para radianos
alpha_rad  = deg2rad(alpha_deg);
beta_rad   = deg2rad(beta_deg);
phi_rad    = deg2rad(phi_deg);
theta_rad  = deg2rad(theta_deg);
psi_rad    = deg2rad(psi_deg);

p_rad_s    = deg2rad(p_deg_s);
q_rad_s    = deg2rad(q_deg_s);
r_rad_s    = deg2rad(r_deg_s);

% -------------------------------------------------------------------------
% Parâmetros da aeronave
m      = aircraft.m;
Ixx    = aircraft.Ixx;
Iyy    = aircraft.Iyy;
Izz    = aircraft.Izz;
Ixz    = aircraft.Ixz;
Gamma  = Ixx*Izz - Ixz^2;   % parâmetro inercial

% -------------------------------------------------------------------------
% Componentes de velocidade no corpo a partir de (V, α, β)
u = V * cos(alpha_rad) * cos(beta_rad);
v = V * sin(beta_rad);
w = V * sin(alpha_rad) * cos(beta_rad);

% -------------------------------------------------------------------------
% Obter Forças e Momentos
% 1. Aerodinâmica
[D, L, m, Y, l, n] = aero_loads(X, U);
M_aero_b = [l; m; n];

% 2. Propulsão (Já em eixos do corpo)
[F_prop_b, M_prop_C_b, T1, T2] = prop_loads(X, U);

% 3. Gravidade (Em eixos do corpo)
F_grav_b = [-m*g*sin(theta_rad);
             m*g*cos(theta_rad)*sin(phi_rad);
             m*g*cos(theta_rad)*cos(phi_rad)];

% 4. Converter Aero (D, L) para eixos do corpo
C_b_a = [cos(alpha_rad)*cos(beta_rad), -cos(alpha_rad)*sin(beta_rad), -sin(alpha_rad);
         sin(beta_rad),                  cos(beta_rad),                   0;
         sin(alpha_rad)*cos(beta_rad), -sin(alpha_rad)*sin(beta_rad),  cos(alpha_rad)];
      
F_aero_wind = [-D; Y; -L]; 
F_aero_b = C_b_a * F_aero_wind;

% -------------------------------------------------------------------------
% Dinâmica de Translação (Eixos do Corpo) 
F_total_b = F_aero_b + F_prop_b + F_grav_b;

udot = F_total_b(1)/m - q_rad_s*w + r_rad_s*v;
vdot = F_total_b(2)/m - r_rad_s*u + p_rad_s*w;
wdot = F_total_b(3)/m - p_rad_s*v + q_rad_s*u;

% Converter de (udot, vdot, wdot) para (Vdot, alphadot, betadot) 
Vdot = (u*udot + v*vdot + w*wdot) / V;
alphadot_rad_s = (u*wdot - w*udot) / (u^2 + w^2);
betadot_rad_s  = (vdot*V - v*Vdot) / (V * sqrt(u^2 + w^2));

% -------------------------------------------------------------------------
% Dinâmica de Rotação (Eixos do Corpo) 
M_total_b = M_aero_b + M_prop_C_b;

% LINHA 77 CORRIGIDA (sinal - no termo qr):
pdot_rad_s2 = ( Izz*M_total_b(1) + Ixz*M_total_b(3) - (Ixz*(Iyy - Ixx - Izz)*p_rad_s*q_rad_s) - (Ixz^2 + Izz*(Izz - Iyy))*q_rad_s*r_rad_s ) / Gamma;
qdot_rad_s2 = ( M_total_b(2) - (Ixx - Izz)*p_rad_s*r_rad_s - Ixz*(p_rad_s^2 - r_rad_s^2) ) / Iyy;
rdot_rad_s2 = ( Ixz*M_total_b(1) + Ixx*M_total_b(3) - (Ixz*(Izz + Ixx - Iyy)*q_rad_s*r_rad_s) + (Ixz^2 + Ixx*(Ixx - Iyy))*p_rad_s*q_rad_s ) / Gamma;

% -------------------------------------------------------------------------
% Cinemática de Euler
phidot_rad_s   = p_rad_s + tan(theta_rad)*( q_rad_s*sin(phi_rad) + r_rad_s*cos(phi_rad) );
thetadot_rad_s = q_rad_s*cos(phi_rad) - r_rad_s*sin(phi_rad);
psidot_rad_s   = ( q_rad_s*sin(phi_rad) + r_rad_s*cos(phi_rad) )/cos(theta_rad);

% -------------------------------------------------------------------------
% Cinemática de translação
xdot =  u*cos(theta_rad)*cos(psi_rad) ...
      + v*( sin(phi_rad)*sin(theta_rad)*cos(psi_rad) - cos(phi_rad)*sin(psi_rad) ) ...
      + w*( cos(phi_rad)*sin(theta_rad)*cos(psi_rad) + sin(phi_rad)*sin(psi_rad) );
      
ydot =  u*cos(theta_rad)*sin(psi_rad) ...
      + v*( cos(phi_rad)*cos(psi_rad) + sin(phi_rad)*sin(theta_rad)*sin(psi_rad) ) ...
      + w*( -sin(phi_rad)*cos(psi_rad) + cos(phi_rad)*sin(theta_rad)*sin(psi_rad) );
      
hdot =  ( u*sin(theta_rad) - v*sin(phi_rad)*cos(theta_rad) - w*cos(phi_rad)*cos(theta_rad) );

% -------------------------------------------------------------------------
% Conversões finais para deg/deg/s (padrão do seu Xdot)
alphadot_deg_s = rad2deg(alphadot_rad_s);
qdot_deg_s2    = rad2deg(qdot_rad_s2);
pdot_deg_s2    = rad2deg(pdot_rad_s2);
rdot_deg_s2    = rad2deg(rdot_rad_s2);
phidot_deg_s   = rad2deg(phidot_rad_s);
thetadot_deg_s = rad2deg(thetadot_rad_s);
psidot_deg_s   = rad2deg(psidot_rad_s);
betadot_deg_s  = rad2deg(betadot_rad_s);

% -------------------------------------------------------------------------
% Empacota Xdot
Xdot = [
    Vdot
    alphadot_deg_s
    qdot_deg_s2
    thetadot_deg_s
    hdot
    xdot
    betadot_deg_s
    phidot_deg_s
    pdot_deg_s2
    rdot_deg_s2
    psidot_deg_s
    ydot
];

% -------------------------------------------------------------------------
% Saídas Y
[rho, ~, ~, a] = ISA(h);
Mach  = V / a;
q_bar = 0.5 * rho * V^2;

[CD, CL, Cm, CY, Cl, Cn] = aero_databank(X, U);

gamma_deg = rad2deg( asin( hdot / V ) );

Y = [
    gamma_deg
    T1
    T2
    Mach
    CD
    CL
    Cm
    CY
    Cl
    Cn
    rho
    q_bar
];
end