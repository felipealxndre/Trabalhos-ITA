function [Xdot, Y] = dynamics(t, X, U)
% Dinâmica completa 6 GDL (eixos do corpo, ângulos de Euler em graus)
% Estados: X = [ V α q θ h x β ϕ p r ψ y ]^T
% Controles: U conforme seu modelo de dois motores + superfícies
% Saídas: Y = [ γ T1 T2 Mach CD CL Cm CY Cl Cn ρ q̄ ]^T

% -------------------------------------------------------------------------
% Leitura de estados (todos em unidades usadas no curso: deg, deg/s)
global g
global aircraft
run("create_aircraft.m");

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

% Conversões úteis
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
x_1    = aircraft.x_1;
x_2    = aircraft.x_2;
z_1    = aircraft.z_1;
z_2    = aircraft.z_2;
iota_1_deg = aircraft.iota_1_deg;
iota_2_deg = aircraft.iota_2_deg;

Gamma  = Ixx*Izz - Ixz^2;   % parâmetro inercial

% -------------------------------------------------------------------------
% Componentes de velocidade no corpo a partir de (V, α, β)
% u = V cosα cosβ ; v = V sinβ ; w = V sinα cosβ
u = V * cos(alpha_rad) * cos(beta_rad);
v = V * sin(beta_rad);
w = V * sin(alpha_rad) * cos(beta_rad);

% -------------------------------------------------------------------------
% Aerodinâmica e propulsão
% aero_loads retorna força lateral Y_for e momentos l_mom (rolamento), n_mom (guinada)
[D, L, M_mom, Y_for, l_mom, n_mom] = aero_loads(X, U);
[T1, T2] = prop_loads(X, U);

% -------------------------------------------------------------------------
% Dinâmica longitudinal (Cap. 6 – caso geral com dois motores)
Vdot = ( -D ...
         + T1*cosd(iota_1_deg + alpha_deg) ...
         + T2*cosd(iota_2_deg + alpha_deg) ...
         - m*g*sind(theta_deg - alpha_deg) ) / m;

alphadot_rad_s = q_rad_s + ( -L ...
         - T1*sind(iota_1_deg + alpha_deg) ...
         - T2*sind(iota_2_deg + alpha_deg) ...
         + m*g*cosd(theta_deg - alpha_deg) ) / (m*V);

% Momento de arfagem (q̇) com braços de empuxo (x_i, z_i)
qdot_rad_s2 = ( M_mom ...
         + z_1*T1*cosd(iota_1_deg) + z_2*T2*cosd(iota_2_deg) ...
         + x_1*T1*sind(iota_1_deg) + x_2*T2*sind(iota_2_deg) ) / Iyy;

% -------------------------------------------------------------------------
% Dinâmica lateral-direcional (Cap. 9 – equações gerais com acoplamentos)
% Se houver termos de momento propulsivo explícitos em x/z no CG, some-os aqui
Mprop_Cx = 0;
Mprop_Cz = 0;

pdot_rad_s2 = ( Izz/Gamma )*( l_mom + Mprop_Cx ) + ( Ixz/Gamma )*( n_mom + Mprop_Cz ) ...
              - ( 1/Gamma )*( (Izz^2 + Ixz^2 - Iyy*Izz)*q_rad_s*r_rad_s ...
                               + Ixz*(Iyy - Ixx - Izz)*p_rad_s*q_rad_s );

% qdot já calculado acima (forma simplificada em torno do eixo y)

rdot_rad_s2 = ( Ixz/Gamma )*( l_mom + Mprop_Cx ) + ( Ixx/Gamma )*( n_mom + Mprop_Cz ) ...
              - ( 1/Gamma )*( Ixz*(Ixx + Izz - Iyy)*q_rad_s*r_rad_s ...
                               + (Ixx*Iyy - Ixx^2 - Ixz^2)*p_rad_s*q_rad_s );

% -------------------------------------------------------------------------
% Cinemática de Euler (Cap. 9)
phidot_rad_s   = p_rad_s + tan(theta_rad)*( q_rad_s*sin(phi_rad) + r_rad_s*cos(phi_rad) );
thetadot_rad_s = q_rad_s*cos(phi_rad) - r_rad_s*sin(phi_rad);
psidot_rad_s   = ( q_rad_s*sin(phi_rad) + r_rad_s*cos(phi_rad) )/cos(theta_rad);

% -------------------------------------------------------------------------
% Cinemática de translação (Cap. 9) – derivadas em RI via Cϕ Cθ Cψ
xdot =  u*cos(theta_rad)*cos(psi_rad) ...
      + v*( sin(phi_rad)*sin(theta_rad)*cos(psi_rad) - cos(phi_rad)*sin(psi_rad) ) ...
      + w*( cos(phi_rad)*sin(theta_rad)*cos(psi_rad) + sin(phi_rad)*sin(psi_rad) );

ydot =  u*cos(theta_rad)*sin(psi_rad) ...
      + v*( cos(phi_rad)*cos(psi_rad) + sin(phi_rad)*sin(theta_rad)*sin(psi_rad) ) ...
      + w*( -sin(phi_rad)*cos(psi_rad) + cos(phi_rad)*sin(theta_rad)*sin(psi_rad) );

hdot =  u*sin(theta_rad) - v*sin(phi_rad)*cos(theta_rad) - w*cos(phi_rad)*cos(theta_rad);

% -------------------------------------------------------------------------
% β̇ (forma prática usando v̇ no corpo)
% vdot_b = (Y_for + m*g*sinϕ cosθ)/m - (r*u - p*w)
vdot_b = ( Y_for + m*g*( sin(phi_rad)*cos(theta_rad) ) )/m - ( r_rad_s*u - p_rad_s*w );
betadot_rad_s = ( V*vdot_b - v*Vdot ) / ( V*sqrt( u^2 + w^2 ) );

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

% Ângulo de trajetória (geral): γ = atan2(ḣ, sqrt(ẋ^2 + ẏ^2))
gamma_deg = rad2deg( atan2( hdot, hypot(xdot, ydot) ) );

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
