function [Xdot, Y]= dynamics(t, X, U)

global g
global aircraft

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
theta_deg = X(4);
h = X(5);
% x = X(6);


q_rad_s = deg2rad(q_deg_s);

m = aircraft.m;
iota_p_deg = aircraft.i_p_deg;
Iyy = aircraft.Iyy;
x_p = aircraft.x_p;
z_p = aircraft.z_p;

[D, L, M] = aero_loads(X, U);
[T1, T2] = prop_loads(X, U);

%Capitulo 6, slide 35:
Vdot = 1/m*(-D + T*cosd(iota_p_deg + alpha_deg)...
    - m*g*sind(theta_deg - alpha_deg));
alphadot_rad_s = q_rad_s...
    + 1/(m*V)*(-L -T*sind(iota_p_deg + alpha_deg) + m*g*cosd(theta_deg - alpha_deg));
qdot_rad_s2 = 1/Iyy*(M + z_p*T*cosd(iota_p_deg) + x_p*T*sind(iota_p_deg));
thetadot_deg_s = q_deg_s;
hdot = V*sind(theta_deg - alpha_deg);
xdot = V*cosd(theta_deg - alpha_deg);

alphadot_deg_s = rad2deg(alphadot_rad_s);
qdot_deg_s2 = rad2deg(qdot_rad_s2);

Xdot = [
    Vdot
    alphadot_deg_s
    qdot_deg_s2
    thetadot_deg_s
    hdot
    xdot
    ];

gamma_deg = theta_deg - alpha_deg;
[rho,~,~,a] = ISA(h);
Mach = V/a;
q_bar = 0.5*rho*V^2;
[CD,CL,Cm,CY,Cl,Cn] = aero_databank(X,U);

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