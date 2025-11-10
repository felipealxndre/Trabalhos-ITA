function [Xdot,Y] = dynamics(t,X,U,i_cond)

global g
global aircraft

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
theta_deg = X(4);
h = X(5);
x = X(6);

beta_deg = X(7);
phi_deg = X(8);
p_deg_s = X(9);
r_deg_s = X(10);
psi_deg = X(11);
y = X(12);

throttle_1 = U(1);
throttle_2 = U(2);
i_t_deg = U(3);
delta_e_deg = U(4);
delta_a_deg = U(5);
delta_r_deg = U(6);

m = aircraft.m;

Ixx = aircraft.Ixx;
Iyy = aircraft.Iyy;
Izz = aircraft.Izz;
Ixz = aircraft.Ixz;


C_psi = DCM(3, deg2rad(psi_deg));
C_theta = DCM(2, deg2rad(theta_deg));
C_phi = DCM(1, deg2rad(phi_deg));

C_bv = C_phi*C_theta*C_psi;
C_vb = transpose(C_bv);

C_mbeta = DCM(3, -deg2rad(beta_deg));
C_alpha = DCM(2, deg2rad(alpha_deg));

C_ba = C_alpha*C_mbeta;
C_ab = transpose(C_ba);

V_a = [ V ; 0 ; 0];
V_b = C_ba*V_a;

u = V_b(1); v = V_b(2); w = V_b(3);

p_rad_s = deg2rad(p_deg_s);
q_rad_s = deg2rad(q_deg_s);
r_rad_s = deg2rad(r_deg_s);

omega_b = [p_rad_s ; q_rad_s ; r_rad_s];

J_C_b = [ Ixx 0 -Ixz
          0 Iyy 0
          -Ixz 0 Izz];
g_v = [0 ; 0 ; g];
g_b = C_bv*g_v;


[F_aero_b, M_aero_C_b] = aero_loads(X,U,i_cond);
[F_prop_b,M_prop_C_b,T1, T2] = prop_loads(X,U);

V_b_dot = 1/m*(F_aero_b + F_prop_b) + g_b - skew(omega_b)*V_b;
omega_b_dot = J_C_b\(M_aero_C_b + M_prop_C_b - skew(omega_b)*J_C_b*omega_b);

udot = V_b_dot(1);
vdot = V_b_dot(2);
wdot = V_b_dot(3);

pdot_rad_s_2 = omega_b_dot(1);
qdot_rad_s_2 = omega_b_dot(2);
rdot_rad_s_2 = omega_b_dot(3);

pdot_deg_s_2 = rad2deg(pdot_rad_s_2);
qdot_deg_s_2 = rad2deg(qdot_rad_s_2);
rdot_deg_s_2 = rad2deg(rdot_rad_s_2);


Vdot = (u*udot + v*vdot + w*wdot)/V;
alphadot_rad_s = (u*wdot - w*udot)/(u^2 + w^2);
betadot_rad_s = (V*vdot - v*Vdot)/(V*sqrt(u^2 + w^2));

alphadot_deg_s = alphadot_rad_s*180/pi;
betadot_deg_s = betadot_rad_s*180/pi;

% Translational kinematics (Chapter 8, slides 32-33):

V_i = C_vb*V_b;
xdot = V_i(1); ydot = V_i(2); zdot = V_i(3);
hdot = -zdot;

% Rotational kinematics (Chapter 8, slides 18-20):
% As derived in class:
I_3 = eye(3);
e_31 = I_3(:,1);
e_32 = I_3(:,2);
e_33 = I_3(:,3);
K_Phi = [e_31 C_phi*e_32 C_bv*e_33];

% % Faster, alternative option:
% I_3 = eye(3);
% K_Phi = [I_3(:,1) C_phi(:,2) C_bv(:,3)];

Phidot_rad_s = K_Phi\omega_b;

phidot_rad_s = Phidot_rad_s(1);
thetadot_rad_s = Phidot_rad_s(2);
psidot_rad_s = Phidot_rad_s(3);

phidot_deg_s = phidot_rad_s*180/pi;
thetadot_deg_s = thetadot_rad_s*180/pi;
psidot_deg_s = psidot_rad_s*180/pi;

Xdot = [
    Vdot
    alphadot_deg_s
    qdot_deg_s_2
    thetadot_deg_s
    hdot
    xdot
    betadot_deg_s
    phidot_deg_s
    pdot_deg_s_2
    rdot_deg_s_2
    psidot_deg_s
    ydot
    ];

[CD,CY,CL,Cl,Cm,Cn] = aero_databank(X,U,i_cond);

[rho,~,~,a] = ISA(h);

Mach = V/a;
q_bar = 0.5*rho*V^2;

C_tv = C_ab*C_bv;
% Chapter 4, slide 43:
gamma_deg = -asind(C_tv(1,3));

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
