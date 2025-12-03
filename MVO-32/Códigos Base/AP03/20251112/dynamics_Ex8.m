function [Xdot,Y] = dynamics_Ex8(t,X,U,i_cond)

global g
global aircraft

V = X(1);
alpha_deg = X(2);
q_deg_s = X(3);
theta_deg = X(4);
h = X(5);
%x = X(6);

beta_deg = X(7);
phi_deg = X(8);
p_deg_s = X(9);
r_deg_s = X(10);
psi_deg = X(11);
%y = X(12);

% throttle = U(1);
% i_t_deg = U(2);
% delta_e_deg = U(3);
% delta_a_deg = U(4);
% delta_r_deg = U(5);

m = aircraft.m;

Ixx = aircraft.Ixx;
Iyy = aircraft.Iyy;
Izz = aircraft.Izz;
Ixz = aircraft.Ixz;


C_psi = DCM(3, psi_deg * pi/180);
C_theta = DCM(2, theta_deg * pi/180);
C_phi = DCM(1, phi_deg * pi/180);

%usar o apostrofe pode dar problema quando a matriz e complexa, pois ele
%transpoe e conjuga. Usa A.' ou transpose(A)

C_bv = C_phi*C_theta*C_psi;
C_vb = C_bv.';

C_mbeta = DCM(3, -beta_deg*pi/180); %rotacao de -beta
C_alpha = DCM(2, alpha_deg*pi/180);

C_ba = C_alpha*C_mbeta;
C_ab = C_ba.';

V_a = [V; 0; 0];
V_b = C_ba * V_a;

u = V_b(1);
v = V_b(2);
w = V_b(3);

p_rad_s = p_deg_s*pi/180;
q_rad_s = q_deg_s*pi/180;
r_rad_s = r_deg_s*pi/180;

omega_b = [p_rad_s; q_rad_s; r_rad_s];

%matriz de inercia em torno do centro de massa no eixo do corpo
J_C_b = [
    Ixx 0 -Ixz
    0 Iyy 0
    -Ixz 0 Izz];

g_v = [0; 0; g];
g_b = C_bv*g_v;

[F_aero_b, M_aero_C_b] = aero_loads_Ex8(X, U, i_cond);
[F_prop_b, M_prop_C_b, T] = prop_loads(X,U);

V_b_dot = 1/m * (F_aero_b + F_prop_b) + g_b - skew(omega_b)*V_b;
omega_b_dot = J_C_b \ (M_aero_C_b + M_prop_C_b - skew(omega_b)*J_C_b*omega_b);

udot = V_b_dot(1);
vdot = V_b_dot(2);
wdot = V_b_dot(3);

Vdot = (u*udot + v*vdot + w*wdot)/V;
alphadot_rad_s = (u*wdot - w*udot)/(u^2 + w^2);
betadot_rad_s = (V*vdot - v*Vdot)/(V*sqrt(u^2 + w^2));


pdot_rad_s_2 = omega_b_dot(1);
qdot_rad_s_2 = omega_b_dot(2);
rdot_rad_s_2 = omega_b_dot(3);

pdot_deg_s_2 = pdot_rad_s_2*180/pi;
qdot_deg_s_2 = qdot_rad_s_2*180/pi;
rdot_deg_s_2 = rdot_rad_s_2*180/pi;

alphadot_deg_s = alphadot_rad_s*180/pi;
betadot_deg_s = betadot_rad_s*180/pi;

% Translational kinematics (Chapter 8, slides 32-33):
%x,y,z,h coordenadas do cg no sistema terrestre fixo
V_i = C_vb * V_b;
xdot = V_i(1);
ydot = V_i(2);
zdot = V_i(3);
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

[CD,CY,CL,Cl,Cm,Cn] = aero_databank_Ex8(X,U,i_cond);

[rho,~,~,a] = ISA(h);

Mach = V/a;
q_bar = 0.5*rho*V^2;

C_tv = C_ab*C_bv;
% Chapter 4, slide 43:
gamma_deg = -asind(C_tv(1,3));

Y = [
    gamma_deg
    T
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
