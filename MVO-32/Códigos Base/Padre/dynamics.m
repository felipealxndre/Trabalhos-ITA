function [Xdot, Y]= dynamics(t, X, U)

global g
global aircraft

V         = X(1);
alpha_deg = X(2);
q_deg_s   = X(3);
theta_deg = X(4);
h         = X(5);
beta_deg  = X(7);
phi_deg   = X(8);
p_deg_s   = X(9);
r_deg_s   = X(10);
psi_deg   = X(11);

d2r = pi/180;
alpha = alpha_deg * d2r;
beta  = beta_deg  * d2r;
phi   = phi_deg   * d2r;
theta = theta_deg * d2r;
psi   = psi_deg   * d2r;

p = p_deg_s * d2r;
q = q_deg_s * d2r;
r = r_deg_s * d2r;

m   = aircraft.m;
Ixx = aircraft.Ixx;
Iyy = aircraft.Iyy;
Izz = aircraft.Izz;
Ixz = aircraft.Ixz;
S   = aircraft.S;
b   = aircraft.b;
c_bar = aircraft.c;

Gamma = Ixx*Izz - Ixz^2;
J_inv = [Izz/Gamma, 0, Ixz/Gamma; 0, 1/Iyy, 0; Ixz/Gamma, 0, Ixx/Gamma];

sa = sin(alpha); ca = cos(alpha);
sb = sin(beta);  cb = cos(beta);
st = sin(theta); ct = cos(theta);
sp = sin(phi);   cp = cos(phi);
ss = sin(psi);   cs = cos(psi);

[rho,~,~,a] = ISA(h);
q_bar = 0.5 * rho * V^2;

[CD, CL, Cm, CY, Cl, Cn] = aero_databank(X, U);

D = q_bar * S * CD;
Y_aero_wind = q_bar * S * CY;
L = q_bar * S * CL;

l_aero = q_bar * S * b * Cl;
m_aero = q_bar * S * c_bar * Cm;
n_aero = q_bar * S * b * Cn;

F_aero_a = [-D; -Y_aero_wind; -L];

C_b_a = [ ca*cb, -ca*sb, -sa;
            sb,     cb,    0;
          sa*cb, -sa*sb,  ca ];
          
F_aero_b = C_b_a * F_aero_a;
X_aero = F_aero_b(1);
Y_aero = F_aero_b(2);
Z_aero = F_aero_b(3);

F_M_prop = prop_loads(X, U);
X_prop = F_M_prop(1);
Y_prop = F_M_prop(2);
Z_prop = F_M_prop(3);
l_prop = F_M_prop(4);
m_prop = F_M_prop(5);
n_prop = F_M_prop(6);

mg = m * g;
X_grav = -mg * st;
Y_grav =  mg * sp * ct;
Z_grav =  mg * cp * ct;

X_total = X_aero + X_prop + X_grav;
Y_total = Y_aero + Y_prop + Y_grav;
Z_total = Z_aero + Z_prop + Z_grav;

L_total = l_aero + l_prop;
M_total = m_aero + m_prop;
N_total = n_aero + n_prop;

u = V * ca * cb;
v = V * sb;
w = V * sa * cb;

u_dot = X_total/m - q*w + r*v;
v_dot = Y_total/m - r*u + p*w;
w_dot = Z_total/m - p*v + q*u;

V_dot = (u*u_dot + v*v_dot + w*w_dot) / V;

if (u^2 + w^2) == 0
    alpha_dot = 0;
else
    alpha_dot = (u*w_dot - w*u_dot) / (u^2 + w^2);
end

if (V * sqrt(u^2 + w^2)) == 0
    beta_dot = 0;
else
    beta_dot  = (V*v_dot - v*V_dot) / (V * sqrt(u^2 + w^2));
end

J = [Ixx, 0, -Ixz; 0, Iyy, 0; -Ixz, 0, Izz];
omega_b = [p; q; r];
omega_tilde = [0, -r, q; r, 0, -p; -q, p, 0];

M_total_vec = [L_total; M_total; N_total];

omega_dot = J_inv * (M_total_vec - omega_tilde * J * omega_b);

p_dot = omega_dot(1);
q_dot = omega_dot(2);
r_dot = omega_dot(3);

x_dot = u*(ct*cs) + v*(sp*st*cs - cp*ss) + w*(cp*st*cs + sp*ss);
y_dot = u*(ct*ss) + v*(sp*st*ss + cp*cs) + w*(cp*st*ss - sp*cs);
z_dot = -u*st     + v*(sp*ct)            + w*(cp*ct);
h_dot = -z_dot;

phi_dot   = p + tan(theta)*(q*sp + r*cp);
theta_dot = q*cp - r*sp;
psi_dot   = (q*sp + r*cp) / ct;

Xdot = [
    V_dot
    rad2deg(alpha_dot)
    rad2deg(q_dot)
    rad2deg(theta_dot)
    h_dot
    x_dot
    rad2deg(beta_dot)
    rad2deg(phi_dot)
    rad2deg(p_dot)
    rad2deg(r_dot)
    rad2deg(psi_dot)
    y_dot
    ];

delta_1 = U(1);
delta_2 = U(2);
Tmax  = aircraft.Tmax;
n_rho = aircraft.n_rho;
rho0  = 1.225;
T1 = delta_1 * Tmax * (rho / rho0)^n_rho;
T2 = delta_2 * Tmax * (rho / rho0)^n_rho;

gamma_deg = theta_deg - alpha_deg;
Mach = V/a;

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
    p_dot
    q_dot
    r_dot
    ];
end