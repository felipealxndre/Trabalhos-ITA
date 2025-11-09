function T = prop_loads(X, U)

global aircraft

h = X(5);
rho = ISA(h);

delta_1 = U(1);
delta_2 = U(2);

Tmax  = aircraft.Tmax;
n_rho = aircraft.n_rho;
rho0  = 1.225;

i_1   = aircraft.iota_1_deg;
tau_1 = aircraft.tau_1_deg;
r1    = [aircraft.x_1; aircraft.y_1; aircraft.z_1];

i_2   = aircraft.iota_2_deg;
tau_2 = aircraft.tau_2_deg;
r2    = [aircraft.x_2; aircraft.y_2; aircraft.z_2];

d2r = pi/180;
i_1_rad   = i_1 * d2r;
tau_1_rad = tau_1 * d2r;
i_2_rad   = i_2 * d2r;
tau_2_rad = tau_2 * d2r;

T1 = delta_1 * Tmax * (rho / rho0)^n_rho;
T2 = delta_2 * Tmax * (rho / rho0)^n_rho;

ci1 = cos(i_1_rad);
si1 = sin(i_1_rad);
ct1 = cos(tau_1_rad);
st1 = sin(tau_1_rad);

X_prop_1 = T1 * ci1 * ct1;
Y_prop_1 = T1 * ci1 * st1;
Z_prop_1 = -T1 * si1;
F_prop_1 = [X_prop_1; Y_prop_1; Z_prop_1];

M_prop_1 = cross(r1, F_prop_1);

ci2 = cos(i_2_rad);
si2 = sin(i_2_rad);
ct2 = cos(tau_2_rad);
st2 = sin(tau_2_rad);

X_prop_2 = T2 * ci2 * ct2;
Y_prop_2 = T2 * ci2 * st2;
Z_prop_2 = -T2 * si2;
F_prop_2 = [X_prop_2; Y_prop_2; Z_prop_2];

M_prop_2 = cross(r2, F_prop_2);

F_total = F_prop_1 + F_prop_2;
M_total = M_prop_1 + M_prop_2;

T = [F_total; M_total];