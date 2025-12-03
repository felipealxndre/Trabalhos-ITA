function [F_prop_b,M_prop_C_b,T] = prop_loads(X,U)

global aircraft

i_p_rad = aircraft.i_p_deg*pi/180;
n_rho = aircraft.n_rho;
Tmax = aircraft.Tmax;
x_p = aircraft.x_p;
z_p = aircraft.z_p;

h = X(5);
rho = ISA(h);

throttle = U(1);
T = throttle*Tmax*(rho/1.225)^n_rho;

C_iota_p = DCM(2,i_p_rad);
C_pb = C_iota_p;

F_prop_b = C_pb.'*[T; 0; 0];

M_prop_C_b = skew([x_p; 0; z_p])*F_prop_b;
