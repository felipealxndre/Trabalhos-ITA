function [D, Y, L, l, m, n] = aero_loads(X,U)

global aircraft

V = X(1);
h = X(5);

[CD, CY, CL, Cl, Cm, Cn] = aero_databank(X, U);

rho = ISA(h);
q_bar = 0.5*rho*V^2;

S = aircraft.S;
c = aircraft.c;
b = aircraft.b;

D = q_bar*S*CD;
L = q_bar*S*CL;
m = q_bar*S*c*Cm;
Y = q_bar*S*CY;
l = q_bar*S*b*Cl;
n = q_bar*S*b*Cn;
