function [D, L, M] = aero_loads(X,U)

global aircraft

V = X(1);
h = X(5);

[CD, CL, Cm] = aero_databank(X, U);

rho = ISA(h);
q_bar = 0.5*rho*V^2;

S = aircraft.S;
c = aircraft.c;

D = q_bar*S*CD;
L = q_bar*S*CL;
M = q_bar*S*c*Cm;