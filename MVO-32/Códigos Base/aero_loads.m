function [D, L, M, Y, l, n] = aero_loads(X,U)

global aircraft
run("create_aircraft.m");

V = X(1);
h = X(5);

[CD,CL,Cm,CY,Cl,Cn] = aero_databank(X,U);

rho = ISA(h);
q_bar = 0.5*rho*V^2;

S = aircraft.S;
c = aircraft.c;
b = aircraft.b;

D = q_bar*S*CD;
L = q_bar*S*CL;
M = q_bar*S*c*Cm;
Y = q_bar*S*CY;
l = q_bar*S*b*Cl;
n = q_bar*S*b*Cn;
