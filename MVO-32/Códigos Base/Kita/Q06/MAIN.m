clear
close all
clc

% MVO-32 ESTABILIDADE E CONTROLE DE AERONAVES
% Dupla: Thiago Akira Missato e Eduardo Pinto Ferreira.
% Código MAIN.m referente a Questão 06 da Lista 03.
% Elaborado por Eduardo Pinto Ferreira.


global aircraft

create_aircraft; 

V_cruise = 228.68;
h_cruise = 11000;
[rho, ~, ~, ~] = ISA(h_cruise);

Ixx = aircraft.Ixx;
S = aircraft.S;
b = aircraft.b;
Cl_p = aircraft.Cl_p;

% Modelo
lambda= (rho * V_cruise * S * b^2 * Cl_p) / (4 * Ixx);

tau_r = -1 / lambda;

fprintf('Autovalor: %.4f (1/s)\n', lambda);
fprintf('tau_r:     %.4f (s)\n', tau_r);