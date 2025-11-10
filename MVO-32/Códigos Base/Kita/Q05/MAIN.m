clear
close all
clc

% MVO-32 ESTABILIDADE E CONTROLE DE AERONAVES
% Dupla: Thiago Akira Missato e Eduardo Pinto Ferreira.
% Código MAIN.m referente a Questão 05 da Lista 03.
% Elaborado por Thiago Akira Missato.

global g
global aircraft

g = 9.80665;
create_aircraft;

V_cruise = 228.68;
h_cruise = 11000;
psidot_deg_s_eq = 0; % Voo reto

trim_par(1) = struct('V',V_cruise,'h',h_cruise,'gamma_deg',0,... % Voo reto nivelado
    'thetadot_deg_s',0,'psidot_deg_s',psidot_deg_s_eq);

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

trim_output(1:3) = struct('X_eq',zeros(12,1),'U_eq',zeros(6,1),...
    'Y_eq',zeros(12,1));

%--------------------------------------------------------------------------
% Ex. 2) Condição de cruzeiro reto e nivelado.

i_cond = 1; % Condição de cruzeiro

x_eq_0 = zeros(12,1);
x_eq_0(1) = trim_par(i_cond).V;
x_eq_0(3) = 0;
x_eq = fsolve(@trim_function,x_eq_0,options,trim_par,i_cond);
[~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par,i_cond);

trim_output(i_cond).X_eq = X_eq;
trim_output(i_cond).U_eq = U_eq;
trim_output(i_cond).Y_eq = Y_eq;

%--------------------------------------------------------------------------
% Ex. 5) Simulação de doublet de empuxo diferencial.

X0 = trim_output(i_cond).X_eq;
U0 = trim_output(i_cond).U_eq;

dt = 0.010;
t_f = 30;
T = 0:dt:t_f;
[X,Y] = ode4xy(@simulate_flight_Ex5,T,X0,U0, i_cond); 

U = Y(:,end-5:end);
Y = Y(:,1:end-6);

%plot_long
plot_latdir
%plot_controls
%plot_path
%plot_outputs
