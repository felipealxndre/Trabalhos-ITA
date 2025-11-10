clear
close all
clc

% MVO-32 ESTABILIDADE E CONTROLE DE AERONAVES
% Dupla: Thiago Akira Missato e Eduardo Pinto Ferreira.
% Código MAIN.m referente a Questão 02 da Lista 03.
% Elaborado em conjunto por Thiago Akira Missato e Eduardo Pinto Ferreira.


global g
global aircraft

g = 9.80665;

create_aircraft;
i_cond = 1;

psidot_deg_s_eq = 0;

Mach = 0.775; % Cruise Mach Number
h = 11000; % Cruise altitude [m]
[~,~,~,a] = ISA(h);

trim_par(1) = struct('V', Mach * a,'h',11000,'gamma_deg',0,...
    'thetadot_deg_s',0,'psidot_deg_s',psidot_deg_s_eq);

options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);

%--------------------------------------------------------------------------
% Ex. 2) Trimmed flight conditions:

trim_output(1:3) = struct('X_eq',zeros(12,1),'U_eq',zeros(6,1),...
    'Y_eq',zeros(12,1));

for i_cond=1:1
    x_eq_0 = zeros(12,1);
    x_eq_0(1) = trim_par(i_cond).V;
    x_eq_0(3) = 0;
    x_eq = fsolve(@trim_function,x_eq_0,options,trim_par,i_cond);
    [~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par,i_cond);

    trim_output(i_cond).X_eq = X_eq;
    trim_output(i_cond).U_eq = U_eq;
    trim_output(i_cond).Y_eq = Y_eq;

    fprintf('----- Cruise Flight Condition -----\n\n');
    fprintf('   %-10s = %10.4f %-4s\n','gamma',trim_par(i_cond).gamma_deg,'deg');
    fprintf('   %-10s = %10.4f %-4s\n','theta_dot',trim_par(i_cond).thetadot_deg_s,'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','psi_dot',trim_par(i_cond).psidot_deg_s,'deg/s');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','V',X_eq(1),'m/s');
    fprintf('   %-10s = %10.4f %-4s\n','alpha',X_eq(2),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','q',X_eq(3),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','theta',X_eq(4),'deg');
    fprintf('   %-10s = %10.1f %-4s\n','h',X_eq(5),'m');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','beta',X_eq(7),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','phi',X_eq(8),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','p',X_eq(9),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','r',X_eq(10),'deg/s');
    fprintf('   %-10s = %10.4f %-4s\n','psi',X_eq(11),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.2f %-4s\n','throttle_1',U_eq(1)*100,'%');
    fprintf('   %-10s = %10.2f %-4s\n','throttle_2',U_eq(1)*100,'%');
    fprintf('   %-10s = %10.2f %-4s\n','Thrust_1',Y_eq(2),'N');
    fprintf('   %-10s = %10.2f %-4s\n','Thrust_2',Y_eq(3),'N');
    fprintf('   %-10s = %10.4f %-4s\n','i_t',U_eq(3),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_e',U_eq(4),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_a',U_eq(5),'deg');
    fprintf('   %-10s = %10.4f %-4s\n','delta_r',U_eq(6),'deg');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','Mach',Y_eq(4),'');
    fprintf('   %-10s = %10.4f %-4s\n','rho',Y_eq(11),'kg/m³');
    fprintf('   %-10s = %10.2f %-4s\n','q bar',Y_eq(12),'N/m²');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','C_D',Y_eq(5),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_L',Y_eq(6),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_m',Y_eq(7),'');
    fprintf('\n');
    fprintf('   %-10s = %10.4f %-4s\n','C_Y',Y_eq(8),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_l',Y_eq(9),'');
    fprintf('   %-10s = %10.4f %-4s\n','C_n',Y_eq(10),'');
    fprintf('\n');
end

save trim_output.mat trim_output
