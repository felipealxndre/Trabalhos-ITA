clear
close all
clc

global g
global aircraft

run("create_aircraft.m");

g = 9.80665;

trim_par = struct('V',228.8138886191133,'h',11000,'gamma_deg',-3,...
    'thetadot_deg_s',0);

% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);
% % In newer MATLAB versions, the following command might be necessary:
% options = optimoptions(@fsolve,'Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10);
% % or:
% options = optimoptions(@fsolve,'Display','iter','TolX',1e-10,'TolFun',1e-10);

%% Ex.2

trim_output = struct('X_eq',zeros(12,1),'U_eq',zeros(6,1),'Y_eq',zeros(12,1));

x_eq_0 = zeros(6,1);
x_eq_0(1) = trim_par.V;
x_eq_0(1) = trim_par.V; 
x_eq_0(2) = 3.0;         
x_eq_0(3) = 0.0;          
x_eq_0(4) = x_eq_0(2);   
x_eq_0(5) = 0.50;        
x_eq_0(6) = 0.0;           
x_eq = fsolve(@trim_function,x_eq_0,options,trim_par);
[~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par);

trim_output.X_eq = X_eq;
trim_output.U_eq = U_eq;
trim_output.Y_eq = Y_eq;

% Impressão resumida
fprintf('\n----- EX.2 -----\n\n');
fprintf('   %-12s = %8.2f %s\n','V',X_eq(1),'m/s');
fprintf('   %-12s = %8.3f %s\n','alpha',X_eq(2),'deg');
fprintf('   %-12s = %8.3f %s\n','q',X_eq(3),'deg/s');
fprintf('   %-12s = %8.3f %s\n','theta',X_eq(4),'deg');
fprintf('   %-12s = %8.1f %s\n','h',X_eq(5),'m');
fprintf('   %-12s = %8.3f %s\n','beta',X_eq(7),'deg');
fprintf('   %-12s = %8.3f %s\n','phi',X_eq(8),'deg');
fprintf('   %-12s = %8.3f %s\n','p',X_eq(9),'deg/s');
fprintf('   %-12s = %8.3f %s\n','r',X_eq(10),'deg/s');
fprintf('\n');
fprintf('   %-12s = %8.2f %s\n','delta1',U_eq(1),'%');
fprintf('   %-12s = %8.2f %s\n','delta2',U_eq(2),'%');
fprintf('   %-12s = %8.3f %s\n','i_t',U_eq(3),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_e',U_eq(4),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_a',U_eq(5),'deg');
fprintf('   %-12s = %8.3f %s\n','delta_r',U_eq(6),'deg');
fprintf('\n');
fprintf('   %-12s = %8.3f %s\n','gamma',Y_eq(1),'deg');
fprintf('   %-12s = %8.1f %s\n','T1',Y_eq(2),'N');
fprintf('   %-12s = %8.1f %s\n','T2',Y_eq(3),'N');
fprintf('   %-12s = %8.3f %s\n','Mach',Y_eq(4),'');
fprintf('   %-12s = %8.4f %s\n','C_D',Y_eq(5),'');
fprintf('   %-12s = %8.4f %s\n','C_L',Y_eq(6),'');
fprintf('   %-12s = %8.4f %s\n','C_m',Y_eq(7),'');
fprintf('   %-12s = %8.4f %s\n','C_Y',Y_eq(8),'');
fprintf('   %-12s = %8.4f %s\n','C_l',Y_eq(9),'');
fprintf('   %-12s = %8.4f %s\n','C_n',Y_eq(10),'');
fprintf('   %-12s = %8.4f %s\n','rho',Y_eq(11),'kg/m^3');
fprintf('   %-12s = %8.1f %s\n','qbar',Y_eq(12),'N/m^2');
fprintf('\n');

save trim_output_L3_Ex2.mat trim_output

%% Ex.3

%% Ex.4

%% Ex.5

%% Ex.6

%% Ex.7