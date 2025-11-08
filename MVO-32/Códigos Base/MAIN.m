clear
close all
clc

global g
global aircraft

run("create_aircraft.m");

g = 9.80665;

trim_par = struct('V',77,'h',600,'gamma_deg',-3,...
    'thetadot_deg_s',0);

% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10);
% % In newer MATLAB versions, the following command might be necessary:
% options = optimoptions(@fsolve,'Display','iter','StepTolerance',1e-10,'FunctionTolerance',1e-10);
% % or:
% options = optimoptions(@fsolve,'Display','iter','TolX',1e-10,'TolFun',1e-10);

%% Ex.2

trim_output = struct('X_eq',zeros(6,1),'U_eq',zeros(3,1),...
    'Y_eq',zeros(8,1));

x_eq_0 = zeros(6,1);
x_eq_0(1) = trim_par.V;
x_eq = fsolve(@trim_function,x_eq_0,options,trim_par);
[~,X_eq,U_eq,Y_eq] = trim_function(x_eq,trim_par);

trim_output.X_eq = X_eq;
trim_output.U_eq = U_eq;
trim_output.Y_eq = Y_eq;

%% Ex.3

%% Ex.4

%% Ex.5

%% Ex.6

%% Ex.7