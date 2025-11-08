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
